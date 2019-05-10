/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2017 James Cussens, Mark Bartlett        *
 *                                                                       *
 *   This program is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU General Public License as      *
 *   published by the Free Software Foundation; either version 3 of the  *
 *   License, or (at your option) any later version.                     *
 *                                                                       *
 *   This program is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   General Public License for more details.                            *
 *                                                                       *
 *   You should have received a copy of the GNU General Public License   *
 *   along with this program; if not, see                                *
 *   <http://www.gnu.org/licenses>.                                      *
 *                                                                       *
 *   Additional permission under GNU GPL version 3 section 7             *
 *                                                                       *
 *   If you modify this Program, or any covered work, by linking or      *
 *   combining it with SCIP (or a modified version of that library),     *
 *   containing parts covered by the terms of the ZIB Academic License,  *
 *   the licensors of this Program grant you additional permission to    *
 *   convey the resulting work.                                          *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** @file
 *  Implements the functions needed to perform model averaging over the n best Bayesian networks.
 */

#include "model_averaging.h"
#include "parent_set_data.h"

/** Find the index of @link ma_info->vars @endlink and @link ma_info->average_scores @endlink
 *  relating to a variable.
 *
 *  @param var The variable of interest.
 *  @param ma_info Model averaging information.
 *  @return The index of the var in @link ma_info->vars @endlink.  If the variable is not
 *  included in the averaging, -1 is returned.
 */
static int indexOf(
   MA_info* ma_info,
   SCIP_VAR* var
   )
{
   int i;
   for( i = 0; i < ma_info->num_vars; i++ )
      if( ma_info->vars[i] == var )
         return i;
   return -1;
}


/** Adds parameters for controlling the model averaging.
 *  @param scip The SCIP instance to which the parameter is to be added.
 *  @param ma_info Model averaging information.
 *  @return SCIP_OKAY if the parameters were added successfully or an error code otherwise.
 */
SCIP_RETCODE MA_addAveragingParameters(
   SCIP* scip,
   MA_info* ma_info
   )
{
   SCIP_CALL(SCIPaddBoolParam(scip, "gobnilp/logaverage", "whether model averaging should assume logarithmic objective function",
         &(ma_info->is_log_score), FALSE, TRUE, NULL, NULL));
   return SCIP_OKAY;
}

/** Allocates memory for the data structures used for model averaging.
 *
 *  @param scip The SCIP instance on which the model averaging will be performed.
 *  @param ma_info Model averaging information.
 *  @return SCIP_OKAY if memory allocation was successful or an appropriate error
*   message otherwise.
 */
SCIP_RETCODE MA_createAverageDataStructure(
   SCIP* scip,
   MA_info* ma_info
   )
{
   int i;
   int j;
   int num_all_vars;
   SCIP_VAR** all_vars;

   assert( scip != NULL );
   assert( ma_info != NULL );

   ma_info->total_score = 0.0;
   ma_info->total_time = 0.0;
   ma_info->is_first_score = TRUE;
   ma_info->first_score = 0.0;
   ma_info->normalising_constant = 0.0;
   
   ma_info->num_vars = SCIPgetNBinVars(scip);
   num_all_vars = SCIPgetNVars(scip);
   all_vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(ma_info->vars), ma_info->num_vars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ma_info->average_scores), ma_info->num_vars) );

   /* Copy the binary variables */
   j = 0;
   for( i = 0; i < num_all_vars; i++ )
      if( SCIPvarIsBinary(all_vars[i]) )
      {
         ma_info->vars[j] = all_vars[i];
         j++;
      }
   assert(j == ma_info->num_vars);

   /* Set the averages to initial values of 0 */
   for( i = 0; i < ma_info->num_vars; i++ )
      ma_info->average_scores[i] = 0.0;

   return SCIP_OKAY;
}
/** Frees memory used for the data structures used for model averaging.
 *
 *  @param scip The SCIP instance on which the model averaging was performed.
 *  @param ma_info Model averaging information.
 *  @return SCIP_OKAY if memory deallocation was successful or an appropriate error
 *  message otherwise.
 */
SCIP_RETCODE MA_destroyAverageDataStructure(
   SCIP* scip,
   MA_info* ma_info
   )
{
   assert( scip != NULL );
   assert( ma_info != NULL );

   SCIPfreeMemoryArray(scip, &(ma_info->vars));
   SCIPfreeMemoryArray(scip, &(ma_info->average_scores));
   return SCIP_OKAY;
}

/** Updates the average scores based on a newly found solution.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param ma_info Model averaging information.
 *  @param sol The new solution to incorporate in to the averages.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error message otherwise.
 */
SCIP_RETCODE MA_updateAverageDataStructure(
   SCIP* scip,
   MA_info* ma_info,
   SCIP_SOL* sol
   )
{
   int i;
   SCIP_Real overall_score = 0;
   SCIP_Real exp_score = 0;

   assert( scip != NULL );
   assert( ma_info != NULL );
   assert( sol != NULL );

   /* Calculate the total solution value */
   for( i = 0; i < ma_info->num_vars; i++ )
      overall_score += SCIPgetSolVal(scip, sol, (ma_info->vars)[i]) * SCIPvarGetObj((ma_info->vars)[i]);
   if( ma_info->is_first_score )
   {
      ma_info->first_score = overall_score;
      ma_info->is_first_score = FALSE;
   }
   exp_score = exp(overall_score - ma_info->first_score);

   /* Update the averages */
   if( ma_info->is_log_score )
   {
      for( i = 0; i < ma_info->num_vars; i++ )
         if( SCIPgetSolVal(scip, sol, ma_info->vars[i]) > 0.5 )
            ma_info->average_scores[i] += exp_score;
   }
   else
   {
      for( i = 0; i < ma_info->num_vars; i++ )
         if( SCIPgetSolVal(scip, sol, ma_info->vars[i]) > 0.5 )
            ma_info->average_scores[i] += overall_score;
   }

   /* Update the totals */
   ma_info->total_score += overall_score;
   ma_info->total_time += SCIPgetSolvingTime(scip);
   ma_info->normalising_constant += exp_score;

   return SCIP_OKAY;
}

/** Returns the model average value of a given variable.
 *
 *  If the variable is not part of this model averaging, the value -1 will be returned.
 *
 *  @param ma_info Model averaging information.
 *  @param variable The variable to get the model averaging score for.
 *  @return The model average score of the variable.
 */
SCIP_Real MA_getAverageValue(
   MA_info* ma_info,
   SCIP_VAR* variable
   )
{
   int index = indexOf(ma_info, variable);
   if( index == -1 )
      return -1;
   else if( ma_info->is_log_score )
      return ma_info->average_scores[index] / ma_info->normalising_constant;
   else
      return ma_info->average_scores[index] / ma_info->total_score;
}
/** Returns the total time spent solving for all the solutions included
 *  in the model average.
 *
 *  @param ma_info Model averaging information.
 *  @return The number of seconds spend on solving all of the solutions
 *  that have been used to make up the average.
 */
SCIP_Real MA_getTotalAveragesTime(
   MA_info* ma_info
   )
{
   return ma_info->total_time;
}

/** Returns the total likelihood of all the solutions included
 *  in the model average.
 *
 *  @param ma_info Model averaging information.
 *  @return The sum of the likelihoods of each of the solutions included.
 */
SCIP_Real MA_getTotalAveragesScore(
   MA_info* ma_info
   )
{
   return ma_info->total_score;
}
