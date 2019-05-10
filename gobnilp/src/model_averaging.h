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
 *  Function declarations for model_averaging.c
 */

#ifndef __MODEL_AVERAGING_H__
#define __MODEL_AVERAGING_H__

#include <scip/scip.h>

typedef struct
{
   int num_vars;                   /**< The number of variables that this modelling averaging is performed on. */
   SCIP_VAR** vars;                /**< The variables which the model averaging applies to. */
   SCIP_Real* average_scores;      /**< The average scores of each of the variables in the program.
                                    * 
                                    *  average_scores[i] is the likelihood weighted average score of vars[i].
                                    */ 
   SCIP_Real total_score;          /**< The sum of the likelihoods of all networks found so far.
                                    *
                                    *  This is used to normalise the scores in @link average_scores @endlink .
                                    */
   SCIP_Real total_time;           /**< The number of seconds spent finding all of the solutions included in the average. */
   SCIP_Bool is_log_score;         /**< Whether the objective function is logarithmic.*/
   SCIP_Bool is_first_score;       /**< Whether the next solution is the first solution found. */
   SCIP_Real first_score;          /**< The score of the first solution found. */
   SCIP_Real normalising_constant; /**< The normalising constant for logarithmic scores. */
} MA_info;

extern SCIP_RETCODE MA_addAveragingParameters(SCIP* scip, MA_info* ma_info);

extern SCIP_RETCODE MA_createAverageDataStructure(SCIP* scip, MA_info* ma_info);
extern SCIP_RETCODE MA_destroyAverageDataStructure(SCIP* scip, MA_info* ma_info);

extern SCIP_RETCODE MA_updateAverageDataStructure(SCIP* scip, MA_info* ma_info, SCIP_SOL* sol);

extern    SCIP_Real MA_getAverageValue(MA_info* ma_info, SCIP_VAR* variable);
extern    SCIP_Real MA_getTotalAveragesTime(MA_info* ma_info);
extern    SCIP_Real MA_getTotalAveragesScore(MA_info* ma_info);

#endif
