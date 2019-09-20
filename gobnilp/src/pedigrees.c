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
 *  Contains all the functions related to the use of pedigrees in the program.
 */

#include "probdata_bn.h"
#include "pedigrees.h"
#include "pedigree_data.h"
#include "parent_set_data.h"
#include "metadata.h"
#include "utils.h"
#include "versiongit.h"
#include <string.h>

/* Useful auxillary functions */

/** @brief Determines whether the program is being run for learning pedigrees.
 *
 *  @param scip The SCIP instance that is running.
 *  @return TRUE if the program is being used to learn pedigrees. FALSE otherwsie.
 */
SCIP_Bool PD_inPedigreeMode(
   SCIP* scip
   )
{
   SCIP_Bool pedigreemode;
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/pedigreemode", &pedigreemode) );
   return pedigreemode;
}
/** Checks whether the program is enforcing sex consistency.
 *
 *  If the deprecated version of rhe sex consistency parameter is being used then
 *  a message will warn the user of this.
 *
 *  @param scip The SCIP instance of which to check the status.
 *  @return TRUE if the program should enforcing sex consistency.  FALSE otherwise.
 */
static
SCIP_Bool usingSexConsistency(
   SCIP* scip
   )
{
   if( PD_inPedigreeMode(scip) )
   {
      SCIP_Bool sexconsistentold;
      SCIP_Bool sexconsistentnew;
      SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/sexconsistent", &sexconsistentold) );
      SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/pedigree/sexconsistent", &sexconsistentnew) );
      if( sexconsistentnew )
      {
         return TRUE;
      }
      else if( !sexconsistentnew && sexconsistentold )
      {
         /* Old sex consistency parameter is set to true and new version is either unused or set to false. */
         printf("WARNING: gobnilp/sexconsistent has been deprecated.\n Please use gobnilp/pedigree/sexconsistent instead.\n");
         return TRUE;
      }
      else
      {
         return FALSE;
      }
   }
   else
   {
      return FALSE;
   }
}
/** Checks whether a SCIP problem is suitable for use as a pedigree.
 *
 *  A problem is acceptable as a pedigree if it has no parent sets of more than size 2.
 *
 *  @param scip The SCIP instance for which to check the data.
 *  @param psd The parent set data associated with this problem.
 *  @return TRUE if this problem represents a valid pedigree problem, FALSE otherwise.
 */
static
SCIP_Bool checkSuitableForPedigree(
   SCIP* scip,
   ParentSetData* psd
   )
{
   int i, j;
   if( psd != NULL )
      for( i = 0; i < psd->n; i++ )
         for( j = 0; j < psd->nParentSets[i]; j++ )
            if( psd->nParents[i][j] > 2 )
               return FALSE;
   return TRUE;
}

/* Functions related to parameters */
/** Makes SCIP recognise parameters related to pedigree reconstruction.
 *
 *  @param scip The SCIP instance to which to add the parameters.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error code otherwise.
 */
SCIP_RETCODE PD_addPedigreeParameters(
   SCIP* scip
   )
{
   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/pedigreemode",
                             "whether to use GOBNILP for pedigrees",
                             FALSE
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/sexconsistent",
                             "whether to enforce sexual consistency in the dag",
                             FALSE
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/pedigree/sexconsistent",
                             "whether to enforce sexual consistency in the dag",
                             FALSE
                            ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxparentchildagegap",
                            "maximum age gap permitted between parent and child (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   
   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxsibagegap",
                            "maximum age gap permitted between full siblings (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxhalfsibagegap",
                            "maximum age gap permitted between half siblings (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxsibsetsize",
                            "maximum number of children a pair of parents can have together (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxchildren",
                            "maximum number of children any individual can have (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxchildrenmother",
                            "maximum number of children a mother can have (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/maxchildrenfather",
                            "maximum number of children a father can have (-1 for no restriction)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/pedigree/singleparents",
                             "whether the pedigree can feature individuals with only one parent",
                             TRUE
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/pedigree/inputfile/columnheaders",
                             "whether the pedigree input file has column headers",
                             TRUE
                            ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/inputfile/namecolumn",
                            "the column with individuals' names (-1 if no such column)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/inputfile/sexcolumn",
                            "the column with individuals' sexes (-1 if no such column)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/inputfile/agecolumn",
                            "the column with individuals' ages (-1 if no such column)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/inputfile/allelestartcolumn",
                            "the first column with individuals' alleles",
                            0, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/inputfile/alleleendcolumn",
                            "the last column with individuals' alleles",
                            1, 1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/inputfile/adultcolumn",
                            "the column stating which individuals are adults (-1 if no such column)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addStringParam(scip,
                            "gobnilp/pedigree/realdata/unknown",
                            "the symbol used for unknown genotypes",
                            "NA"
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/pedigree/realdata/maxmismatches",
                            "the most genotype mismatches that are allowed between a child and parent",
                            0, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addRealParam(scip,
                            "gobnilp/pedigree/realdata/mutationprob",
                            "the value to use for mutation probability",
                            0, 0, 1
                           ));

   SCIP_CALL( UT_addStringParam(scip,
         "gobnilp/pedigree/constraintsfile",
         "file containing constraints on pedigree structure",
         ""
         ) );

   return SCIP_OKAY;
}

/* Functions related to variables */
/** Creates any needed variables needed specifically for learning pedigrees.
 *
 *  @param scip The SCIP instance in which to create the variables.
 *  @param psd The parent set data associated with this problem.
 *  @param peddata The pedigree data to store the variables in.
 *  @return SCIP_OKAY if the variable was successful, or an approriate error code otherwise.
 */
static
SCIP_RETCODE PD_addPedigreeVariables(
   SCIP* scip,
   ParentSetData* psd,
   PedigreeData* peddata
   )
{
   if( !checkSuitableForPedigree(scip, psd) )
   {
      SCIPerrorMessage("The specified problem is not a valid pedigree problem.\n");
      return SCIP_ERROR;
   }
   else
   {
      if( usingSexConsistency(scip) )
      {
         if( peddata->SexVars == NULL )
         {
            int i;
            char s[SCIP_MAXSTRLEN];
            SCIP_PROBDATA* probdata;

            /* get problem data */
            probdata = SCIPgetProbData(scip);
            assert( probdata != NULL );

            SCIP_CALL( SCIPallocMemoryArray(scip, &(peddata->SexVars), peddata->n) );
            SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->isFemale), peddata->n) );
            for( i = 0; i < peddata->n; i++ )
            {
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "isFemale(%s)", psd->nodeNames[i]);
               SCIP_CALL( SCIPcreateVar(scip, &(peddata->SexVars[i]), s, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, peddata->SexVars[i]) );
               /* also save var in probdata for eventual release */
               probdata->isFemale[i] =  peddata->SexVars[i];
            }
         }
         if( peddata->sexes != NULL )
         {
            int i;
            SCIP_Bool infeasible, fixed;
            for( i = 0; i < peddata->n; i++ )
            {
               if( peddata->sexes[i] == 'M' )
                  SCIP_CALL( SCIPfixVar(scip, peddata->SexVars[i], 0, &infeasible, &fixed) );
               if( peddata->sexes[i] == 'F' )
                  SCIP_CALL( SCIPfixVar(scip, peddata->SexVars[i], 1, &infeasible, &fixed) );
            }

         }
      }
      return SCIP_OKAY;
   }
}
/** Gets the values of the pedigree specific variables in the current solution.
 *
 *  If no pedigree specific variables are being used, NULL is returned.
 *
 *  @param scip The SCIP instance for which to find the variable values.
 *  @return A list of SCIP_Bools which are TRUE for all variables included in the solution and FALSE otherwise.
 */
SCIP_Bool* PD_getCurrentPedigreeVarValues(
   SCIP* scip
   )
{
   if( usingSexConsistency(scip) )
   {
      int i;
      SCIP_Bool* vals;
      SCIP_SOL* sol = SCIPgetBestSol(scip);
      PedigreeData* peddata = MD_getPedigreeData(scip);
      SCIPallocMemoryArray(scip, &vals, peddata->n);
      for( i = 0; i < peddata->n; i++ )
         if( SCIPgetSolVal(scip, sol, peddata->SexVars[i]) > 0.5 )
            vals[i] = 1;
         else
            vals[i] = 0;
      return vals;
   }
   else
   {
      return NULL;
   }
}

/* Functions that create constraints */
/** Creates constraints stating that individuals cannot have more than a given number of children.
 *
 *  @param scip The SCIP instance this applies in.
 *  @param psd The parent set data associated with this problem.
 *  @param max_size The maximum number of children an indiivdual may have.
 *  @param mother Whether the constraint should apply to mothers.
 *  @param father Whether the constraint should apply to fathers.
 *  @param peddata The pedigree data to use to construct the constraint.
 *  @return SCIP_OKAY if adding the constraints worked.  An appropriate error code otherwise.
 */
static
SCIP_RETCODE addMaximumNumberOfOffspringConstraint(SCIP* scip,
   ParentSetData* psd,
   int max_size,
   SCIP_Bool mother,
   SCIP_Bool father,
   PedigreeData* peddata
   )
{
   SCIP_CONS* cons;
   char s[SCIP_MAXSTRLEN];
   int i;
   int j;
   int k;
   int l;

   /* Constraints must apply to either mothers or fathers or both. */
   if( !mother && !father )
   {
      SCIPerrorMessage("Maximum number of offspring constraint must apply to mothers, fathers or both.\n");
      return SCIP_ERROR;
   }

   for( i = 0; i < peddata->n; i++ )
   {
      int slack = max_size;
      if( mother && father )
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum offspring of %s", psd->nodeNames[i]);
      else if( mother && !father )
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum offspring of %s (as mother)", psd->nodeNames[i]);
      else
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum offspring of %s (as father)", psd->nodeNames[i]);
      SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons, s, max_size) );
      for( j = 0; j < peddata->n; j++ )
         for( k = 0; k < psd->nParentSets[j]; k++ )
            for( l = 0; l < psd->nParents[j][k]; l++ )
               if( i == psd->ParentSets[j][k][l] )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[j][k], 1) );
                  slack--;
               }
      if( mother && !father )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, peddata->SexVars[i], peddata->n) );
         SCIP_CALL( SCIPchgRhsLinear(scip, cons, SCIPgetRhsLinear(scip, cons) + peddata->n) );
      }
      else if( !mother && father )
      {
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, peddata->SexVars[i], -peddata->n) );
      }
      if( slack < 0 )
      {
         /* Only add constraint if it is not always satisfied */
         SCIP_CALL( SCIPaddCons(scip, cons) );
      }
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** Creates constraints stating that a pair of individuals cannot have more than a given number of children together.
 *
 *  @param scip The SCIP instance this applies in.
 *  @param psd The parent set data associated with this problem.
 *  @param max_size The maximum number of children an indiivdual may have.
 *  @param peddata The pedigree data to use to construct the constraint.
 *  @return SCIP_OKAY if adding the constraints worked.  An appropriate error code otherwise.
 */
static
SCIP_RETCODE addMaximumSibsetSizeConstraint(
   SCIP* scip,
   ParentSetData* psd,
   int max_size,
   PedigreeData* peddata
   )
{
   SCIP_CONS* cons;
   char s[SCIP_MAXSTRLEN];
   int par1;
   int par2;
   int i;
   int j;

   for( par1 = 0; par1 < peddata->n - 1; par1++ )
      for( par2 = par1 + 1; par2 < peddata->n; par2++ )
      {
         int terms = 0;
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum family size of_%s and %s", psd->nodeNames[par1], psd->nodeNames[par2]);
         SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons, s, max_size) );
         for( i = 0; i < peddata->n; i++ )
            for( j = 0; j < psd->nParentSets[i]; j++ )
               if( psd->nParents[i][j] == 2 )
                  if((par1 == psd->ParentSets[i][j][0] && par2 == psd->ParentSets[i][j][1]) ||
                        (par1 == psd->ParentSets[i][j][1] && par2 == psd->ParentSets[i][j][0]))
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][j], 1) );
                     terms++;
                  }
         if( terms > max_size )
         {
            /* Constraint isn't trivially satisfied */
            SCIP_CALL( SCIPaddCons(scip, cons) );
         }
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

   return SCIP_OKAY;
}
/** Creates constraints which prevent the age gap better a pair of full siblings being too great.
 *
 *  @param scip The SCIP instance in which to apply the constraint.
 *  @param psd The parent set data associated with this problem.
 *  @param max_age_gap The maximum permissible age gap.
 *  @param peddata The pedigree data to use to construct the constraint.
 *  @return SCIP_OKAY if all needed constraints were created successfully.  An appropriate error code otherwise.
 */
static
SCIP_RETCODE addFullSiblingAgeGapConstraint(
   SCIP* scip,
   ParentSetData* psd,
   int max_age_gap,
   PedigreeData* peddata
   )
{
   SCIP_CONS* cons;
   char s[SCIP_MAXSTRLEN];
   int par1;
   int par2;
   int child1;
   int child2;
   int i;
   int j;

   for( child1 = 0; child1 < peddata->n - 1; child1++ )
      for( i = 0; i < psd->nParentSets[child1]; i++ )
         if( psd->nParents[child1][i] == 2 )
         {
            par1 = psd->ParentSets[child1][i][0];
            par2 = psd->ParentSets[child1][i][1];
            for( child2 = child1 + 1; child2 < peddata->n; child2++ )
               for( j = 0; j < psd->nParentSets[child2]; j++ )
                  if( psd->nParents[child2][j] == 2 )
                     if((par1 == psd->ParentSets[child2][j][0] && par2 == psd->ParentSets[child2][j][1]) ||
                           (par1 == psd->ParentSets[child2][j][1] && par2 == psd->ParentSets[child2][j][0]))
                     {
                        /* Found a matching parent set */
                        (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum sibling age difference for %s and %s as children of %s and %s", psd->nodeNames[child1], psd->nodeNames[child2], psd->nodeNames[par1], psd->nodeNames[par2]);
                        SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons, s, max_age_gap + peddata->ages[child1] + peddata->ages[child2]) );
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[child1][i], peddata->ages[child1]) );
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[child2][j], peddata->ages[child1]) );
                        SCIP_CALL( SCIPaddCons(scip, cons) );
                        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                        (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum sibling age difference for %s and %s as children of %s and %s", psd->nodeNames[child2], psd->nodeNames[child1], psd->nodeNames[par1], psd->nodeNames[par2]);
                        SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons, s, max_age_gap + peddata->ages[child1] + peddata->ages[child2]) );
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[child1][i], peddata->ages[child2]) );
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[child2][j], peddata->ages[child2]) );
                        SCIP_CALL( SCIPaddCons(scip, cons) );
                        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                     }

         }

   return SCIP_OKAY;
}

/** Creates constraints which prevent the age gap better a pair of half siblings with a common mother being too great.
 *
 *  @param scip The SCIP instance in which to apply the constraint.
 *  @param psd The parent set data associated with this problem.
 *  @param max_age_gap The maximum permissible age gap.
 *  @param peddata The pedigree data to use to construct the constraint.
 *  @return SCIP_OKAY if all needed constraints were created successfully.  An appropriate error code otherwise.
 */
static
SCIP_RETCODE addHalfSiblingAgeGapConstraint(
   SCIP* scip,
   ParentSetData* psd,
   int max_age_gap,
   PedigreeData* peddata
   )
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   char s[SCIP_MAXSTRLEN];
   int child1;
   int child2;
   int par;
   int i;
   int j;

   for( child1 = 0; child1 < peddata->n - 1; child1++ )
      for( child2 = child1 + 1; child2 < peddata->n; child2++ )
         for( par = 0; par < peddata->n; par++ )
         {
            SCIP_Bool isEverParentOfChild1 = FALSE;
            SCIP_Bool isEverParentOfChild2 = FALSE;

            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum half sibling age difference for %s and %s with %s as mother", psd->nodeNames[child1], psd->nodeNames[child2], psd->nodeNames[par]);
            SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons1, s, max_age_gap + peddata->ages[child1] + peddata->ages[child2]) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons1, peddata->SexVars[par], peddata->ages[child1]) );

            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "maximum half sibling age difference for %s and %s with %s as mother", psd->nodeNames[child2], psd->nodeNames[child1], psd->nodeNames[par]);
            SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons2, s, max_age_gap + peddata->ages[child1] + peddata->ages[child2]) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons2, peddata->SexVars[par], peddata->ages[child2]) );

            for( i = 0; i < psd->nParentSets[child1]; i++ )
            {
               SCIP_Bool parInThisSet = FALSE;
               for( j = 0; j < psd->nParents[child1][i]; j++ )
                  if( par == psd->ParentSets[child1][i][j] )
                  {
                     /* parent is in this parent set of child 1 */
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons1, psd->PaVars[child1][i], peddata->ages[child1]) );
                     parInThisSet = TRUE;
                     isEverParentOfChild1 = TRUE;
                  }
               if( parInThisSet == FALSE )
               {
                  /* parent is not in this parent set of child 1 */
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons2, psd->PaVars[child1][i], peddata->ages[child2]) );
               }
            }

            for( i = 0; i < psd->nParentSets[child2]; i++ )
            {
               SCIP_Bool parInThisSet = FALSE;
               for( j = 0; j < psd->nParents[child2][i]; j++ )
                  if( par == psd->ParentSets[child2][i][j] )
                  {
                     /* parent is in this parent set of child 2 */
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons2, psd->PaVars[child2][i], peddata->ages[child2]) );
                     parInThisSet = TRUE;
                     isEverParentOfChild2 = TRUE;
                  }
               if( parInThisSet == FALSE )
               {
                  /* parent is not in this parent set of child 2 */
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons1, psd->PaVars[child2][i], peddata->ages[child1]) );
               }
            }

            /* Only add the constraints if the children can be half-sibs */
            if( isEverParentOfChild1 && isEverParentOfChild2 )
            {
               SCIP_CALL( SCIPaddCons(scip, cons1) );
               SCIP_CALL( SCIPaddCons(scip, cons2) );
            }

            SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
         }

   return SCIP_OKAY;
}

/** Creates constraints enforcing sexual consistency on a pedigree.
 *
 *  @param scip The SCIP instance in which the constraints should be added.
 *  @param psd The parent set data associated with this problem.
 *  @param peddata The pedigree data to use to construct the constraint.
 *  @return SCIP_OKAY if the constraints were added successfully, or an appropriate error code otherwsie.
 */
static
SCIP_RETCODE addSexConsistencyConstraint(
   SCIP* scip,
   ParentSetData* psd,
   PedigreeData* peddata
   )
{
   SCIP_CONS* cons;
   char s[SCIP_MAXSTRLEN];
   int i;
   int j;

   for( i = 0; i < peddata->n; i++ )
      for( j = 0; j < psd->nParentSets[i]; j++ )
         if( psd->nParents[i][j] == 2 )
         {
            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sex_consistency_1_for_%s", SCIPvarGetName(psd->PaVars[i][j]));
            SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons, s, 2) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][j], 1) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, peddata->SexVars[psd->ParentSets[i][j][0]], 1) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, peddata->SexVars[psd->ParentSets[i][j][1]], 1) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sex_consistency_2_for_%s", SCIPvarGetName(psd->PaVars[i][j]));
            SCIP_CALL( UT_createEmptyLTEConstraint(scip, &cons, s, 0) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][j], 1) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, peddata->SexVars[psd->ParentSets[i][j][0]], -1) );
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, peddata->SexVars[psd->ParentSets[i][j][1]], -1) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
   return SCIP_OKAY;
}

/** Adds all appropriate pedigree based constraints to the problem.
 *
 *  @param scip The SCIP instance in which to add the constraints.
 *  @param psd The parent set data associated with this problem.
 *  @param peddata The pedigree data to use to construct the constraints.
 *  @return SIP_OKAY if the constraints were added successfully, or an appropriate error code otherwise.
 */
static
SCIP_RETCODE PD_addPedigreeConstraints(
   SCIP* scip,
   ParentSetData* psd,
   PedigreeData* peddata
   )
{
   if( !checkSuitableForPedigree(scip, psd) )
   {
      SCIPerrorMessage("The specified problem is not a valid pedigree problem.\n");
      return SCIP_ERROR;
   }
   else
   {
      int maxsibagegap;
      int maxhalfsibagegap;
      int maxsibsetsize;
      int maxchildren;
      int maxchildrenmother;
      int maxchildrenfather;
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pedigree/maxsibagegap",      &maxsibagegap) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pedigree/maxhalfsibagegap",  &maxhalfsibagegap) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pedigree/maxsibsetsize",     &maxsibsetsize) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pedigree/maxchildren",       &maxchildren) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pedigree/maxchildrenmother", &maxchildrenmother) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pedigree/maxchildrenfather", &maxchildrenfather) );
      if( maxsibagegap != -1 )
      {
         if( peddata->sexes == NULL )
         {
            SCIPerrorMessage("Can't enforce maximum sibling age gap when no age data is given.\n");
            return SCIP_ERROR;
         }
         else
         {
            SCIP_CALL( addFullSiblingAgeGapConstraint(scip, psd, maxsibagegap, peddata) );
         }
      }
      if( maxhalfsibagegap != -1 )
      {
         if( peddata->sexes == NULL )
         {
            SCIPerrorMessage("Can't enforce maximum maternal half-sibling age gap when no age data is given.\n");
            return SCIP_ERROR;
         }
         else if( !usingSexConsistency(scip) )
         {
            SCIPerrorMessage("Can't enforce maximum maternal half-sibling age gap without sexual consistency.\n");
            return SCIP_ERROR;
         }
         else
         {
            SCIP_CALL( addHalfSiblingAgeGapConstraint(scip, psd, maxhalfsibagegap, peddata) );
         }
      }
      if( maxsibsetsize != -1 )
         SCIP_CALL( addMaximumSibsetSizeConstraint(scip, psd, maxsibsetsize, peddata) );
      if( maxchildren != -1 )
         SCIP_CALL( addMaximumNumberOfOffspringConstraint(scip, psd, maxchildren, TRUE, TRUE, peddata) );
      if( maxchildrenmother != -1 )
      {
         if( usingSexConsistency(scip) )
            SCIP_CALL( addMaximumNumberOfOffspringConstraint(scip, psd, maxchildrenmother, TRUE, FALSE, peddata) );
         else
         {
            SCIPerrorMessage("Can't use a limit on the maximum number of children a mother can have without using sexual consistency.\n");
            return SCIP_ERROR;
         }
      }
      if( maxchildrenfather != -1 )
      {
         if( usingSexConsistency(scip) )
            SCIP_CALL( addMaximumNumberOfOffspringConstraint(scip, psd, maxchildrenfather, FALSE, TRUE, peddata) );
         else
         {
            SCIPerrorMessage("Can't use a limit on the maximum number of children a father can have without using sexual consistency.\n");
            return SCIP_ERROR;
         }
      }
      if( usingSexConsistency(scip) )
         SCIP_CALL( addSexConsistencyConstraint(scip, psd, peddata) );
      return SCIP_OKAY;
   }
}

/** Adds any constraints and variables needed for pedigree based constraints.
 *
 *  @param scip The SCIP instance in which to add the constraints.
 *  @param psd The parent set data associated with this problem.
 *  @return SIP_OKAY if the constraints were added successfully, or an appropriate error code otherwise.
 */
SCIP_RETCODE PD_addPedigreeSpecificConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   PedigreeData* peddata = MD_getPedigreeData(scip);
   SCIP_Bool data_already_existed = (peddata != NULL);

   if( !data_already_existed )
   {
      SCIP_CALL( SCIPallocMemory(scip, &peddata) );
      peddata->n = psd->n;
      peddata->SexVars = NULL;
      peddata->ages = NULL;
      peddata->sexes = NULL;
   }

   SCIP_CALL( PD_addPedigreeVariables(scip, psd, peddata) );
   SCIP_CALL( PD_addPedigreeConstraints(scip, psd, peddata) );

   if( !data_already_existed )
   {
      SCIP_CALL( MD_setPedigreeData(scip, peddata) );
      SCIP_CALL( PE_deallocatePedigreeData(scip, &peddata) );
   }
   return SCIP_OKAY;
}

/* Functions related to the primal heuristic */
/** Determines a possible assignment of the sex variables for a given primal solution.
 *
 *  If there is no possible assignment of sex variables (i.e. the primal is sexually
 *  inconsistent) then the function returns an appropriate error message.
 *
 *  @param scip The SCIP instance on which the heuristic is running.
 *  @param sol The heuristic solution being worked on.
 *  @param psd The heursitic data related to this primal heursitic.
 *  @param possible Will be set to TRUE if an assignment was possible.
 *  @return SCIP_OKAY if a consistent labelling could be found and was made to sol.
 *  An appropriate error code otherwise.
 */
static
SCIP_RETCODE assignSexVariables(
   SCIP* scip,
   SCIP_SOL* sol,
   ParentSetData* psd,
   SCIP_Bool* possible
   )
{
   int i;
   int k;
   PedigreeData* peddata = MD_getPedigreeData(scip);
   int numTwoParents = 0;              /* The number of individuals with two parents who haven't yet had sexes assigned */
   SCIP_Bool* isFinished;              /* TRUE iff we have dealt with all the parents of the node */
   SCIP_Bool* isSexed;                 /* TRUE iff a sex has been assigned to the node */
   SCIP_Bool* isFemale;                /* TRUE if the node has been made female, FALSE if it has been made male and undefined if isSexed[i] == FALSE */
   int* chosenSet;                     /* The parent set of i that is selected */
   SCIP_Bool no_parents;
   SCIP_Bool assignment_made;
   int val;

   SCIP_CALL( SCIPallocMemoryArray(scip, &isFinished, peddata->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &isSexed, peddata->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &isFemale, peddata->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &chosenSet, peddata->n) );

   /* Initialise the data structures */
   for( i = 0; i < peddata->n; i++ )
   {
      isFinished[i] = FALSE;
      isSexed[i] = FALSE;
   }

   /* Check if any sexes are already fixed */
   for( i = 0; i < peddata->n; i++ )
      if( SCIPvarGetLbLocal(peddata->SexVars[i]) > 0.5 )
      {
         isSexed[i] = TRUE;
         isFemale[i] = TRUE;
      }
      else if( SCIPvarGetUbLocal(peddata->SexVars[i]) < 0.5 )
      {
         isSexed[i] = TRUE;
         isFemale[i] = FALSE;
      }

   /* Quick pass through to find individuals with two parents */
   for( i = 0; i < peddata->n; ++i )
   {
      no_parents = TRUE;
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         val = SCIPgetSolVal(scip, sol, psd->PaVars[i][k]);
         assert(SCIPisIntegral(scip, val));
         if( val > 0.5 )
         {
            chosenSet[i] = k;
            if( psd->nParents[i][k] == 2 )
               numTwoParents++;
            else
               isFinished[i] = TRUE;
            no_parents = FALSE;
            break;
         }
      }
      if( no_parents == TRUE )
         isFinished[i] = TRUE;
   }
   /* Main loop */
   /* We assign sexes to any parent whose mate is already sexed */
   /* If that is not possible, we randomly assign sex to any individual */
   while( numTwoParents > 0 )
   {
      assignment_made = FALSE;
      for( i = 0; i < peddata->n; i++ )
      {
         if( isFinished[i] == FALSE )
         {
            int       parent1 = psd->ParentSets[i][chosenSet[i]][0];
            int       parent2 = psd->ParentSets[i][chosenSet[i]][1];
            SCIP_Bool known1  = isSexed[parent1];
            SCIP_Bool known2  = isSexed[parent2];
            if( known1 && known2 )
            {
               if( isFemale[parent1] == isFemale[parent2] )
               {
                  /* SEXUAL INCONSISTENT ASSIGNMENT */
                  /* Just exit and let the calling procedure deal with it */
                  *possible = FALSE;
                  return SCIP_OKAY;
               }
               else
               {
                  /* Just make it as done */
                  isFinished[i] = TRUE;
                  numTwoParents--;
               }
            }
            else if( known1 && !known2 )
            {
               /* Assign parent2 the opposite sex */
               if( isFemale[parent1] )
                  SCIP_CALL( SCIPsetSolVal(scip, sol, peddata->SexVars[parent2], 0) );
               else
                  SCIP_CALL( SCIPsetSolVal(scip, sol, peddata->SexVars[parent2], 1) );
               isFemale[parent2] = !isFemale[parent1];
               isSexed[parent2] = TRUE;
               isFinished[i] = TRUE;
               numTwoParents--;
               assignment_made = TRUE;
            }
            else if( !known1 && known2 )
            {
               /* Assign parent1 the opposite sex */
               if( isFemale[parent2] )
                  SCIP_CALL( SCIPsetSolVal(scip, sol, peddata->SexVars[parent1], 0) );
               else
                  SCIP_CALL( SCIPsetSolVal(scip, sol, peddata->SexVars[parent1], 1) );
               isFemale[parent1] = !isFemale[parent2];
               isSexed[parent1] = TRUE;
               isFinished[i] = TRUE;
               numTwoParents--;
               assignment_made = TRUE;
            }
            else
            {
               /* Can't do anything */
            }
         }
      }
      if( assignment_made == FALSE && numTwoParents > 0 )
      {
         /* Need to make a random assignment */
         /* Make the first unsexed individual female */
         i = 0;
         while( assignment_made == FALSE )
         {
            if( isSexed[i] == FALSE )
            {
               SCIP_CALL( SCIPsetSolVal(scip, sol, peddata->SexVars[i], 1) );
               isSexed[i] = TRUE;
               isFemale[i] = TRUE;
               assignment_made = TRUE;
            }
            i++;
         }
      }
   }

   /*Go through and assign sexes to anyone still unsexed now */
   for( i = 0; i < peddata->n; i++ )
      if( isSexed[i] == FALSE )
         SCIP_CALL( SCIPsetSolVal(scip, sol, peddata->SexVars[i], 1) );

   SCIPfreeMemoryArray(scip, &isFinished);
   SCIPfreeMemoryArray(scip, &isSexed);
   SCIPfreeMemoryArray(scip, &isFemale);
   SCIPfreeMemoryArray(scip, &chosenSet);

   /* All nodes have been assigned sexes in a consistent way */
   *possible = TRUE;
   return SCIP_OKAY;
}

/** Tries to find values for any pedigree specific variables in a given primal solution.
 *
 *  If there is no possible assignment of pedigree variables then the function returns
 *  an appropriate error message.
 *
 *  @param scip The SCIP instance on which the heuristic is running.
 *  @param sol The heuristic solution being worked on.
 *  @param psd The heursitic data related to this primal heursitic.
 *  @param possible Returns whether a valid assignment was possible.
 *  @return SCIP_OKAY if a consistent labelling could be found and was made to sol.
 *  An appropriate error code otherwise.
 */
SCIP_RETCODE PD_assignPedigreeVariables(
   SCIP* scip,
   SCIP_SOL* sol,
   ParentSetData* psd,
   SCIP_Bool* possible
   )
{
   if( !checkSuitableForPedigree(scip, psd) )
   {
      SCIPerrorMessage("The specified problem is not a valid pedigree problem.\n");
      *possible = FALSE;
      return SCIP_ERROR;
   }
   else
   {
      *possible = TRUE;
      if( usingSexConsistency(scip) )
      {
         SCIP_Bool sex_possible = TRUE;
         SCIP_CALL( assignSexVariables(scip, sol, psd, &sex_possible) );
         if( sex_possible == FALSE )
            *possible = FALSE;
      }
      return SCIP_OKAY;
   }
}

/* Functions for output */
/** Prints the solution as a pedigree.
 *
 *  The pedigree consists of four columns.  The first is the individual, thes second is the sex of the individual, the third
 *  is its father and the fourth is its mother.  If either the father or the mother is not present in the sample, a '-' is
 *  printed instead.  If sex consistency is enforced, then individuals will not appear as father of one individual and mother
 *  of another.  The pedigree is sorted such that all individuals are declared on a earlier line than those in which they
 *  appear as parents.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @param pedvars Whether each of the pedigree specific variables is selected in the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
SCIP_RETCODE PD_printSolutionPedigreeFormat(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_Real** Scores,
   SCIP_Bool** selected,
   SCIP_Real total_score,
   FILE* stream,
   SCIP_Bool* pedvars
   )
{
   if( !checkSuitableForPedigree(scip, psd) )
   {
      SCIPerrorMessage("The specified problem is not a valid pedigree problem.\n");
      return SCIP_ERROR;
   }
   else
   {
      int i, k, k2;
      SCIP_Bool no_parents;
      SCIP_Bool **Done;
      SCIP_Bool *allDone;
      SCIP_Bool *someDone;
      int numDone = 0;
      char sex;
      PedigreeData* peddata = MD_getPedigreeData(scip);
      SCIP_CALL( SCIPallocMemoryArray(scip, &someDone, peddata->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &allDone, peddata->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &Done, peddata->n) );
      for( i = 0; i < peddata->n; i++ )
      {
         allDone[i] = TRUE;
         someDone[i] = FALSE;
         SCIP_CALL( SCIPallocMemoryArray(scip, &(Done[i]), psd->nParentSets[i]) );
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            if( selected[i][k] )
            {
               Done[i][k] = FALSE;
               allDone[i] = FALSE;
            }
            else
               Done[i][k] = TRUE;
         }
      }

      while( numDone < peddata->n )
      {
         for( i = 0; i < peddata->n; ++i )
         {
            if( allDone[i] == FALSE )
            {
               no_parents = TRUE;
               if( usingSexConsistency(scip) == FALSE || pedvars == NULL )
                  sex = 'U';
               else if( pedvars[i] )
                  sex = 'F';
               else
                  sex = 'M';
               for( k = 0; k < psd->nParentSets[i]; ++k )
               {
                  if( selected[i][k] && Done[i][k] == FALSE )
                  {
                     no_parents = FALSE;
                     if( psd->nParents[i][k] == 0 )
                     {
                        fprintf(stream, "%s\t%c\t-\t-\n", psd->nodeNames[i], sex);
                        someDone[i] = TRUE;
                        Done[i][k] = TRUE;
                        allDone[i] = TRUE;
                        for( k2 = 0; k2 < psd->nParentSets[i]; ++k2 )
                           if( !Done[i][k2] )
                              allDone[i] = FALSE;
                        if( allDone[i] )
                           numDone++;
                     }
                     else if( psd->nParents[i][k] == 1 )
                     {
                        if( someDone[psd->ParentSets[i][k][0]] == TRUE )
                        {
                           if( !usingSexConsistency(scip) || pedvars == NULL || !pedvars[psd->ParentSets[i][k][0]] )
                              fprintf(stream, "%s\t%c\t%s\t-\n", psd->nodeNames[i], sex, psd->nodeNames[psd->ParentSets[i][k][0]]);
                           else
                              fprintf(stream, "%s\t%c\t-\t%s\n", psd->nodeNames[i], sex, psd->nodeNames[psd->ParentSets[i][k][0]]);
                           someDone[i] = TRUE;
                           Done[i][k] = TRUE;
                           allDone[i] = TRUE;
                           for( k2 = 0; k2 < psd->nParentSets[i]; ++k2 )
                              if( !Done[i][k2] )
                                 allDone[i] = FALSE;
                           if( allDone[i] )
                              numDone++;
                        }
                     }
                     else
                     {
                        if( someDone[psd->ParentSets[i][k][0]] == TRUE && someDone[psd->ParentSets[i][k][1]] == TRUE )
                        {
                           if( !usingSexConsistency(scip) || pedvars == NULL || pedvars[psd->ParentSets[i][k][0]] )
                              fprintf(stream, "%s\t%c\t%s\t%s\n", psd->nodeNames[i], sex, psd->nodeNames[psd->ParentSets[i][k][1]], psd->nodeNames[psd->ParentSets[i][k][0]]);
                           else
                              fprintf(stream, "%s\t%c\t%s\t%s\n", psd->nodeNames[i], sex, psd->nodeNames[psd->ParentSets[i][k][0]], psd->nodeNames[psd->ParentSets[i][k][1]]);
                           someDone[i] = TRUE;
                           Done[i][k] = TRUE;
                           allDone[i] = TRUE;
                           for( k2 = 0; k2 < psd->nParentSets[i]; ++k2 )
                              if( !Done[i][k2] )
                                 allDone[i] = FALSE;
                           if( allDone[i] )
                              numDone++;
                        }
                     }
                  }
               }
               if( no_parents )
               {
                  fprintf(stream, "%s\t%c\t-\t-\n", psd->nodeNames[i], sex);
                  someDone[i] = TRUE;
                  allDone[i] = TRUE;
                  for( k2 = 0; k2 < psd->nParentSets[i]; ++k2 )
                     if( !Done[i][k2] )
                        allDone[i] = FALSE;
                  if( allDone[i] )
                     numDone++;
               }
            }
         }
      }

      SCIPfreeMemoryArray(scip, &someDone);
      SCIPfreeMemoryArray(scip, &allDone);
      for( i = 0; i < peddata->n; i++ )
         SCIPfreeMemoryArray(scip, &(Done[i]));
      SCIPfreeMemoryArray(scip, &Done);

      return SCIP_OKAY;
   }
}
