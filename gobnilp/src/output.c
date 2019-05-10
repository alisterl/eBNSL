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
 *  Contains all the functions related to writing data to files.
 */

#include "output.h"
#include "parent_set_data.h"
#include "pedigrees.h"
#include "utils.h"
#include "probdata_bn.h"
#include "versiongit.h"
#include <string.h>
#include <time.h>
#include "scip/scipdefplugins.h"
#include <scip/cons_countsols.h>

/** Adds the parameters to SCIP related to outputing results.
 *
 *  These parameters consist of all those of the form @c gobnilp/outputfile/\<param\> .
 *  Some other output parameters exist, but these can be currently found in probdata_bn.c .
 *
 * @param scip The SCIP instance to which parameters should be added.
 * @return SCIP_OKAY if all the parameters were added successfully, or an error otherwise.
 */
SCIP_RETCODE IO_addOutputParameters(SCIP* scip)
{
   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/solution",
                               "file which solution should be printed to (stdout for standard out, empty string for nowhere)",
                               "stdout"
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/dot",
                               "file which dot output should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/pedigree",
                               "file which pedigree output should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/scoreandtime",
                               "file which additional score and time data should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/adjacencymatrix",
                               "file which adjacency matrix output should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/problem",
                               "file which problem representation should be printed to (stdout for standard out, empty string for nowhere). File extension determines format.",
                               ""
                              ));

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/printmeta",
         "whether to print out dummy metadata constraints when printing the problem to file",
         TRUE
         ) );


   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/avgoutputoffset",
                            "how many iterations to skip before first model averaging output (-1 to suppress always)",
                            -1, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/avgoutputstep",
                            "how many iterations between outputting model averaging information",
                            1, 1, INT_MAX
                           ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/solutionavg",
                               "file which model averging solution should be printed to (stdout for standard out, empty string for nowhere)",
                               "stdout"
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/dotavg",
                               "file which model averging dot output should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/pedigreeavg",
                               "file which model averging pedigree output should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/scoreandtimeavg",
                               "file which model averging additional score and time data should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/adjacencymatrixavg",
                               "file which model averging adjacency matrix output should be printed to (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/scores",
                               "file which parent set scores should be output (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
                               "gobnilp/outputfile/mec",
                               "file which Markov equivalence class should be output (stdout for standard out, empty string for nowhere)",
                               ""
                              ));

   SCIP_CALL(UT_addStringParam(scip,
         "gobnilp/outputfile/pss",
         "file to which parent set scores should be output (stdout for standard out, empty string for nowhere)",
         ""
         ));


   SCIP_CALL( UT_addStringParam(scip,
         "gobnilp/outputfile/countsols",
         "collected feasible solutions are sent to this file",
         ""
         ) );


   return SCIP_OKAY;
}

/** Prints the solution in the traditional GOBNILP format.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionLegacyFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, SCIP_Real total_score, FILE* stream)
{
   int i, k, l;
   SCIP_Bool no_parents;

   for( i = 0; i < psd->n; ++i )
   {
      no_parents = TRUE;
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         if( selected[i][k] )
         {
            fprintf(stream, "%s<-", psd->nodeNames[i]);
            for( l = 0; l < psd->nParents[i][k]; ++l )
               fprintf(stream, "%s,", psd->nodeNames[psd->ParentSets[i][k][l]]);
            fprintf(stream, " %f\n", Scores[i][k]);
            no_parents = FALSE;
         }
      }
      if( no_parents )
      {
         fprintf(stream, " ,");
         fprintf(stream, " %f\n", Scores[i][psd->nParentSets[i] - 1]);
      }
   }
   fprintf(stream, "BN score is %f\n", total_score);
   return SCIP_OKAY;
}

/** Prints the solution in a Bayesian network format.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionBNFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, SCIP_Real total_score, FILE* stream)
{
   int i, k, l;
   SCIP_Bool no_parents;

   for( i = 0; i < psd->n; ++i )
   {
      SCIP_Bool first_line = TRUE;
      no_parents = TRUE;
      fprintf(stream, "%s\t<-", psd->nodeNames[i]);
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         if( selected[i][k] )
         {
            if( first_line )
               first_line = FALSE;
            else
               fprintf(stream, "\t");
            fprintf(stream, "\t{");
            for( l = 0; l < psd->nParents[i][k] - 1; ++l )
            {
               fprintf(stream, "%s,", psd->nodeNames[psd->ParentSets[i][k][l]]);
            }
            if( psd->nParents[i][k] > 0 )
               fprintf(stream, "%s", psd->nodeNames[psd->ParentSets[i][k][psd->nParents[i][k] - 1]]);
            fprintf(stream, "}\t%f\n", Scores[i][k]);
            no_parents = FALSE;
         }
      }
      if( no_parents )
         fprintf(stream, "{}\t%f\n", Scores[i][psd->nParentSets[i] - 1]);
   }
   fprintf(stream, "BN score is %f\n", total_score);
   return SCIP_OKAY;
}

/** Prints the solution as a file suitable for plotting using the dot command from graphviz.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionDotFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, SCIP_Real total_score, FILE* stream)
{
   int i, j, k, l;
   SCIP_Real **matrix;
   SCIP_Bool one_per_row;

   SCIP_CALL( SCIPallocMemoryArray(scip, &matrix, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(matrix[i]), psd->n) );
      for( j = 0; j < psd->n; ++j )
         matrix[i][j] = 0;
   }

   one_per_row = TRUE;
   for( j = 0; j < psd->n; ++j )
   {
      int entries_this_row = 0;
      for( k = 0; k < psd->nParentSets[j]; ++k )
         if( selected[j][k] )
            entries_this_row += 1;
      if( entries_this_row != 1 )
         one_per_row = FALSE;
   }

   for( j = 0; j < psd->n; ++j )
      for( k = 0; k < psd->nParentSets[j]; ++k )
         if( selected[j][k] )
         {
            for( l = 0; l < psd->nParents[j][k]; ++l )
            {
               if( one_per_row )
                  matrix[psd->ParentSets[j][k][l]][j] = 1 ;
               else
                  matrix[psd->ParentSets[j][k][l]][j] += Scores[j][k] ;
            }
         }

   fprintf(stream, "digraph {\n");
   for( i = 0; i < psd->n; ++i )
      fprintf(stream, "   %s;\n", psd->nodeNames[i]);
   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->n; ++k )
         if( matrix[i][k] != 0 )
         {
            fprintf(stream, "   %s -> %s", psd->nodeNames[i], psd->nodeNames[k]);
            if( !one_per_row )
               fprintf(stream, "[penwidth=%d]", (int)(10 * matrix[i][k] + 0.4999999));
            fprintf(stream, ";\n");
         }
   fprintf(stream, "}\n");


   for( i = 0 ; i < psd->n ; ++i )
      SCIPfreeMemoryArray(scip, &(matrix[i]));
   SCIPfreeMemoryArray(scip, &matrix);

   return SCIP_OKAY;
}

/** Prints the solution as an adjacency matrix.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionAdjacencyMatrixFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, SCIP_Real total_score, FILE* stream)
{
   int i, j, k, l;
   SCIP_Real **matrix;
   SCIP_Bool one_per_row;

   SCIP_CALL( SCIPallocMemoryArray(scip, &matrix, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(matrix[i]), psd->n) );
      for( j = 0; j < psd->n; ++j )
         matrix[i][j] = 0;
   }

   one_per_row = TRUE;
   for( j = 0; j < psd->n; ++j )
   {
      int entries_this_row = 0;
      for( k = 0; k < psd->nParentSets[j]; ++k )
         if( selected[j][k] )
            entries_this_row += 1;
      if( entries_this_row != 1 )
         one_per_row = FALSE;
   }

   for( j = 0; j < psd->n; ++j )
      for( k = 0; k < psd->nParentSets[j]; ++k )
         if( selected[j][k] )
         {
            for( l = 0; l < psd->nParents[j][k]; ++l )
            {
               if( one_per_row )
                  matrix[psd->ParentSets[j][k][l]][j] = 1 ;
               else
                  matrix[psd->ParentSets[j][k][l]][j] += Scores[j][k] ;
            }
         }

   for( i = 0; i < psd->n; ++i )
   {
      for( j = 0; j < (psd->n) - 1; ++j )
         if( one_per_row )
            fprintf(stream, "%d ", (int)matrix[i][j]);
         else
            fprintf(stream, "%f ", matrix[i][j]);
      if( one_per_row )
         fprintf(stream, "%d\n", (int)matrix[i][j]);
      else
         fprintf(stream, "%f\n", matrix[i][j]);
   }

   for( i = 0 ; i < psd->n ; ++i )
      SCIPfreeMemoryArray(scip, &(matrix[i]));
   SCIPfreeMemoryArray(scip, &matrix);

   return SCIP_OKAY;
}

/** Prints just the objective value of the given solution and the time taken to find the solution.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @param soltime The time taken to find the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionScoreAndTimeFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, SCIP_Real total_score, FILE* stream, SCIP_Real soltime)
{
   fprintf(stream, "%f\t%f\n", total_score, soltime);
   return SCIP_OKAY;
}


static SCIP_RETCODE printSCIPSolution(SCIP* scip)
{
   SCIP_SOL* sol;
   SCIP_Bool printscipsol;
   SCIPgetBoolParam(scip, "gobnilp/printscipsol", &printscipsol);
   sol = SCIPgetBestSol(scip);
   if( printscipsol )
      SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
   return SCIP_OKAY;
}


/** Prints the Markov equivalance class of the solution.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionMECFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, SCIP_Real total_score, FILE* stream)
{
   int i, j, k, l;

   SCIP_Real val = 0.0;

   int** edges;
   int parent1;
   int parent2;
   int ll;

   SCIP_CALL( SCIPallocMemoryArray(scip, &edges, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(edges[i]), psd->n) );
      for( j = 0; j < psd->n; ++j )
         edges[i][j] = FALSE;
   }

   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         val = Scores[i][k];
         if( val != 0 )
         {
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               j = psd->ParentSets[i][k][l];
               edges[i][j] = TRUE;
               edges[j][i] = TRUE;
            }
         }
      }


   for( i = 0; i < psd->n; ++i )
      for( j = i + 1; j < psd->n; ++j )
         if( edges[i][j] )
            fprintf(stream, "%s-%s\n", psd->nodeNames[i], psd->nodeNames[j]);

   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         val = Scores[i][k];
         if( val != 0 )
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               parent1 = psd->ParentSets[i][k][l];
               for( ll = l + 1; ll < psd->nParents[i][k]; ++ll )
               {
                  parent2 = psd->ParentSets[i][k][ll];
                  if( !edges[parent1][parent2] )
                     fprintf(stream, "%s->%s<-%s\n", psd->nodeNames[parent1], psd->nodeNames[i], psd->nodeNames[parent2]);
               }
            }
      }

   for( i = 0 ; i < psd->n ; ++i )
      SCIPfreeMemoryArray(scip, &(edges[i]));
   SCIPfreeMemoryArray(scip, &edges);

   return SCIP_OKAY;
}

/** Prints a solution to the problem.
 *
 *  @param scip The SCIP instance for which to print the solution.
 *  @param psd The problem data used by the solution.
 *  @param Scores The score data to use for the solution.
 *  @param selected Whether each of the variables is selected in the solution.
 *  @param total_score The overall score of this solution.
 *  @param filename The filename to output to, "stdout" for stdout or "" for nowhere.
 *  @param format The format in which to print the solution.  Recognised values are dot, pedigree, scoreandtime, legacy, adjacencymatrx and normal.
 *  @param soltime The time taken to find the solution.
 *  @param pedVals Whether each of the pedigree specific variables is selected in the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printToFile(SCIP* scip, ParentSetData* psd, SCIP_Real** Scores, SCIP_Bool** selected, char* filename, char* format, SCIP_Real soltime, SCIP_Real total_score, SCIP_Bool* pedVals)
{
   /* Open the file for writing */
   FILE* file;
   if( strcmp(filename, "") == 0 )
      return SCIP_OKAY;
   else if( strcmp(filename, "stdout") == 0 )
      file = stdout;
   else
   {
      file = fopen(filename, "w");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Writing output to %s\n", filename);
   }
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file %s for writing.\n", filename);
      return SCIP_WRITEERROR;
   }

   /* Print the solution to the file */
   if( strcmp(format, "dot") == 0 )
      printSolutionDotFormat(scip, psd, Scores, selected, total_score, file);
   else if( strcmp(format, "pedigree") == 0 )
      PD_printSolutionPedigreeFormat(scip, psd, Scores, selected, total_score, file, pedVals);
   else if( strcmp(format, "scoreandtime") == 0 )
      printSolutionScoreAndTimeFormat(scip, psd, Scores, selected, total_score, file, soltime);
   else if( strcmp(format, "legacy") == 0 )
      printSolutionLegacyFormat(scip, psd, Scores, selected, total_score, file);
   else if( strcmp(format, "adjacencymatrix") == 0 )
      printSolutionAdjacencyMatrixFormat(scip, psd, Scores, selected, total_score, file);
   else if( strcmp(format, "mec") == 0 )
      printSolutionMECFormat(scip, psd, Scores, selected, total_score, file);
   else
      printSolutionBNFormat(scip, psd, Scores, selected, total_score, file);

   /* Close the file */
   if( file != stdout )
      fclose(file);

   return SCIP_OKAY;
}

/** Prints the solution of the problem after solving.
 *
 *  @param scip The SCIP instance for which to print the solution.
 *  @param psd The parentage data for the problem.vvv
 *  @param filename The filename to output to, "stdout" for stdout or "" for nowhere.
 *  @param format The format in which to print the solution.  Recognised values are dot, pedigree, scoreandtime, legacy, adjacencymatrx and normal.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolution(SCIP* scip, ParentSetData* psd, char* filename, char* format)
{
   SCIP_SOL* sol = SCIPgetBestSol(scip);
   SCIP_Real** Scores;
   SCIP_Bool** selected;
   SCIP_Bool* pedVals = NULL;
   int i, j;
   SCIP_Real total_score = 0;

   if( psd != NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &selected, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &Scores, psd->n) );
      for( i = 0; i < psd->n; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(selected[i]), psd->nParentSets[i]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(Scores[i]), psd->nParentSets[i]) );
      }

      for( i = 0; i < psd->n; i++ )
         for( j = 0; j < psd->nParentSets[i]; j++ )
            if( SCIPgetSolVal(scip, sol, psd->PaVars[i][j]) > 0.5 )
            {
               selected[i][j] = TRUE;
               Scores[i][j] = SCIPvarGetObj(psd->PaVars[i][j]);
               total_score += Scores[i][j];
            }
            else
            {
               selected[i][j] = FALSE;
               Scores[i][j] = 0;
            }

      if( PD_inPedigreeMode(scip) )
         pedVals = PD_getCurrentPedigreeVarValues(scip);

      SCIP_CALL( printToFile(scip, psd, Scores, selected, filename, format, SCIPgetSolvingTime(scip), total_score, pedVals) );

      for( i = 0; i < psd->n; i++ )
      {
         SCIPfreeMemoryArray(scip, &(Scores[i]));
         SCIPfreeMemoryArray(scip, &(selected[i]));
      }
      SCIPfreeMemoryArray(scip, &Scores);
      SCIPfreeMemoryArray(scip, &selected);
      if( pedVals != NULL )
         SCIPfreeMemoryArray(scip, &pedVals);
   }

   return SCIP_OKAY;
}

/** Prints the model average data.
 *
 *  @param scip The SCIP instance for which to print the averages.
 *  @param psd The parentage data for the problem.
 *  @param filename The filename to output to, "stdout" for stdout or "" for nowhere.
 *  @param format The format in which to print the averages.  Recognised values are dot, pedigree, scoreandtime, legacy, adjacencymatrx and normal.
 *  @return SCIP_OKAY if the averages were printed correctly or an appropriate error message otherwise.
 */
static
SCIP_RETCODE printAverages(
   SCIP* scip,
   MA_info* ma_info,
   ParentSetData* psd,
   char* filename,
   char* format
   )
{
   int i;
   int j;
   SCIP_Real** Scores;
   SCIP_Bool** selected;

   if( psd != NULL )
   {
      /* Allocate memory */
      SCIP_CALL( SCIPallocMemoryArray(scip, &Scores, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &selected, psd->n) );
      for( i = 0; i < psd->n; i++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(Scores[i]), psd->nParentSets[i]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(selected[i]), psd->nParentSets[i]) );
      }

      /* Get data and print it */
      for( i = 0; i < psd->n; i++ )
         for( j = 0; j < psd->nParentSets[i]; j++ )
         {
            Scores[i][j] = MA_getAverageValue(ma_info, psd->PaVars[i][j]);
            if( Scores[i][j] != 0 )
               selected[i][j] = TRUE;
            else
               selected[i][j] = FALSE;
         }
      SCIP_CALL( printToFile(scip, psd, Scores, selected, filename, format, MA_getTotalAveragesTime(ma_info), MA_getTotalAveragesScore(ma_info), NULL) );

      /* Deallocate memory */
      for( i = 0; i < psd->n; i++ )
      {
         SCIPfreeMemoryArray(scip, &(Scores[i]));
         SCIPfreeMemoryArray(scip, &(selected[i]));
      }
      SCIPfreeMemoryArray(scip, &Scores);
      SCIPfreeMemoryArray(scip, &selected);
   }

   return SCIP_OKAY;
}

/** Finds the location in a string of the file extension.
 *
 *  @param filename The filename to inspect.
 *  @return The location of the "." marking the beginning of the file extension, or -1 if there is no extension.
 */
static int findExtension(char* filename)
{
   int i;
   for( i = strlen(filename); i >= 0; i-- )
      if( filename[i] == '.' )
         return i;
      else if( filename[i] == '/' )
         return -1;
   return -1;
}

/** Constructs a filename for output for a particular iteration of the program.
 *
 *  The filename constructed will be the input filename the iteration number appearing before the file extension
 *  or appended if there is no extension.  If there is only a single network to find, then the iteration number is
 *  not inserted.  For special values ("" and "stdout"), the value returned is the same value as given as input.
 *  Any occurrence of @c <probname> will be replaced with the basename of the problem.
 *
 *  @param filename The filename to insert the string in to.
 *  @param probname The basename of the problem.
 *  @param nbns The number of Bayesian networks that the program is trying to find.
 *  @param iteration The current iteration of the program.
 *  @return The filename that should be used for output.
 */
static
char* createFilename(
   char* filename,
   const char* probname,
   int nbns,
   int iteration
   )
{
   if( strcmp(filename, "") == 0 || strcmp(filename, "stdout") == 0 )
      return strdup(filename);
   else
   {
      /* Replace all occurances of <basename> */
      char* replaced_filename;
      if( strstr(filename, "<probname>") == NULL )
      {
         replaced_filename = strdup(filename);
      }
      else
      {
         char* ans;
         int i;
         int to_pos = 0;
         int from_pos = 0;
         char* replace_pos = filename;
         int num_occurances = 0;
         char* tmp = strstr(filename, "<probname>");
         while( tmp != NULL )
         {
            tmp = strstr(tmp + 1, "<probname>");
            num_occurances++;
         }
         ans = malloc((strlen(filename) + num_occurances * (strlen(probname) - 10) + 1) * sizeof(char));
         for( i = 0; i < num_occurances; i++ )
         {
            replace_pos = strstr(&(filename[from_pos]), "<probname>");
            while( &(filename[from_pos]) != replace_pos )
            {
               ans[to_pos] = filename[from_pos];
               from_pos++;
               to_pos++;
            }
            strcpy(&(ans[to_pos]), probname);
            to_pos += strlen(probname);
            from_pos += strlen("<probname>");
         }
         strcpy(&(ans[to_pos]), &(filename[from_pos]));
         replaced_filename = ans;
      }

      /* Add the interation numbers */
      if( nbns == 1 )
         /* If only one BN, there is no need to add iteration numbers */
         return replaced_filename;
      else
      {
         /* Need to add iteration numbers for each BN */
         int i;
         int extpos;
         char* ans;
         char insertion[SCIP_MAXSTRLEN];
         sprintf(insertion, "_%d", iteration + 1);
         ans = malloc((strlen(replaced_filename) + strlen(insertion) + 1) * sizeof(char));
         extpos = findExtension(replaced_filename);
         if( extpos == -1 )
         {
            for( i = 0; i < (int)strlen(replaced_filename); i++ )
               ans[i] = replaced_filename[i];
            for( i = 0; i < (int)strlen(insertion); i++ )
               ans[i + strlen(replaced_filename)] = insertion[i];
            ans[strlen(replaced_filename) + strlen(insertion)] = '\0';
         }
         else
         {
            for( i = 0; i < extpos; i++ )
               ans[i] = replaced_filename[i];
            for( i = 0; i < (int)strlen(insertion); i++ )
               ans[i + extpos] = insertion[i];
            for( i = extpos; i < (int)strlen(replaced_filename); i++ )
               ans[i + strlen(insertion)] = replaced_filename[i];
            ans[strlen(replaced_filename) + strlen(insertion)] = '\0';
         }
         return ans;
      }
   }
}

/** Prints the current problem in cip format.
 *
 *  @param scip The SCIP instance for which to print the problem.
 *  @param run The iteration of the main loop that the problem is to be solved on.
 *  @return SCIP_OKAY if the problem was printed correctly or an appropriate error message otherwise.
 */
SCIP_RETCODE IO_printProblem(
   SCIP* scip,
   int run
   )
{
   int nbns;
   char* probfile;
   char* outputfile;
   SCIP_Bool presolveprob;

   SCIP* targetscip;
   SCIP_CONS* cons;
   SCIP_Bool valid;
   SCIP_Bool printmeta;

   SCIP_CALL( SCIPgetBoolParam(scip,"gobnilp/presolveprob",&presolveprob) );

   SCIPgetIntParam(scip, "gobnilp/nbns", &nbns);
   SCIPgetStringParam(scip, "gobnilp/outputfile/problem", &probfile);
   SCIPgetBoolParam(scip, "gobnilp/printmeta", &printmeta);
   outputfile = createFilename(probfile, SCIPgetProbName(scip), nbns, run);

   if( strcmp(outputfile, "") == 0 )
   {
      free(outputfile);
      return SCIP_OKAY;
   }
   else
   {
      if( strcmp(outputfile, "stdout") == 0 )
         outputfile = NULL;
      else
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Writing problem to %s\n", outputfile);

      if( !printmeta )
      {
         SCIP_CALL( SCIPcreate(&targetscip) );
         SCIP_CALL( SCIPcopyOrig(scip, targetscip, NULL, NULL, "nometa", FALSE, FALSE, &valid) );
         /* if(!valid) */
         /*    return SCIP_ERROR; */
         cons = SCIPfindCons(targetscip, "parent set data");
         if( cons != NULL )
            SCIP_CALL( SCIPdelCons(targetscip, cons) );
         cons = SCIPfindCons(targetscip, "pedigree data");
         if( cons != NULL )
            SCIP_CALL( SCIPdelCons(targetscip, cons) );
         cons = SCIPfindCons(targetscip, "property data");
         if( cons != NULL )
            SCIP_CALL( SCIPdelCons(targetscip, cons) );
      }
      else
      {
         targetscip = scip;
      }

      if( presolveprob )
      {
         SCIPpresolve(targetscip);
         SCIP_CALL( SCIPwriteTransProblem(targetscip, outputfile, NULL, FALSE) );
      }
      else
         SCIP_CALL( SCIPwriteOrigProblem(targetscip, outputfile, NULL, FALSE) );
      free(outputfile);

      if( !printmeta )
         SCIP_CALL( SCIPfree(&targetscip) );

      return SCIP_OKAY;
   }

   free(outputfile);
}

/** Prints appropriate information about each optimal solution obtained.
 *
 *  @param scip The SCIP instance for which the solution has been found.
 *  @param run The iteration of the main loop that the solution was found on.
 *  @param psd The parentage data for the problem.
 *  @return SCIP_OKAY if printing succeeded or an appropriate error code otherwise.
 */
SCIP_RETCODE IO_doIterativePrint(
   SCIP* scip,
   MA_info* ma_info,
   ParentSetData* psd,
   int run
   )
{
   char* solfile;
   char* dotfile;
   char* pedfile;
   char* satfile;
   char* matfile;
   char* mecfile;

   char* run_solfile;
   char* run_dotfile;
   char* run_pedfile;
   char* run_satfile;
   char* run_matfile;
   char* run_mecfile;


   int avgoutputoffset;
   int avgoutputstep;
   char* avgsolfile;
   char* avgdotfile;
   char* avgpedfile;
   char* avgsatfile;
   char* avgmatfile;

   char* run_avgsolfile;
   char* run_avgdotfile;
   char* run_avgpedfile;
   char* run_avgsatfile;
   char* run_avgmatfile;


   int nbns;
   SCIPgetIntParam(scip, "gobnilp/nbns", &nbns);

   SCIPgetStringParam(scip, "gobnilp/outputfile/solution", &solfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/dot", &dotfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/pedigree", &pedfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/scoreandtime", &satfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/adjacencymatrix", &matfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/mec", &mecfile);

   SCIPgetIntParam(scip, "gobnilp/avgoutputoffset", &avgoutputoffset);
   SCIPgetIntParam(scip, "gobnilp/avgoutputstep", &avgoutputstep);
   SCIPgetStringParam(scip, "gobnilp/outputfile/solutionavg", &avgsolfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/dotavg", &avgdotfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/pedigreeavg", &avgpedfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/scoreandtimeavg", &avgsatfile);
   SCIPgetStringParam(scip, "gobnilp/outputfile/adjacencymatrixavg", &avgmatfile);

   if( SCIPgetVerbLevel(scip) >= 5 )
   {
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );
      SCIP_CALL( SCIPprintBranchingStatistics(scip, NULL) );
   }

   SCIP_CALL( printSCIPSolution(scip) );

   run_solfile = createFilename(solfile, SCIPgetProbName(scip), nbns, run);
   run_dotfile = createFilename(dotfile, SCIPgetProbName(scip), nbns, run);
   run_pedfile = createFilename(pedfile, SCIPgetProbName(scip), nbns, run);
   run_satfile = createFilename(satfile, SCIPgetProbName(scip), nbns, run);
   run_matfile = createFilename(matfile, SCIPgetProbName(scip), nbns, run);
   run_mecfile = createFilename(mecfile, SCIPgetProbName(scip), nbns, run);

   SCIP_CALL( printSolution(scip, psd, run_solfile, (char*)"legacy") );
   SCIP_CALL( printSolution(scip, psd, run_dotfile, (char*)"dot") );
   SCIP_CALL( printSolution(scip, psd, run_pedfile, (char*)"pedigree") );
   SCIP_CALL( printSolution(scip, psd, run_satfile, (char*)"scoreandtime") );
   SCIP_CALL( printSolution(scip, psd, run_matfile, (char*)"adjacencymatrix") );
   SCIP_CALL( printSolution(scip, psd, run_mecfile, (char*)"mec") );

   free(run_solfile);
   free(run_dotfile);
   free(run_pedfile);
   free(run_satfile);
   free(run_matfile);
   free(run_mecfile);

   if( avgoutputoffset > -1 && (run + 1) >= avgoutputoffset )
   {
      if( (run - avgoutputoffset + 1) % avgoutputstep == 0 )
      {
         run_avgsolfile = createFilename(avgsolfile, SCIPgetProbName(scip), nbns, run);
         run_avgdotfile = createFilename(avgdotfile, SCIPgetProbName(scip), nbns, run);
         run_avgpedfile = createFilename(avgpedfile, SCIPgetProbName(scip), nbns, run);
         run_avgsatfile = createFilename(avgsatfile, SCIPgetProbName(scip), nbns, run);
         run_avgmatfile = createFilename(avgmatfile, SCIPgetProbName(scip), nbns, run);

         SCIP_CALL( printAverages(scip, ma_info, psd, run_avgsolfile, (char*)"legacy") );
         SCIP_CALL( printAverages(scip, ma_info, psd, run_avgdotfile, (char*)"dot") );
         SCIP_CALL( printAverages(scip, ma_info, psd, run_avgpedfile, (char*)"pedigree") );
         SCIP_CALL( printAverages(scip, ma_info, psd, run_avgsatfile, (char*)"scoreandtime") );
         SCIP_CALL( printAverages(scip, ma_info, psd, run_avgmatfile, (char*)"adjacencymatrix") );

         free(run_avgsolfile);
         free(run_avgdotfile);
         free(run_avgpedfile);
         free(run_avgsatfile);
         free(run_avgmatfile);

      }
   }

   if( strcmp(solfile, "") != 0 || strcmp(dotfile, "") != 0 || strcmp(pedfile, "") != 0 || strcmp(satfile, "") != 0 || strcmp(matfile, "") != 0 || strcmp(mecfile, "") != 0 )
      SCIPinfoMessage(scip, NULL, "\n");
   else if((avgoutputoffset > -1 && (run + 1) >= avgoutputoffset) && ((run - avgoutputoffset + 1) % avgoutputstep == 0) && (
              strcmp(avgsolfile, "") != 0 || strcmp(avgdotfile, "") != 0 || strcmp(avgpedfile, "") != 0 || strcmp(avgsatfile, "") != 0 || strcmp(avgmatfile, "") != 0))
      SCIPinfoMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}

/** Prints any of the current SCIP or GOBNILP parameters not at their default value.
 *
 *  @param scip The SCIP instance to consult the parameters of.
 *  @return SCIP_OKAY if the parameters were printed correctly, or an error code otherwise.
 */
SCIP_RETCODE IO_printParameters(SCIP* scip)
{
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "START Parameters not at default value\n");
   /* verblevel 4 is SCIP_VERBLEVEL_FULL */
   if( SCIPgetVerbLevel(scip) >= 5 )
      SCIP_CALL( SCIPwriteParams(scip, NULL, FALSE, TRUE) );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "END Parameters not at default value\n");
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "\n");
   fflush(stdout);

   return SCIP_OKAY;
}

/** Prints a header which describes the GOBNILP and SCIP systems being used.
 *
 *  @param scip The SCIP instance that is being used.
 *  @return SCIP_OKAY if printing succeeded or an error code otherwise.
 */
SCIP_RETCODE IO_printHeader(SCIP* scip)
{
   /* output version information */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "GOBNILP version %s [GitHash: %s ]\n", GOBNILP_VERSION, GOBNILP_GITHASH);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Solving the BN structure learning problem using SCIP.\n\n");

   /* verblevel 3 is SCIP_VERBLEVEL_MINIMAL */
   if( SCIPgetVerbLevel(scip) >= 2 )
   {
      SCIPprintVersion(scip, NULL);
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "\n");
   return SCIP_OKAY;
}

/** Conditionally, prints the scores to file.
 *
 *
 *  @param scip The SCIP instance the scores belong to.
 *  @param psd The parent set data to print the data for.
 *  @param scores Local scores (only used if SCIP variables have not been created)
 *  @return SCIP_OKAY if printing succeeded, or an error code otherwsie.
 */
SCIP_RETCODE IO_printScoresInJKLFormat(
   SCIP* scip, 
   ParentSetData* psd,
   SCIP_Real** scores
   )
{
   char* filename;
   FILE* file;
   int i;
   int j;
   int k;

   assert(scip != NULL);
   assert(psd != NULL);

   SCIPgetStringParam(scip, "gobnilp/outputfile/scores", &filename);
   if( scores == NULL )
      filename = createFilename(filename, SCIPgetProbName(scip), 1, 1);
   else
      filename = createFilename(filename, "", 1, 1);

   if( strcmp(filename, "") == 0 )
    {
      free(filename);
      return SCIP_OKAY;
   }
   else if( strcmp(filename, "stdout") == 0 )
      file = stdout;
   else
   {
      file = fopen(filename, "w");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Writing scores to %s\n", filename);
   }

   fprintf(file, "%d\n", psd->n);
   assert(psd->nodeNames != NULL);
   assert(psd->nParentSets != NULL);
   assert(psd->PaVars == NULL || scores == NULL);
   for( i =  0; i < psd->n; i++ )
   {
      assert(psd->nParents[i] != NULL);
      fprintf(file, "%s %d\n", psd->nodeNames[i], psd->nParentSets[i]);
      for( j = 0; j < psd->nParentSets[i]; j++ )
      {
         assert(psd->ParentSets[i][j] != NULL);
         if( psd->PaVars == NULL )
            fprintf(file, "%.8f %d", scores[i][j], psd->nParents[i][j]);
         else
            fprintf(file, "%.8f %d", SCIPvarGetObj(psd->PaVars[i][j]), psd->nParents[i][j]);
         for( k = 0; k < psd->nParents[i][j]; k++ )
            fprintf(file, " %s", psd->nodeNames[psd->ParentSets[i][j][k]]);
         fprintf(file, "\n");
      }
   }

   if( file != stdout )
      fclose(file);

   free(filename);

   return SCIP_OKAY;
}


static SCIP_Bool isRecognisedMETA(const char* name)
{
   if((strcmp(name, "pss_version") == 0) ||
         (strcmp(name, "time_created") == 0) ||
         (strcmp(name, "scorer_name") == 0) ||
         (strcmp(name, "scorer_version") == 0) ||
         (strcmp(name, "scorer_url") == 0) ||
         (strcmp(name, "score_type") == 0) ||
         (strcmp(name, "ess") == 0) ||
         (strcmp(name, "parent_limit") == 0) ||
         (strcmp(name, "pruning") == 0) ||
         (strcmp(name, "input_file") == 0) ||
         (strcmp(name, "num_records") == 0))
      return TRUE;
   return FALSE;
}
static void printCommentLine(FILE* file, const char* comment)
{
   fprintf(file, "# %s\n", comment);
}
static void printBlankLine(FILE* file)
{
   fprintf(file, "\n");
}
static void printVariableLine(FILE* file, const char* name)
{
   fprintf(file, "VAR %s\n", name);
}
static void printMetaLine(FILE* file, const char* name, const char* value)
{
   fprintf(file, "META %s = %s\n", name, value);
}
static void lookupAndPrintMetaLine(SCIP* scip, FILE* file, PropertyData* prop, int individual, const char* name)
{
   char* value = NULL;
   if( individual == -1 && PR_hasGlobalProperty(scip, prop, name) )
      value = PR_getGlobalProperty(scip, prop, name);
   else if( individual != -1 && PR_hasProperty(scip, prop, individual, name) )
      value = PR_getProperty(scip, prop, individual, name);

   if( value != NULL )
      fprintf(file, "META %s = %s\n", name, value);
}
static void printScoreLine(FILE* file, SCIP_VAR* var, int num_parents, int* parents, char** names)
{
   int i;
   fprintf(file, "%f", SCIPvarGetObj(var));
   for( i = 0; i < num_parents; i++ )
      fprintf(file, " %s", names[parents[i]]);
   fprintf(file, "\n");
}
SCIP_RETCODE IO_printScoresInPSSFormat(SCIP* scip, ParentSetData* psd, PropertyData* prop)
{
   char* filename;
   FILE* file;
   int i, j;
   time_t time_now = time(NULL);
   char* current_time = asctime(gmtime(&time_now));
   current_time[strlen(current_time) - 1] = '\0';

   SCIPgetStringParam(scip, "gobnilp/outputfile/pss", &filename);
   filename = createFilename(filename, SCIPgetProbName(scip), 1, 1);
   if( strcmp(filename, "") == 0 )
   {
      free(filename);
      return SCIP_OKAY;
   }
   else if( strcmp(filename, "stdout") == 0 )
      file = stdout;
   else
   {
      file = fopen(filename, "w");
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Writing scores to %s\n", filename);
   }

   printMetaLine(file, "pss_version", "0.1");
   printMetaLine(file, "time_created", current_time);
   printBlankLine(file);

   printCommentLine(file, "Info on scoring program");
   lookupAndPrintMetaLine(scip, file, prop, -1, "scorer_name");
   lookupAndPrintMetaLine(scip, file, prop, -1, "scorer_version");
   lookupAndPrintMetaLine(scip, file, prop, -1, "scorer_url");
   printBlankLine(file);

   printCommentLine(file, "Info on scoring algorithm");
   lookupAndPrintMetaLine(scip, file, prop, -1, "score_type");
   lookupAndPrintMetaLine(scip, file, prop, -1, "ess");
   lookupAndPrintMetaLine(scip, file, prop, -1, "parent_limit");
   lookupAndPrintMetaLine(scip, file, prop, -1, "pruning");
   printBlankLine(file);

   printCommentLine(file, "Info on data set");
   lookupAndPrintMetaLine(scip, file, prop, -1, "num_records");
   lookupAndPrintMetaLine(scip, file, prop, -1, "input_file");
   printBlankLine(file);

   printCommentLine(file, "Unrecognised global METAs from input");
   for( i = 0; i < prop->num_global; i++ )
      if( !isRecognisedMETA(prop->global_property_names[i]) )
         printMetaLine(file, prop->global_property_names[i], prop->global_property_values[i]);
   printBlankLine(file);

   for( i =  0; i < psd->n; i++ )
   {
      printVariableLine(file, psd->nodeNames[i]);
      for( j = 0; j < prop->num_properties[i]; j++ )
         printMetaLine(file, prop->property_names[i][j], prop->property_values[i][j]);
      for( j = 0; j < psd->nParentSets[i]; j++ )
         printScoreLine(file, psd->PaVars[i][j], psd->nParents[i][j], psd->ParentSets[i][j], psd->nodeNames);
      printBlankLine(file);
   }
   printBlankLine(file);

   if( file != stdout )
      fclose(file);

   free(filename);

   return SCIP_OKAY;
}




/** Prints solutions found using countsols
 *
 *  @param scip The SCIP instance to print solutions for.
 *  @param filename The name of the file where solutions will be put.
 *  @return SCIP_OKAY if printing succeeded or an error otherwise.
 */
SCIP_RETCODE IO_printcountsols(
   SCIP* scip,
   ParentSetData* psd,
   char* filename
   )
{
   SCIP_VAR** activevars; /* active variables */
   int nactivevars;
   SCIP_SPARSESOL** sparsesols;
   int nsparsesols; 
   FILE* file;
   SCIP_VAR** vars;
   SCIP_Real* scalars;
   SCIP_Longint* sol;
   int s;
   int i;
   int j;
   int k;
   SCIP_VAR* transvar;
   SCIP_Bool** selected;
   SCIP_Real** scores;
   SCIP_Real total_score;
   
   SCIP_Bool valid;
   SCIP_Longint nsols;
   SCIP_SPARSESOL* sparsesol;
   /* SCIP_RETCODE retcode; */
   /* SCIP_VAR* transvar; */

   /* get number of integral variables */
   const int nallvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   if( strcmp(filename, "") == 0 )
   {
      SCIPerrorMessage("Cannot use empty filename for collecting solutions\n");
      return SCIP_ERROR;
   }      

   file = fopen(filename, "w");
   if( file == NULL )
   {
      SCIPerrorMessage("error creating file <%s>\n", filename);
      return SCIP_ERROR;
   }

   valid = FALSE;
   nsols = SCIPgetNCountedSols(scip, &valid);
   if( !valid )
   {
      /* don't bother using method of cons_countsols to print out number of solutions */
      SCIPerrorMessage("Too many solutions (to be counted by a long int)\n");
      return SCIP_ERROR;
   }

   if( nsols == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "There are no counted solutions.\n");
      return SCIP_OKAY;
   }

   SCIPgetCountedSparseSols(scip, &activevars, &nactivevars, &sparsesols, &nsparsesols);

   if( nsparsesols == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No solutions were collected.\n");
      return SCIP_OKAY;
   }

   /* get memory to store active solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &sol, nactivevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nactivevars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scalars, nactivevars) );

   SCIP_CALL( SCIPallocBufferArray(scip, &selected, psd->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(selected[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(scores[i]), psd->nParentSets[i]) );
      for( k = 0; k < psd->nParentSets[i]; ++k )
         scores[i][k] = SCIPvarGetObj(psd->PaVars[i][k]);	
   }
   
   /* loop over all sparse solutions */
   for( s = 0; s < nsparsesols; ++s )
   {
      
      sparsesol = sparsesols[s]; /*lint !e613*/
      assert(sparsesol != NULL);
      assert(SCIPsparseSolGetNVars(sparsesol) == nactivevars);


      /* get first solution of the sparse solution */
      SCIPsparseSolGetFirstSol(sparsesol, sol, nactivevars);

      do
      {
         total_score = 0.0;
         for( i = 0; i < psd->n; ++i )
         {
            for( k = 0; k < psd->nParentSets[i]; ++k )
               selected[i][k] = FALSE;
            for( k = 0; k < psd->nParentSets[i]; ++k )
            {
               /* compute value for PaVars[i][k] in this solution
                  and mark as selected if > 0.5 (i.e. set to 1 )
               */

               SCIP_Real constant;
               SCIP_Real realvalue;
               int requiredsize;
               int nvars;
               int idx;

               SCIP_CALL( SCIPgetTransformedVar(scip, psd->PaVars[i][k], &transvar) );

               /* get representation of transvar as a linear function of the active vars */

               vars[0] = transvar;
               scalars[0] = 1.0;
               nvars = 1;
               constant = 0.0;

               SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, scalars, &nvars, nallvars, &constant, &requiredsize, TRUE) );
               assert(requiredsize <= nallvars);
               assert(nvars <= nactivevars);

               /* compute objective value of transvar from objective values of active variables 
                  defining it (which are now in vars array) */

               realvalue = constant;
               for( j = 0; j < nvars; ++j )
               {
                  /* no hashmap so have to scan */
                  for( idx = 0; idx < nactivevars; idx++ )
                  {
                     if( activevars[idx] == vars[j] )
                     {
                        realvalue += scalars[j] * sol[idx];
                        break;
                     }
                  }
               }
               
               /* mark kth parent set as selected if obj value in this solution above 0.5 */ 
               if( realvalue > 0.5 )
               {
                     selected[i][k] = TRUE;
                     total_score += scores[i][k];
                     break;
               }
            }
         }
         SCIP_CALL( printSolutionLegacyFormat(scip, psd, scores, selected, total_score, file) );
         fprintf(file,"\n\n");
      }
      while( SCIPsparseSolGetNextSol(sparsesol, sol, nactivevars) );
   }

   fprintf(file,"%d\n",nsparsesols);
   fclose(file);

   SCIPfreeBufferArray(scip, &sol);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &scalars);

   for( i = 0; i < psd->n; ++i )
   {
      SCIPfreeBufferArray(scip, &(selected[i]));
      SCIPfreeBufferArray(scip, &(scores[i]));
   }
   SCIPfreeBufferArray(scip, &selected);
   SCIPfreeBufferArray(scip, &scores);
   
   return SCIP_OKAY;
}
