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
 *  Provides the core functionality for the Bayesian network learning problem.
 */

/*#define SCIP_DEBUG*/

/*
This file was created by editing the file psd_lop.c that comes with the linear ordering example
in SCIP
*/
#include <string.h>
#include <ctype.h>
#include <float.h>
#include "probdata_bn.h"
#include "pedigrees.h"
#include "pedigree_scorer.h"
#include "cons_dagcluster.h"
#include "cons_ci.h"
#include "cons_chordal.h"
#include "cons_distinguishedrep.h"
#include "heur_sinks.h"
#include "utils.h"
#include "output.h"
#include "metadata.h"
#include "scoring.h"
#include "scip/misc.h"
#include "disp_clearcols.h"
#include "event_splitdag.h"
#include "cons_lop.h"
#include "cons_partialordering.h"
#include "cons_vanilla.h"

#define DEFAULT_GOBNILP_INPUT_FORMAT "jkl"


/** delete problem data */
static
SCIP_DECL_PROBDELORIG(probdelorigBN)
{  /*lint --e{831} */

   int i;
   int j;

   assert( probdata != NULL );
   assert( *probdata != NULL );

   assert( (*probdata)->psd != NULL );
   assert( (*probdata)->scores != NULL );

   if( (*probdata)->ancestorvars != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
      {
         for( j = 0; j < (*probdata)->psd->n; j++ )
            if( i != j )
               SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->ancestorvars)[i][j])) );
         SCIPfreeMemoryArray(scip, &((*probdata)->ancestorvars[i]));
      }
      SCIPfreeMemoryArray(scip, &((*probdata)->ancestorvars));
   }

   if( (*probdata)->totalordervars != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
      {
         for( j = 0; j < (*probdata)->psd->n; j++ )
            if( i != j )
               SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->totalordervars)[i][j])) );
         SCIPfreeMemoryArray(scip, &((*probdata)->totalordervars[i]));
      }
      SCIPfreeMemoryArray(scip, &((*probdata)->totalordervars));
   }

   if( (*probdata)->imsetvars != NULL )
   {
      int i0;
      int i1;
      int i2;
      
      for( i0 = 0; i0 < (*probdata)->psd->n; i0++ )
      {
         for( i1 = i0 + 1; i1 < (*probdata)->psd->n; i1++ )
         {
            int ii1 = i1 - i0 - 1;
            for( i2 = i1 + 1; i2 < (*probdata)->psd->n; i2++ )
            {
               int ii2 = i2 - i1 - 1;
               if( (*probdata)->imsetvars[i0][ii1][ii2] != NULL ) 
                  SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->imsetvars)[i0][ii1][ii2])) );
            }
            SCIPfreeMemoryArray(scip, &((*probdata)->imsetvars[i0][ii1]));
         }
         SCIPfreeMemoryArray(scip, &((*probdata)->imsetvars[i0]));
      }
      SCIPfreeMemoryArray(scip, &((*probdata)->imsetvars));
   }


   
   if( (*probdata)->posind != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
      {
         for( j = 0; j < (*probdata)->psd->n; j++ )
            SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->posind)[i][j])) );
         SCIPfreeMemoryArray(scip, &((*probdata)->posind[i]));
      }
      SCIPfreeMemoryArray(scip, &((*probdata)->posind));
   }

   
   if( (*probdata)->posvars != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
         SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->posvars)[i])) );
      SCIPfreeMemoryArray(scip, &((*probdata)->posvars));
   }

   if( (*probdata)->ispa != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
         SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->ispa)[i])) );
      SCIPfreeMemoryArray(scip, &((*probdata)->ispa));
   }

   if( (*probdata)->nch != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
         SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->nch)[i])) );
      SCIPfreeMemoryArray(scip, &((*probdata)->nch));
   }

   
   if( (*probdata)->pasize != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
         SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->pasize)[i])) );
      SCIPfreeMemoryArray(scip, &((*probdata)->pasize));
   }

   
   if( (*probdata)->isFemale != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
         SCIP_CALL( SCIPreleaseVar(scip, &(((*probdata)->isFemale)[i])) );
      SCIPfreeMemoryArray(scip, &((*probdata)->isFemale));
   }

   
   if( (*probdata)->im_vars != NULL )
   {
      for( i = 0; i < (*probdata)->psd->n; i++ )
      {
         for( j = 0; j < (*probdata)->psd->n; j++ )
            SCIPfreeMemoryArray(scip, &((*probdata)->im_vars[i][j]));
         SCIPfreeMemoryArray(scip, &((*probdata)->im_vars[i]));
      }
      SCIPfreeMemoryArray(scip, &((*probdata)->im_vars));
   }

   for( i = 0; i < (*probdata)->psd->n; i++ )
      SCIPfreeMemoryArray(scip, &((*probdata)->scores[i]));
   SCIPfreeMemoryArray(scip, &((*probdata)->scores));
   SCIP_CALL( PS_deallocateParentSetData(scip, &((*probdata)->psd), TRUE) );
   
   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** Includes all plugins needed by the problem.
 *
 *  @param scip The SCIP instance to add the plugins to.
 *  @return SCIP_OKAY if all plugins were added successfully, or an error otherwise.
 */
SCIP_RETCODE BN_includePlugins(
   SCIP* scip
   )
{
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( DC_includeConshdlr(scip) );
   SCIP_CALL( SCIPincludeConshdlrLOP(scip) );
   SCIP_CALL( SCIPincludeConshdlrPartialOrdering(scip) );
   /*SCIP_CALL(  SCIPincludeConshdlrIndicator(scip)  ); */
   SCIP_CALL( SCIPincludeConshdlrCi(scip) );
   SCIP_CALL( SCIPincludeConshdlrChordal(scip) );
   SCIP_CALL( HS_includePrimal(scip) );
   SCIP_CALL( MD_initialiseMetadata(scip) );
   SCIP_CALL( SCIPincludeDispClearCols(scip) );
   SCIP_CALL( SCIPincludeEventHdlrSplitDAG(scip) );
   return SCIP_OKAY;
}



/** Adds GOBNILP specific parameters to those recognised by SCIP.
 *
 *  @param scip The SCIP instance to add to the parameters to.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error coede otherwise.
 */
SCIP_RETCODE BN_addParameters(
   SCIP* scip
   )
{
   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/verblevelsetscols",
         "whether a verbosity level of at most 3 suppresses columns in the display",
         TRUE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/vanilla",
         "whether the problem is 'vanilla' BN learning",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/onlyscores",
         "whether to only compute local scores without creating SCIP variables",
         FALSE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/noimmoralities",
         "whether to disallow immoralities",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/decomposable",
         "whether to assume a decomposable model is being learned with a likelihood-equivalent score",
         FALSE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/orderedcoveredarcs",
         "whether to only allow a covered arc i<-j if i<j",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/printscipsol",
         "whether to (additionally) print BNs in SCIP solution format",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/presolveprob",
         "whether to presolve the problem before printing it out CIP",
         TRUE
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/maxnchildren",
         "the maximum number of children a node can have (-1 for no limit)",
         -1, -1, INT_MAX
         ) );



   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/nbns",
         "gobnilp to find the 'nbns' best BNs (in decreasing order of score )",
         1, 1, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/minfounders",
         "minimum number of founders",
         0, 0, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/maxfounders",
         "maximum number of founders (-1 for no upper bound )",
         -1, -1, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/minparents",
         "minimum number of parents",
         0, 0, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/maxparents",
         "maximum number of parents (-1 for no upper bound )",
         -1, -1, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/minedges",
         "minimum number of edges",
         0, 0, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/maxedges",
         "maximum number of edges (-1 for no upper bound )",
         -1, -1, INT_MAX
         ) );

   SCIP_CALL( UT_addStringParam(scip,
         "gobnilp/dagconstraintsfile",
         "file containing constraints on dag structure",
         ""
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/mergedelimiters",
         "whether to treat consecutive delimitors in the input file as one",
         TRUE
         ) );

   SCIP_CALL( UT_addStringParam(scip,
         "gobnilp/delimiter",
         "the delimiter for fields in the input file (special values - whitespace, tab)",
         "whitespace"
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/totalordervars",
         "whether to create totalorder variables",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/ancestorvars",
         "whether to create ancestor variables",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/looseancestors",
         "whether ancestor variables encode an arbitrary consistent partial order rather than the ancestor relation",
         FALSE
         ) );

   
   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/nchildrenvars",
         "whether to create variables counting the number of children a node has",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/bestparentsfororder",
         "whether to only allow best scoring parents for a given total order",
         TRUE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/useconslinearordering",
         "whether to use SCIP's linear ordering constraint (when total order vars exist)",
         TRUE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/useconspartialordering",
         "whether to use the partial ordering constraint (when ancestor vars exist)",
         TRUE
         ) );
   
   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/useconsdagcluster",
         "whether to use the dagcluster acyclicity constraint",
         TRUE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/inittransitivity",
         "whether to add transitivity constraints on totalorder variables initially (rather than later as cutting planes)",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/imsetvars",
         "whether to create imset variables",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/posindvars",
         "whether to create position indicator variables",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/posvars",
         "whether to create position variables",
         FALSE
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/posvarpriority",
         "branching priority of posind variables",
         0, -INT_MAX, INT_MAX
         ) );


   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/maxpos",
         "maximum value of a position variable (-1 for no restriction)",
         -1, -1, INT_MAX
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/postotal",
         "whether to create position variables indicate position in a total order of BN variables",
         FALSE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/pasizevars",
         "whether to create parent set size variables",
         FALSE
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/pasizepriority",
         "branching priority of parent set size variables",
         0, 0, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/pavarspriority",
         "branching priority of family variables",
         0, -INT_MAX, INT_MAX
         ));

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/ancestorvarpriority",
         "branching priority of edge variables",
         0, -INT_MAX, INT_MAX
         ));

   
   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/edgevarpriority",
         "branching priority of edge variables",
         10, -INT_MAX, INT_MAX
         ));

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/nchildrenvarpriority",
         "branching priority of nchildrenvar variables",
         0, -INT_MAX, INT_MAX
         ));


   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/arrowvarpriority",
         "branching priority of arrow variables",
         10, -INT_MAX, INT_MAX
         ));


   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/imsetvarpriority",
         "branching priority of imset variables",
         0, -INT_MAX, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/totalordervarpriority",
         "branching priority of total order variables",
         0, -INT_MAX, INT_MAX
         ) );

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/posindvarpriority",
         "branching priority of posind variables",
         0, -INT_MAX, INT_MAX
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/splitdags",
         "whether to split dagcluster constraints into strongly connected components",
         TRUE
         ) );

   SCIP_CALL( UT_addRealParam(scip,
         "gobnilp/edge_penalty",
         "edge_penalty*|edges| is subtracted from objective",
         0, SCIP_REAL_MIN, SCIP_REAL_MAX
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/kbestMEC",
         "whether GOBNILP should return exemplars of the k best MECs instead of the k best BNs",
         FALSE
         ) );

   SCIP_CALL( UT_addRealParam(scip,
         "gobnilp/objlimit",
         "sets objective limit, only BNs with a strictly better score are considered",
         SCIP_REAL_MIN, SCIP_REAL_MIN, SCIP_REAL_MAX
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/countsols",
         "whether to count (and perhaps collect) feasible solution rather than search for an optimal solution",
         FALSE
         ) );

   SCIP_CALL( UT_addLongintParam(scip,
         "gobnilp/countsols/sollimit",
         "counting stops, if the given number of solutions were found (-1: no limit)",
         -1LL, -1LL, SCIP_LONGINT_MAX
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/countsols/collect",
         "whether to collect feasible solutions (as opposed to just counting them)",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/distinguishedrep",
         "whether to only allow DAGs which are the distinguished representative of their Markov equivalence class",
         FALSE
         ));

   SCIP_CALL( UT_addIntParam(scip,
         "gobnilp/arrowfamilylink",
         "how arrow variables and family variables are linked",
         4, 1, 4
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/allowmultiaggregation",
         "whether to allow family, arrow and edge variables to be multi-aggregated",
         FALSE
         ));

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/gobnilpparams",
         "whether to set SCIP parameters to values optimised for GOBNILP (else use defaults)",
         TRUE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/ancestralbranching",
         "whether to impose 'ancestral branching'",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/arrowvarsinitial",
         "should variables representing arrows be in the initial root LP?",
         TRUE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/arrowvarsremovable",
         "should variables representing arrows be removable?",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/edgevarsinitial",
         "should variables representing undirected edges be in the initial root LP?",
         TRUE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/edgevarsremovable",
         "should variables representing undirected edges be removable?",
         FALSE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/pavarsinitial",
         "should family variables be in the initial root LP?",
         TRUE
         ) );
   
   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/pavarsremovable",
         "should family variables be removable?",
         FALSE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/arrowconsinitial",
         "should the constraint linking arrow to family variables be in the initial root LP?",
         TRUE
         ) );

   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/edgeconsinitial",
         "should the constraint linking edge to arrow variables be in the initial root LP?",
         TRUE
         ) );


   SCIP_CALL( UT_addBoolParam(scip,
         "gobnilp/convexconsinitial",
         "should the constraint that each BN node has exactly one parent be in the initial root LP?",
         TRUE
         ) );


   
   SCIP_CALL( IO_addOutputParameters(scip) );
   SCIP_CALL( PD_addPedigreeParameters(scip) );
   SCIP_CALL( SC_addScoringParameters(scip) );

   return SCIP_OKAY;
}

/** Gets the format in which the input file has been given.
 *
 * If the format has been set on the command line, this will be given.
 * If this is not set, then the suffix of the filename will be examined.
 * If this does not work, then DEFAULT_GOBNILP_PARAMS_FILE is used.
 *
 * @param scip The SCIP instance this applies to.
 * @param filename The name of the file to read.
 * @return The format of the input file.
 */
static
char* getInputFormat(
   SCIP*       scip,
   char*   inputformat,
   const char* filename
   )
{
   char* result;
   SCIPallocMemoryArray(scip, &result, SCIP_MAXSTRLEN);
   strcpy(result, DEFAULT_GOBNILP_INPUT_FORMAT);
   if( strcmp(inputformat, "") != 0 )
      strcpy(result, inputformat);
   else
   {
      char suffix[SCIP_MAXSTRLEN];
      int i = 0;
      int j = 0;
      while( filename[i] != 0 )
         i++;
      while( (filename[i] != '.') && (i > 0) )
         i--;
      if( filename[i] == '.' )
      {
         i++;
         while( filename[i] != 0 )
            suffix[j++] = filename[i++];
      }
      suffix[j] = 0;
      if( strcmp(suffix, "") != 0 )
         strcpy(result, suffix);
   }

   if(!(strcmp(result, "jkl") == 0 || strcmp(result, "cip") == 0 ||
         strcmp(result, "dat") == 0 || strcmp(result, "gen") == 0 ||
         strcmp(result, "pss") == 0))
   {
      SCIPwarningMessage(scip, "Input file format not recognised - assuming it is Jaakkola.\n");
      strcpy(result, DEFAULT_GOBNILP_INPUT_FORMAT);
   }

   return result;
}


/** Sets various built-in SCIP parameters to suitable GOBNILP values.
 *
 *  These values can still be overriden using the settings file.
 *
 *  @param scip The SCIP instance the parameters relate to.
 *  @return SCIP_OKAY if the setting worked or an appropriate error otherwise.
 */
SCIP_RETCODE BN_setParamaterDefaults(
   SCIP* scip
   )
{
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/coefdiving/freq"        , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/crossover/freq"         , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq"          , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq"        , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/guideddiving/freq"      , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq"  , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/nlpdiving/freq"         , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/freq"            , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq"   , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/pscostdiving/freq"      , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rens/freq"              , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rootsoldiving/freq"     , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/undercover/freq"        , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/veclendiving/freq"      , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq" , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/randrounding/freq"      , -1) );

   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/vbounds/freq"      , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/locks/freq"      , -1) );

   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/farkasdiving/freq"      , -1) );

   SCIP_CALL( SCIPsetIntParam(scip, "separating/clique/freq"            , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/flowcover/freq"         , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/impliedbounds/freq"     , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/intobj/freq"            , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/mcf/freq"               , -1) );
   SCIP_CALL( SCIPsetCharParam(scip, "separating/efficacynorm"          , 'd') );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxcutsroot"            , 300) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxstallrounds"            , 10) );
   SCIP_CALL( SCIPsetRealParam(scip, "separating/minefficacyroot"        , 0.0005) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/poolfreq"               ,  1) );

   

   SCIP_CALL( SCIPsetCharParam(scip, "nodeselection/childsel"           , 'u') );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts"            ,  0) );
   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype"                  ,  2) );

   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/maxlookahead"            , 1) );
   SCIP_CALL( SCIPsetRealParam(scip, "branching/relpscost/inferenceweight"            , 0.1) );

   SCIP_CALL( SCIPsetIntParam(scip, "display/depth/active"              ,  2) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/vars/active"               ,  0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/poolsize/active"           ,  2) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/lpavgiterations/active"    ,  0) );

   SCIP_CALL( SCIPsetIntParam(scip, "separating/closecuts/freq"         ,  0) );
   SCIP_CALL( SCIPsetBoolParam(scip, "separating/closecuts/separelint"  , FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/oddcycle/freq"          ,  0) );

   SCIP_CALL( SCIPsetRealParam(scip, "separating/gomory/maxbounddist"   ,  1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/gomory/maxrounds"       , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/gomory/maxroundsroot"   , -1) );

   SCIP_CALL( SCIPsetIntParam(scip, "separating/zerohalf/freq"          ,  1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/zerohalf/maxrounds"     , -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/zerohalf/maxroundsroot" , -1) );
   SCIP_CALL( SCIPsetRealParam(scip, "separating/zerohalf/maxbounddist" ,  1) );
   return SCIP_OKAY;
}

/** Suppresses output columns according to the value of parameter gobnilp/verblevelsetcols
 *  @return SCIP_OKAY assuming all is well
 */
SCIP_RETCODE BN_suppresscols(
   SCIP* scip    /*< SCIP instance */
   )
{
   SCIP_Bool verblevelsetscols;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/verblevelsetscols", &verblevelsetscols) );
   /* Simplify the LP output for lower verbosity levels */
   if( verblevelsetscols && SCIPgetVerbLevel(scip) <= 3 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "display/solfound/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/nnodes/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/nodesleft/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/lpiterations/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/lpavgiterations/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/lpcond/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/memused/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/depth/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/maxdepth/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/plungedepth/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/nfrac/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/nexternbranchcands/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/vars/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/conss/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/curconss/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/curcols/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/currows/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/cuts/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/separounds/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/poolsize/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/conflicts/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/strongbranchs/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/pseudoobj/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/lpobj/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/curdualbound/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/estimate/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/avgdualbound/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/cutoffbound/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/primalgap/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/nsols/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/dualbound/active", 0) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/primalbound/active", 0) );

      SCIP_CALL( SCIPsetIntParam(scip, "display/time/active", 2) );
      SCIP_CALL( SCIPsetIntParam(scip, "display/gap/active", 2) );
   }

   return SCIP_OKAY;
}

/** If gobnilp/pasizevars == TRUE adds and sets the branching priorities of variables representing parent set size
 *
 * Creates pasize[i] for each node i which is constrained to be equal to the parent set size of node i
 * @return SCIP_OKAY as long as all is well
 */
static
SCIP_RETCODE addPaSizeVariables(
   SCIP*          scip,         /**< The SCIP instance */
   ParentSetData* psd           /**< Parent set data */
   )
{
   int i;
   int k;
   char s[SCIP_MAXSTRLEN];
   int maxsize;
   SCIP_CONS* cons;

   SCIP_Bool pasizevars;
   int pasizepriority;

   SCIP_PROBDATA* probdata;
   
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/pasizevars", &pasizevars) );

   if( !pasizevars )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pasizepriority", &pasizepriority) );

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->pasize), psd->n) );

   for( i = 0; i < psd->n; ++i )
   {
      maxsize = 0;
      for( k = 0; k < psd->nParentSets[i]; ++k )
         if( psd->nParents[i][k] > maxsize )
            maxsize = psd->nParents[i][k];

      SCIPsnprintf(s, SCIP_MAXSTRLEN, "pasize#%s", psd->nodeNames[i]);
      SCIP_CALL(SCIPcreateVar(scip, &(probdata->pasize[i]), s, 0.0, maxsize,
                              0,
                              SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
      SCIP_CALL( SCIPaddVar(scip, probdata->pasize[i]) );

      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->pasize[i]) );

      SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->pasize[i], pasizepriority) );

      SCIPsnprintf(s, SCIP_MAXSTRLEN, "pasize_cons#%s", psd->nodeNames[i]);
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL,
                                     0, 0,
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->pasize[i], -1) );

      for( k = 0; k < psd->nParentSets[i]; ++k )
         SCIP_CALL( SCIPaddCoefLinear(scip, cons,  psd->PaVars[i][k], psd->nParents[i][k]) );

      SCIP_CALL( SCIPaddCons(scip, cons) );
      /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   return SCIP_OKAY;
}

/** If gobnilp/posindvars == TRUE, adds binary variables indicating the position of nodes in a total order consistent
 *  with a DAG solution
 *
 * Creates n^2 variables posind[i][j], where posind[i][j]==1 iff variable i is in position j in some total order consistent with a DAG solution
 * This is done by the following constraints: pos(j) - pos(i) + n*I(j is a parent of i) =< n-1
 * pos(i) and pos(j) are integers denoted the index of i and j in a total order.
 * They are defined as a linear combination of posind variables in the obvious way.
 * @return SCIP_OKAY as long as all is well
 */
static
SCIP_RETCODE addPosindVariables(
   SCIP* scip,          /**< SCIP data structure */
   ParentSetData* psd   /**< parent sets data structure */
   )
{
   int i;
   int j;
   int k;
   /* int pos; */
   /* int p; */

   SCIP_VAR** vars;
   SCIP_CONS* cons;

   SCIP_VAR** indicvars;
   SCIP_Real* indicvals;

   char s[SCIP_MAXSTRLEN];

   SCIP_Bool posindvars;
   int posindvarpriority;

   SCIP_VAR* arrow_i_j;

   int nvals;
   int nvars;

   SCIP_PROBDATA* probdata;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/posindvars", &posindvars) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/posindvarpriority", &posindvarpriority) );

   if( !posindvars )
      return SCIP_OKAY;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->posind), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, psd->n) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &indicvars, 2*psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &indicvals, 2*psd->n) );

   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->posind[i]), psd->n) );
      for( j = 0; j < psd->n; ++j )
      {
         /* posind[i][j]==1 if variable i is in position j on 'the' total order */
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "Ipos#%d#%d", i, j);
         SCIP_CALL(SCIPcreateVar(scip, &(probdata->posind[i][j]), s, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
                                 TRUE,  /* if TRUE then used in the root */
                                 TRUE,  /* if TRUE then removable */
                                 NULL, NULL, NULL, NULL, NULL));
         SCIP_CALL( SCIPaddVar(scip, probdata->posind[i][j]) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip,probdata->posind[i][j],posindvarpriority - j)  );
      }
   }

   for( i = 0; i < psd->n; ++i )
   {
      /* each variable in exactly one position */
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "oneposd#%d", i);
      SCIP_CALL(SCIPcreateConsSetpart(scip, &cons, s, psd->n, probdata->posind[i],
                                      TRUE, /*initial*/
                                      TRUE, /*separate*/
                                      TRUE, TRUE,/*enforce, check */
                                      TRUE,
                                      FALSE, FALSE, FALSE, FALSE, FALSE));
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* i is here indexing position, only one variable for each position */
      for( j = 0; j < psd->n; ++j )
         vars[j] = probdata->posind[j][i];
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "onevar#%d", i);
      SCIP_CALL(SCIPcreateConsSetpart(scip, &cons, s, psd->n, vars,
                                      TRUE,
                                      TRUE,
                                      TRUE, TRUE,
                                      TRUE,
                                      FALSE, FALSE, FALSE, FALSE, FALSE));
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
            

   /* if j is a parent of i then pos(j) - pos(i) =< -1 */

   for( i = 0; i < psd->n; ++i )
   {
      for( j = 0; j < psd->n; ++j )
      {

         if( i == j )
            continue;

         arrow_i_j = get_arrow(psd,i,j);
         if( arrow_i_j == NULL )
            continue;

         nvars = 0;
         nvals = 0;

         /* k indexes position */
         for( k = 0; k < psd->n; ++k )
         {
            indicvars[nvars++] = probdata->posind[i][k];
            indicvals[nvals++] = -k;
            indicvars[nvars++] = probdata->posind[j][k];
            indicvals[nvals++] = k;
         }

         assert( nvals == 2*psd->n );
         assert( nvals == nvars );

         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "gen#%d#%d", i, j);
         SCIP_CALL(SCIPcreateConsIndicator(
               scip, &cons, s, arrow_i_j, nvars, indicvars, indicvals, -1,
               TRUE, TRUE,
               TRUE, TRUE,
               TRUE,
               FALSE,FALSE,FALSE,FALSE));

         SCIP_CALL( SCIPaddCons(scip, cons) );
         /* SCIP_CALL( SCIPprintCons(scip, cons, NULL)  );*/
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      }
   }


   SCIPfreeMemoryArray(scip, &vars);

   SCIPfreeMemoryArray(scip, &indicvars);
   SCIPfreeMemoryArray(scip, &indicvals);

   return SCIP_OKAY;
}

/** If gobnilp/posvars == TRUE, creates integer variables giving a generation number for each node
 * If a node is a founder its generation number is 0
 * If j is a parent of i its generation number is strictly less than i's
 * An upper bound on any such variable can be set using gobnilp/maxpos
 * Uses an Indicator constraint
 * @return SCIP_OKAY as long as all is well
 */
static
SCIP_RETCODE addPosVariables(
   SCIP* scip,           /**< SCIP data structure */
   ParentSetData* psd    /**< parent sets data structure */
   )
{
   int i;
   int j;
   int k;

   SCIP_CONS* cons;

   char s[SCIP_MAXSTRLEN];

   SCIP_Bool posvars;

   int maxpos;

   SCIP_Real indicvals[2] = {1, -1};
   SCIP_VAR* indicvars[2];

   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real* vals;

   SCIP_Bool total;
   SCIP_VAR* arrow_i_j;

   int posvarpriority;

   SCIP_PROBDATA* probdata;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/posvars", &posvars) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/posvarpriority", &posvarpriority) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/maxpos", &maxpos) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/postotal", &total) );


   if( !posvars )
      return SCIP_OKAY;

   if( total || maxpos == -1 )
      maxpos = (psd->n) - 1;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->posvars, psd->n) );

   for( i = 0; i < psd->n; ++i )
   {
      /* pos[i] is position of i in a partial (or total) order */
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "pos#%d", i);
      SCIP_CALL(SCIPcreateVar(scip, &(probdata->posvars[i]), s, 0.0, maxpos, 0.0, SCIP_VARTYPE_INTEGER,
                              TRUE,  /*if TRUE then used in the root */
                              TRUE,   /*if TRUE then removable */
                              NULL, NULL, NULL, NULL, NULL));
      SCIP_CALL( SCIPaddVar(scip, probdata->posvars[i]) );
      SCIP_CALL( SCIPchgVarBranchPriority(scip,probdata->posvars[i],posvarpriority)  );
   }


   /* founders at position 0 if order not total */

   if( !total )
      for( i = 0; i < psd->n; ++i )
      {
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "founderpos#%d", i);
         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL, -SCIPinfinity(scip), maxpos,
               TRUE,
               TRUE,
               TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->posvars[i], 1) );
         

         for( k = 0; k < psd->nParentSets[i]; ++k )
            if( psd->nParents[i][k] == 0 )
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][k], maxpos) );
               break;
            }
         SCIP_CALL( SCIPaddCons(scip, cons) );
         /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }

   /* if a non-empty parent set then pos[i] > 0 */
   for( i = 0; i < psd->n; ++i )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "nonfounderpos#%d", i);
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL, 1, SCIPinfinity(scip),
            TRUE,
            TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->posvars[i], 1) );
      for( k = 0; k < psd->nParentSets[i]; ++k )
         if( psd->nParents[i][k] == 0 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][k], 1) );
            break;
         }
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   
   /* we have a bound on the sum of all pos variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, psd->n) );
   for( i = 0; i < psd->n; ++i )
      vals[i] = 1;
   if( total )
   {
      lb = maxpos*(maxpos+1)/2;
      ub = maxpos*(maxpos+1)/2;
   }
   else
   {
      lb = -SCIPinfinity(scip);
      ub = maxpos*(maxpos-1)/2 + ((psd->n)-maxpos)*maxpos;
   }
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "posvarssumcons", psd->n, probdata->posvars, vals, lb, ub) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPfreeBufferArray(scip, &vals);

   /* if j is a parent of i then pos(j) - pos(i) <= -1 */

   for( i = 0; i < psd->n; ++i )
   {
      for( j = 0; j < psd->n; ++j )
      {
         if( i == j )
            continue;
         
         /* add indicator constraint */

         indicvars[0] = probdata->posvars[j];
         indicvars[1] = probdata->posvars[i];
         arrow_i_j = get_arrow(psd,i,j);
         if( arrow_i_j != NULL )
         {
            SCIP_CALL(SCIPcreateConsBasicIndicator(
                  scip, &cons, "indic", arrow_i_j, 2 , indicvars, indicvals, -1));
            SCIP_CALL( SCIPaddCons(scip, cons) );
            /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
   }

   return SCIP_OKAY;
}


/** If gobnilp/ancestorvars == TRUE, creates binary variables encoding the ancestor relation
 * Creates binary variables anc[i][j] where anc[i][j]==1 iff j is an ancestor of i
 * If gobnilp/useconspartialordering == TRUE, the partial order constraint is used on these variables
 * @return SCIP_OKAY as long as all is well
 */
static
SCIP_RETCODE addAncestorVariables(
   SCIP* scip,           /**< SCIP data structure */
   ParentSetData* psd    /**< parent sets data structure */
   )
{

   int i;
   int j;

   SCIP_Bool ancestorvars;
   int ancestorvarpriority;
   SCIP_Bool useconspartialordering;
   SCIP_Bool looseancestors;
   char s[SCIP_MAXSTRLEN];
   
   SCIP_CONS* cons;

   SCIP_RETCODE retcode = SCIP_OKAY;

   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   assert(psd != NULL);
   
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/ancestorvars", &ancestorvars) );

   if( !ancestorvars )
      return SCIP_OKAY;


   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/ancestorvarpriority",&ancestorvarpriority) );
   
   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   
   /* create variables */

   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->ancestorvars, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->ancestorvars[i]), psd->n) );
      for( j = 0; j < psd->n; ++j )
      {
         if( i == j )
            continue;

         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "anc#%d#%d", i, j);
         SCIP_CALL(SCIPcreateVar(scip, &(probdata->ancestorvars[i][j]), s, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
                                 FALSE,  /*if TRUE then used in the root */
                                 TRUE,   /*if TRUE then removable */
                                 NULL, NULL, NULL, NULL, NULL));
         /* SCIP_CALL(SCIPcreateVarBasic(scip, &(probdata->ancestorvars[i][j]), s, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) ); */
         SCIP_CALL( SCIPaddVar(scip, probdata->ancestorvars[i][j]) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->ancestorvars[i][j], ancestorvarpriority) );
      }
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/useconspartialordering", &useconspartialordering) );

   if( useconspartialordering )
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/looseancestors", &looseancestors) );

      /* SCIP_CALL(SCIPcreateConsBasicPartialOrdering(scip, &cons, "ANC", psd->n, probdata->ancestorvars, psd, !looseancestors) ); */
      SCIP_CALL(SCIPcreateConsPartialOrdering(scip, &cons, "ANC", psd->n, probdata->ancestorvars, psd, !looseancestors,
            FALSE, /*not initial*/
                                             FALSE,
                                             TRUE,
                                             TRUE,
                                             TRUE,
                                             FALSE,
                                             FALSE,
                                             FALSE,
                                             FALSE,
                                             FALSE));
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   else
   {
      SCIPerrorMessage("Currently must have gobnilp/useconspartialordering == TRUE if gobnilp/ancestorvars == TRUE");
      retcode = SCIP_ERROR;
   }
   

   return retcode;

}
   

/** If gobnilp/totalordervars == TRUE, creates, and sets the branching priority for, binary variables encoding a total order consistent with a DAG solution
 * Creates binary variables order[i][j] where order[i][j]==1 if i comes after j in the total ordering.
 * If gobnilp/useconslinearordering == TRUE, the linear ordering constraint handler that comes as an example with SCIP is used
 * If gobnilp/inittransitivity ==  TRUE then all n^3 transitivity constraints are added to the problem
 * If gobnilp/bestparentsfororder == TRUE then, for each node, a constraint is added so that
 *     the highest scoring parent set consistent with the total order must be selected
 * @return SCIP_OKAY as long as all is well
 */
static
SCIP_RETCODE addTotalorderVariables(
   SCIP* scip,           /**< SCIP data structure */
   ParentSetData* psd    /**< parent sets data structure */
   )
{
   int i;
   int j;
   int k;
   int kk;
   int l;

   SCIP_VAR* vars[2];
   SCIP_VAR** bestpavars;

   char s[SCIP_MAXSTRLEN];
   int nvars;

   int biggest;

   SCIP_CONS* cons;

   SCIP_Bool totalordervars;
   SCIP_Bool useconslinearordering;
   SCIP_Bool inittransitivity;
   SCIP_Bool bestparentsfororder;
   int totalordervarpriority;

   SCIP_VAR* initvars[3] = {NULL, NULL, NULL};
   SCIP_Longint weights[3] = {1, 1, 1};
   
   SCIP_VAR* arrow_i_j;

   SCIP_PROBDATA* probdata;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/totalordervars", &totalordervars) );

   if( !totalordervars )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/totalordervarpriority", &totalordervarpriority) );

   /* create variables */

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->totalordervars), psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->totalordervars[i]), psd->n) );
      for( j = 0; j < psd->n; ++j )
      {
         if( i == j )
            continue;

         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "ord#%d#%d", i, j);
         SCIP_CALL(SCIPcreateVar(scip, &(probdata->totalordervars[i][j]), s, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY,
                                 TRUE,  /*if TRUE then used in the root */
                                 TRUE,   /*if TRUE then removable */
                                 NULL, NULL, NULL, NULL, NULL));
         SCIP_CALL( SCIPaddVar(scip, probdata->totalordervars[i][j]) );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->totalordervars[i][j], totalordervarpriority) );
      }
   }

   /* parents come before children
      link to arrow variables */

   for( i = 0; i < psd->n; ++i )
   {
      for( j = 0; j < psd->n; ++j )
      {
         if( i == j )
            continue;

         arrow_i_j = get_arrow(psd,i,j);
         if( arrow_i_j != NULL )
         {
            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "tocons#%d#%d", i, j);
            vars[0] = arrow_i_j;
            SCIP_CALL( SCIPgetNegatedVar(scip,probdata->totalordervars[i][j],&(vars[1])) );
            SCIP_CALL( SCIPcreateConsSetpack(
                  scip,
                  &cons,
                  s,
                  2,
                  vars,   /*vars,*/
                  TRUE,   /*initial,*/
                  TRUE,   /* separate, */
                  TRUE,   /* enforce */
                  TRUE,   /* check */
                  TRUE,   /* propagate */
                  FALSE,  /* local */
                  FALSE,  /* modifiable */
                  FALSE,     /* dynamic */
                  FALSE,  /* removable */
                  FALSE   /* stickingatnode */
                  ));
            SCIP_CALL( SCIPaddCons(scip, cons) );
            /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
   }

   /* always add anti-symmetry constraints */

   for( i = 0; i < psd->n; ++i )
   {
      for( j = 0; j < psd->n; ++j )
      {
         if( i == j )
            continue;

         initvars[0] = probdata->totalordervars[i][j];
         initvars[1] = probdata->totalordervars[j][i];

         /* add anti-symmetry constraint */

         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "asym#%d#%d", i, j);
         SCIP_CALL(SCIPcreateConsSetpart(
                      scip,
                      &cons,
                      s,
                      2,
                      initvars,     /*vars,*/
                      TRUE,   /*initial,*/
                      TRUE,   /* separate, */
                      TRUE,   /* enforce */
                      TRUE,   /* check */
                      TRUE,   /* propagate */
                      FALSE,  /* local */
                      FALSE,  /* modifiable */
                      FALSE,     /* dynamic */
                      FALSE,  /* removable */
                      FALSE   /* stickingatnode */
                   ));
         SCIP_CALL( SCIPaddCons(scip, cons) );
         /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }


   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/inittransitivity", &inittransitivity) );

   if( inittransitivity )
   {
      for( i = 0; i < psd->n; ++i )
      {
         for( j = 0; j < psd->n; ++j )
         {
            if( i == j )
               continue;

            initvars[0] = probdata->totalordervars[i][j];

            for( k = 0; k < psd->n; ++k )
            {
               if( i == k || j == k )
                  continue;

               initvars[1] = probdata->totalordervars[j][k];
               initvars[2] = probdata->totalordervars[k][i];

               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "trans#%d#%d#%d", i, j, k);
               SCIP_CALL(SCIPcreateConsKnapsack(
                            scip,
                            &cons,
                            s,
                            3,
                            initvars,    /*vars,*/
                            weights,
                            2,
                            TRUE,     /*initial,*/
                            TRUE,     /* separate, */
                            TRUE,     /* enforce */
                            TRUE,     /* check */
                            TRUE,     /* propagate */
                            FALSE,    /* local */
                            FALSE,    /* modifiable */
                            FALSE,    /* dynamic */
                            FALSE,    /* removable */
                            FALSE     /* stickingatnode */
                         ));
               SCIP_CALL( SCIPaddCons(scip, cons) );
               /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }
      }
   }


   SCIPgetBoolParam(scip, "gobnilp/bestparentsfororder", &bestparentsfororder);

   if( bestparentsfororder )
   {
      for( i = 0; i < psd->n; ++i )
      {

         /* find size of biggest parent set */
         biggest = 0;
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            if( psd->nParents[i][k] > biggest )
               biggest = psd->nParents[i][k];
         }

         /* allocate enough space for the variables in the constraint */
         SCIP_CALL( SCIPallocBufferArray(scip, &bestpavars, biggest - 1 + psd->nParentSets[i]) );

         for( k = 0; k < psd->nParentSets[i] - 1; ++k )
         {
            nvars = 0;
            /* either one of the parents comes later in the ordering ... */
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               SCIP_CALL(SCIPgetNegatedVar(scip,
                                           probdata->totalordervars[i][psd->ParentSets[i][k][l]],
                                           &bestpavars[nvars++]
                                          ));
            }
            /* ... or some 'better' parent set chosen ... */
            for( kk = 0; kk < k; ++kk )
            {
               bestpavars[nvars++] = psd->PaVars[i][kk];
            }
            /* ... or this one is chosen */
            assert(kk == k);
            bestpavars[nvars++] = psd->PaVars[i][kk];

            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "select%s", SCIPvarGetName(psd->PaVars[i][k]));
            SCIP_CALL(SCIPcreateConsSetcover(
                         scip,
                         &cons,
                         s,
                         nvars,    /*nvars,*/
                         bestpavars,     /*vars,*/
                         TRUE,     /*initial,*/
                         TRUE,     /* separate, */
                         TRUE,     /* enforce */
                         TRUE,     /* check */
                         TRUE,     /* propagate */
                         FALSE,    /* local */
                         FALSE,    /* modifiable */
                         FALSE,    /* dynamic */
                         FALSE,    /* removable */
                         FALSE     /* stickingatnode */
                      ));
            SCIP_CALL( SCIPaddCons(scip, cons) );
            /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
         SCIPfreeBufferArray(scip, &bestpavars);
      }
   }




   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/useconslinearordering", &useconslinearordering) );

   if( useconslinearordering )
   {
      /* constrain order variables with linear ordering constraint */

      SCIP_CALL(SCIPcreateConsLOP(scip, &cons, "LOP", psd->n, probdata->totalordervars,
                                             TRUE,
                                             TRUE,
                                             TRUE,
                                             TRUE,
                                             TRUE,
                                             FALSE,
                                             FALSE,
                                             FALSE,
                                             FALSE,
                                             FALSE));
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;

}

/* Add any variables for k best MEC learning */
static
SCIP_RETCODE addKBestMECVariables(
   SCIP* scip,           /**< SCIP data structure */
   ParentSetData* psd    /**< parent sets data structure */
   )
{
   int i;
   int j;
   int k;
   int l;
   int m;

   SCIP_Bool kbestMEC;

   char s[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   SCIP_VAR* var;

   SCIP_Bool needed;
   SCIP_Bool ifound;
   SCIP_Bool jfound;

   SCIP_VAR* vars[2];
   SCIP_VAR* edge_i_j;

   SCIP_PROBDATA* probdata;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/kbestMEC", &kbestMEC) );

   if( !kbestMEC )
      return SCIP_OKAY;
   
   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   

   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->im_vars, psd->n) );
   for( i = 0; i < psd->n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->im_vars[i]), psd->n) );
      for( j = 0; j < psd->n; j++ )
            SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->im_vars[i][j]), psd->n) );
   }

   /* Add the immorality variables */
   for( i = 0; i < psd->n - 1; i++ )
      for( j = i + 1; j < psd->n; j++ )
         for( k = 0; k < psd->n; k++ )
            if( i != k && j != k )
            {
               /* Create the variable representing the relationship i->k<-j */
               SCIPsnprintf(s, SCIP_MAXSTRLEN, "PA#%s#%s#%s", psd->nodeNames[i], psd->nodeNames[j], psd->nodeNames[k]);
               SCIP_CALL( SCIPcreateVarBasic(scip, &var, s, 0.0, 1.0, 0, SCIP_VARTYPE_BINARY) );
               SCIP_CALL( SCIPaddVar(scip, var) );
               SCIP_CALL( SCIPchgVarBranchPriority(scip, var, -1) );
               
               /* Create the constraint linking it to the other vars */
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, s, 0, NULL, NULL, 0, 0) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, var, -1) );
               needed = FALSE;
               for( l = 0; l < psd->nParentSets[k]; l++ )
               {
                  ifound = FALSE;
                  jfound = FALSE;
                  for( m = 0; m < psd->nParents[k][l]; m++ )
                  {
                     if( psd->ParentSets[k][l][m] == i )
                        ifound = TRUE;
                     if( psd->ParentSets[k][l][m] == j )
                        jfound = TRUE;
                  }
                  if( ifound && jfound )
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[k][l], 1) );
                     needed = TRUE;
                  }
               }
               if( needed )
               {
                  SCIP_CALL( SCIPaddCons(scip, cons) );
               }
               else
               {
                  SCIP_Bool feasible;
                  SCIP_Bool fixed;
                  SCIP_CALL( SCIPfixVar(scip, var, 0, &feasible, &fixed) );
               }
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               
               /* Create the variable representing the immorality i->k<-j */
               SCIPsnprintf(s, SCIP_MAXSTRLEN, "IM#%s#%s#%s", psd->nodeNames[i], psd->nodeNames[j], psd->nodeNames[k]);
               SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->im_vars[i][j][k]), s, 0.0, 1.0, 0, SCIP_VARTYPE_BINARY) );
               SCIP_CALL( SCIPaddVar(scip, probdata->im_vars[i][j][k]) );
               SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->im_vars[i][j][k], -1) );
               
               /* Create the constraint linking it to the other variables */
               if( needed )
               {
                  vars[0] = var;
                  edge_i_j = get_edge(psd,i,j);
                  if( edge_i_j != NULL )
                  {
                     SCIP_CALL( SCIPgetNegatedVar(scip, edge_i_j, &(vars[1])) );
                     SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, s, probdata->im_vars[i][j][k], 2, vars) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPcreateConsBasicAnd(scip, &cons, s, probdata->im_vars[i][j][k], 1, vars) );
                  }
                  SCIP_CALL( SCIPaddCons(scip, cons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               }
               else
               {
                  SCIP_Bool feasible;
                  SCIP_Bool fixed;
                  SCIP_CALL( SCIPfixVar(scip, probdata->im_vars[i][j][k], 0, &feasible, &fixed) );
               }
            }
   
   return SCIP_OKAY;
}


/** Adds variables representing directed and undirected edges
 *  The branching priority for these variables is set at 10, so typically SCIP will branch on these variables rather than parent set variables
 * @return SCIP_OKAY if all is well
 */
static
SCIP_RETCODE addArrowVariables(
   SCIP* scip,          /**< SCIP data structure */
   ParentSetData* psd   /**< parent sets data structure */
   )
{
   char s[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   int i;
   int j;
   int k;
   int l;
   SCIP_Bool** ParentSetsWithArrowFlag;
   SCIP_VAR** include;
   SCIP_VAR** exclude;
   int n_include;
   int n_exclude;
   SCIP_VAR* neg_var;
   SCIP_VAR* edgeconsvars[3];

   SCIP_Real edgepenalty;

   SCIP_VAR* arrow_i_j;
   SCIP_VAR* arrow_j_i;
   SCIP_VAR* edge_i_j;

   int arrowfamilylink;

   int num_arrows = 0;
   int num_edges = 0;

   SCIP_Bool arrowvarsinitial;
   SCIP_Bool arrowvarsremovable;
   SCIP_Bool edgevarsinitial;
   SCIP_Bool edgevarsremovable;

   SCIP_Bool arrowconsinitial;
   SCIP_Bool edgeconsinitial;

   int nvars;
   
   SCIP_CALL( SCIPgetRealParam(scip, "gobnilp/edge_penalty", &edgepenalty) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/arrowfamilylink", &arrowfamilylink) );

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/arrowvarsinitial", &arrowvarsinitial) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/arrowvarsremovable", &arrowvarsremovable) );

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/edgevarsinitial", &edgevarsinitial) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/edgevarsremovable", &edgevarsremovable) );

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/arrowconsinitial", &arrowconsinitial) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/edgeconsinitial", &edgeconsinitial) );
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &ParentSetsWithArrowFlag, psd->n) );

   /* set up hash table for arrow and edge variables */
   SCIP_CALL( hashtableCreateArrow(scip, psd) );

   for( i = 0; i < psd->n; ++i )
   {
      for( j = 0; j < psd->n; ++j )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSetsWithArrowFlag[j]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &include, psd->nParentSets[i] + 1) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &exclude, psd->nParentSets[i] + 1) );

      /* record parents in each parent set */
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         for( j = 0; j < psd->n; ++j )
            if( j != i )
               ParentSetsWithArrowFlag[j][k] = FALSE;
         for( l = 0; l < psd->nParents[i][k]; ++l )
            ParentSetsWithArrowFlag[psd->ParentSets[i][k][l]][k] = TRUE;
      }

      /* for each parent j (of i) compute in which parent sets it exists, and which not */
      for( j = 0; j < psd->n; ++j )
      {
         
         if( j == i )
            continue;

         n_include = 0;
         n_exclude = 0;
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            if( ParentSetsWithArrowFlag[j][k] )
               include[n_include++] = psd->PaVars[i][k];
            else
               exclude[n_exclude++] = psd->PaVars[i][k];
         }

         /* never create arrow variables for arrows which can never occur */
         if( n_include == 0 )
            continue;

         /* /\* if arrow only appears in one family, just set arrow variable to be the family variable *\/ */
         /* if( arrowfamilylink != 1 && n_include == 1 ) */
         /* { */
         /*    SCIP_CALL( put_arrow(scip, psd, i, j, include[0]) ); */
         /*    continue; */
         /* } */

         /* /\* if arrow missing from only one family, just set arrow variable to be the negation of that family variable *\/ */
         /* if( arrowfamilylink != 1 && n_exclude == 1 ) */
         /* { */
         /*    SCIP_CALL( SCIPgetNegatedVar(scip, exclude[0], &neg_var) ); */
         /*    SCIP_CALL( put_arrow(scip, psd, i, j, neg_var) ); */
         /*    continue; */
         /* } */

         /* create an arrow variable only ( we know j can be a parent of i )*/
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "arrow#%s#%s", psd->nodeNames[i],psd->nodeNames[j]);
         SCIP_CALL( SCIPcreateVar(scip, &arrow_i_j, s, 0.0, 1.0, 0,
               SCIP_VARTYPE_BINARY,
               arrowvarsinitial,  /* initial? */
               arrowvarsremovable,   /* removable? */
               NULL, NULL, NULL, NULL, NULL));
         SCIP_CALL( SCIPaddVar(scip, arrow_i_j) );
         num_arrows++;
         SCIP_CALL( put_arrow(scip, psd, i, j, arrow_i_j) );

               
         SCIP_CALL( SCIPgetNegatedVar(scip, arrow_i_j, &neg_var) );

         include[n_include++] = neg_var;
         exclude[n_exclude++] = arrow_i_j;

         if( arrowfamilylink == 1 )
         {
            SCIPsnprintf(s, SCIP_MAXSTRLEN, "setpart1#%s#%s", psd->nodeNames[i],psd->nodeNames[j]);
            SCIPcreateConsSetpack(scip, &cons, s, n_include, include,
               arrowconsinitial, /* SCIP_Bool  initial,*/
               TRUE, /* SCIP_Bool separate, */
               TRUE, /* SCIP_Bool  enforce, */
               TRUE, /* SCIP_Bool  check, */
               TRUE,  /* SCIP_Bool  propagate, */
               FALSE, /* SCIP_Bool  local, */
               FALSE, /* SCIP_Bool  modifiable, */
               FALSE, /* SCIP_Bool  dynamic, */
               FALSE, /* SCIP_Bool removable, might change this */
               FALSE  /* SCIP_Bool  stickingatnode */
               );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            /* SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); */
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }

         
         if( arrowfamilylink == 2 || arrowfamilylink == 4 )
         {
            SCIPsnprintf(s, SCIP_MAXSTRLEN, "setpart1#%s#%s", psd->nodeNames[i],psd->nodeNames[j]);
            SCIPcreateConsSetpart(scip, &cons, s, n_include, include,
               arrowconsinitial, /* SCIP_Bool  initial,*/
               TRUE, /* SCIP_Bool separate, */
               TRUE, /* SCIP_Bool  enforce, */
               TRUE, /* SCIP_Bool  check, */
               TRUE,  /* SCIP_Bool  propagate, */
               FALSE, /* SCIP_Bool  local, */
               FALSE, /* SCIP_Bool  modifiable, */
               FALSE, /* SCIP_Bool  dynamic, */
               FALSE, /* SCIP_Bool removable, might change this */
               FALSE  /* SCIP_Bool  stickingatnode */
               );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            /* SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); */
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }

         if( arrowfamilylink == 3 || arrowfamilylink == 4 )
         {
            SCIPsnprintf(s, SCIP_MAXSTRLEN, "setpart2#%s#%s", psd->nodeNames[i],psd->nodeNames[j]);
            SCIPcreateConsSetpart(scip, &cons, s, n_exclude, exclude,
               arrowconsinitial, /* SCIP_Bool   initial,*/
               TRUE, /* SCIP_Bool separate, */
               TRUE, /* SCIP_Bool   enforce, */
               TRUE, /* SCIP_Bool   check, */
               TRUE,  /* SCIP_Bool  propagate, */
               FALSE, /* SCIP_Bool  local, */
               FALSE, /* SCIP_Bool  modifiable, */
               FALSE, /* SCIP_Bool  dynamic, */
               FALSE,  /* SCIP_Bool removable, might change this */
               FALSE  /* SCIP_Bool  stickingatnode */
               );
            
            SCIP_CALL( SCIPaddCons(scip, cons) );
            /* SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); */
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }

      for( j = 0; j < psd->n; ++j )
         SCIPfreeMemoryArray(scip, &(ParentSetsWithArrowFlag[j]));

      SCIPfreeMemoryArray(scip, &include);
      SCIPfreeMemoryArray(scip, &exclude);

   }

   for( i = 0; i < psd->n; ++i )
   {
      for( j = i+1; j < psd->n; ++j )
      {
         arrow_i_j = get_arrow(psd,i,j);
         arrow_j_i = get_arrow(psd,j,i);

         if( arrow_i_j == NULL && arrow_j_i == NULL )
            continue;

         SCIPsnprintf(s, SCIP_MAXSTRLEN, "edge#%s#%s", psd->nodeNames[i],psd->nodeNames[j]);
         SCIP_CALL( SCIPcreateVar(scip, &edge_i_j, s, 0.0, 1.0, -edgepenalty,
               SCIP_VARTYPE_BINARY,
               edgevarsinitial,     /* initial?    */
               edgevarsremovable,   /* removable ? */
               NULL, NULL, NULL, NULL, NULL));
         SCIP_CALL( SCIPaddVar(scip, edge_i_j) );
         num_edges++;
         SCIP_CALL( put_edge(scip, psd, i, j, edge_i_j) );
         
         
         if( arrow_i_j != NULL)
         {
            edgeconsvars[1] = arrow_i_j;
            
            if( arrow_j_i != NULL)
            {
               edgeconsvars[2] = arrow_j_i;
               nvars = 3;
            }
            else
            {
               nvars = 2;
            }
         }
         else
         {
            assert(arrow_j_i != NULL);
            edgeconsvars[1] = arrow_j_i;
            nvars = 2;
         }
      
         SCIP_CALL( SCIPgetNegatedVar(scip, edge_i_j, &neg_var) );
         edgeconsvars[0] = neg_var;
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "setpart3#%s#%s", psd->nodeNames[i],psd->nodeNames[j]);
         SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, s, nvars, edgeconsvars,
               edgeconsinitial,
               TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE, FALSE, FALSE, FALSE) );

         
         SCIP_CALL( SCIPaddCons(scip, cons) );
         /* SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); */
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   SCIPfreeMemoryArray(scip, &ParentSetsWithArrowFlag);

   /* for( i = 0; i < psd->n; ++i ) */
   /*    for( j = 0; j < psd->n; ++j ) */
   /*    { */
   /*       printf("%d %d ",i,j); */
   /*       arrow_i_j = get_arrow(psd,i,j); */
   /*       if( arrow_i_j == NULL ) */
   /*          printf("NULL\n"); */
   /*       else */
   /*          printf("%s\n",SCIPvarGetName(arrow_i_j)); */
   /*    } */

   SCIPdebugMessage("%d arrow and %d edge variables added\n", num_arrows, num_edges);

   return SCIP_OKAY;
}

/** The imsetvar for a subset of nodes = 1 iff one of the elements in the subset has all the others as parents. It is binary due to acyclicity.
 *  If gobnilp/imsetvars == TRUE, creates and sets the branching priority for, variables representing imsets for subsets of nodes of size 3
 *  Note that the imsetvar for a subset of size 2 is an 'edge' variable - these are always created.
 *  @return SCIP_OKAY, if all is well
*/
static
SCIP_RETCODE addImsetVariables(
   SCIP* scip,           /**< SCIP data structure */
   ParentSetData* psd    /**< parent sets data structure */
   )
{
   int k;
   int l;
   char s[SCIP_MAXSTRLEN];
   int i0;
   int i1;
   int i2;
   SCIP_CONS* cons;
   int found;

   SCIP_Bool imsetvars;
   int imsetvarpriority;

   int child;
   int current_parent;
   int parent1;
   int parent2;
   int tmp;
   int iter;
   int nvars;

   SCIP_Bool decomposable;
   SCIP_VAR* edge_var;
   int j0;
   int j1;

   SCIP_PROBDATA* probdata;
   
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/decomposable", &decomposable) );

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/imsetvars", &imsetvars) );

   if( !imsetvars )
      return SCIP_OKAY;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/imsetvarpriority", &imsetvarpriority) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->imsetvars), psd->n) );
   for( i0 = 0; i0 < psd->n; ++i0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->imsetvars[i0]), psd->n - (i0+1)) );
      for( i1 = i0 + 1; i1 < psd->n; ++i1 )
      {
         int ii1 = i1 - i0 - 1;
         SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->imsetvars[i0][ii1]), psd->n - (i1+1)) );
         for( i2 = i1 + 1; i2 < psd->n; ++i2 )
         {
            int ii2 = i2 - i1 - 1;

            SCIPsnprintf(s, SCIP_MAXSTRLEN, "imcons#%s#%s#%s", psd->nodeNames[i0], psd->nodeNames[i1], psd->nodeNames[i2]);
            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL, 0, 0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
            child = i0;
            parent1 = i1;
            parent2 = i2;
            nvars = 0;

            for( iter = 0; iter < 3; ++iter )
            {
               for( k = 0; k < psd->nParentSets[child]; ++k )
               {
                  found = 0;
                  for( l = 0; l < psd->nParents[child][k]; ++l )
                  {
                     current_parent = psd->ParentSets[child][k][l];
                     if( current_parent == parent1 || current_parent == parent2 )
                        found++;
                     if( found == 2 )
                     {
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[child][k], 1) );
                        nvars++;
                        break;
                     }
                  }
               }
               
               tmp = child;
               child = parent1;
               parent1 = parent2;
               parent2 = tmp;
            }

            if( nvars > 0 )
            {
               
               SCIPsnprintf(s, SCIP_MAXSTRLEN, "im#%s#%s#%s", psd->nodeNames[i0], psd->nodeNames[i1], psd->nodeNames[i2]);
               SCIP_CALL(SCIPcreateVar(scip, &(probdata->imsetvars[i0][ii1][ii2]), s, 0.0, 1.0, 0,
                     SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
               SCIP_CALL( SCIPaddVar(scip, probdata->imsetvars[i0][ii1][ii2]) );
               SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip,probdata->imsetvars[i0][ii1][ii2]) );
               SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->imsetvars[i0][ii1][ii2], imsetvarpriority) );
               
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->imsetvars[i0][ii1][ii2], -1) );
               
               SCIP_CALL( SCIPaddCons(scip, cons) );
               /*SCIP_CALL( SCIPprintCons(scip, cons, NULL)  ); */
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               
               if( decomposable )
               {
                  j0 = i0;
                  j1 = i1;
                  
                  for( iter = 0; iter < 3; ++iter )
                  {
                     edge_var = get_edge(psd,j0,j1);
                     if( edge_var != NULL )
                     {
                        SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "todo", 0, NULL, NULL, 0, SCIPinfinity(scip),
                              TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->imsetvars[i0][ii1][ii2], -1) );
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, edge_var, 1) );
                        SCIP_CALL( SCIPaddCons(scip, cons) );
                        /*SCIP_CALL( SCIPprintCons(scip, cons, NULL)  ); */
                        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                     }
                     if( iter == 0 )
                        j1 = i2;
                     else
                        j0 = i1;
                  }
               }
            }
            else
            {
               probdata->imsetvars[i0][ii1][ii2] = NULL;
            }
         }
      }
   }
   return SCIP_OKAY;
}


/** Adds variables which count how many children a node has */
static
SCIP_RETCODE addNChildrenVars(
   SCIP* scip,
   const ParentSetData* psd
   )
{
   char s[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   int i;        /* indexes variables */
   int j;        /* indexes variables */
   int maxnchildren;
   int nchildrenvarpriority;
   SCIP_Bool nchildrenvars;
   SCIP_VAR* arrow_j_i;

   SCIP_PROBDATA* probdata = NULL;
   
   assert(scip != NULL);
   assert(psd != NULL);

   SCIPgetIntParam(scip, "gobnilp/maxnchildren", &maxnchildren);
   SCIPgetIntParam(scip, "gobnilp/nchildrenvarpriority", &nchildrenvarpriority);
   SCIPgetBoolParam(scip, "gobnilp/nchildrenvars", &nchildrenvars);

   /* just return if no variables to create, and no constraint to post */
   if( !nchildrenvars && maxnchildren == -1 )
      return SCIP_OKAY;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   
   if( maxnchildren == -1 )
      maxnchildren = (psd->n)-1;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->nch), psd->n) );
   
   
   for( i = 0; i < psd->n; ++i )
   {

      /* variables created even if nchildrenvars=FALSE, since they are needed for the constraint 
         and they will be aggregated away
      */
      SCIPsnprintf(s, SCIP_MAXSTRLEN, "nchildren#%d", i);
      SCIP_CALL(SCIPcreateVar(scip, &(probdata->nch[i]), s, 0.0, maxnchildren, 0, SCIP_VARTYPE_INTEGER,
                              TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
      SCIP_CALL( SCIPaddVar(scip, probdata->nch[i]) );
      SCIPdebugMessage("adding variable %s with obj coefficient %f\n", SCIPvarGetName(probdata->nch[i]),
         SCIPvarGetObj(probdata->nch[i]));
      /* if !nchildrenvars then it is OK for them to be aggregated away */
      if( nchildrenvars )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip,probdata->nch[i]) ); 
         SCIP_CALL( SCIPchgVarBranchPriority(scip, probdata->nch[i], nchildrenvarpriority) ); 
      }
      
      SCIPsnprintf(s, SCIP_MAXSTRLEN, "nchildrencons#%d", i);
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL,
            0, 0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->nch[i], -1) );

      for( j = 0; j < psd->n; ++j )
      {
         if( j == i )
            continue;

         arrow_j_i = get_arrow(psd,j,i);
         if( arrow_j_i != NULL )
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, arrow_j_i, 1) );
      }

      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   return SCIP_OKAY;
}


/* Functions for creating the fundamental problem */
/** Creates and adds the parent variables to the problem.
 *
 *  @param scip The SCIP instance to create the variables for.
 *  @param psd The parent set data to base the variables on.
 *  @param scores The objective function coefficients for each variable.
 *  @return SCIP_OKAY if all variables were successfully created.
 */
static
SCIP_RETCODE addParentVariables(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_Real** scores
   )
{
   int i;
   int k;
   int l;
   char s[SCIP_MAXSTRLEN];
   char tmp[SCIP_MAXSTRLEN];

   SCIP_Bool pavarsinitial;
   SCIP_Bool pavarsremovable;
   
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/pavarsinitial", &pavarsinitial) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/pavarsremovable", &pavarsremovable) );

   
   SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->PaVars), psd->n) );
   for( i = 0; i < psd->n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->PaVars[i]), psd->nParentSets[i]) );

      /* sort, best parent set first for primal heuristic */
      for( k = 0; k < psd->nParentSets[i]; k++ )
         scores[i][k] = -scores[i][k];
      SCIPsortRealPtrPtrInt(scores[i], (void**)psd->PaVars[i], (void**)psd->ParentSets[i], psd->nParents[i], psd->nParentSets[i]);
      for( k = 0; k < psd->nParentSets[i]; ++k )
         scores[i][k] = -scores[i][k];

      for( k = 0; k < psd->nParentSets[i]; k++ )
      {
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "I#%s", psd->nodeNames[i]);
         for( l = 0; l < (psd->nParents[i][k]); l++ )
         {
            SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "#%s", psd->nodeNames[psd->ParentSets[i][k][l]]);
            strcat(s, tmp);
         }

         SCIP_CALL(SCIPcreateVar(scip, &(psd->PaVars[i][k]), s, 0.0, 1.0,
                                 scores[i][k],
                                 SCIP_VARTYPE_BINARY, pavarsinitial, pavarsremovable, NULL, NULL, NULL, NULL, NULL));

         SCIP_CALL( SCIPaddVar(scip, psd->PaVars[i][k]) );
         SCIPdebugMessage("adding variable %s with obj coefficient %f\n", SCIPvarGetName(psd->PaVars[i][k]), SCIPvarGetObj(psd->PaVars[i][k]));


      }

   }


   return SCIP_OKAY;
}

/** Adds the acyclity constraints to the problem.
 *
 *  @param scip The SCIP instance to add the constraints to.
 *  @param psd The parent set data to use for generating the acyclity constraints.
 *  @return SCIP_OKAY if the consraints were all added successfully.
 */
static
SCIP_RETCODE addAcyclicityConstraints(
   SCIP* scip,
   ParentSetData* psd
)
{
   SCIP_Bool useconsdagcluster;
   SCIP_Bool splitdags;
   SCIP_CONS* cons;
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/splitdags", &splitdags) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/useconsdagcluster", &useconsdagcluster) );

   if( !useconsdagcluster )
      return SCIP_OKAY;

   /* Only want to check for components if gobnilp/splitdags = TRUE */
   if( splitdags )
   {
      int i;
      int num_components;
      ParentSetData** components;
      SCIP_CALL( PS_splitToComponents(scip, psd, &num_components, &components) );
      for( i = 0; i < num_components; i++ )
      {
         if( components[i]->n > 2 )
         {
            SCIP_CALL( DC_createCons(scip, &cons, "DagCluster", components[i],  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
      for( i = 0; i < num_components; i++ )
         SCIP_CALL( PS_deallocateParentSetData(scip, &(components[i]), FALSE) );
      SCIPfreeMemoryArray(scip, &components);
   }
   else
   {
      SCIP_CALL( DC_createCons(scip, &cons, "DagCluster", psd,  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIPdebugMessage("Dagcluster constraint(s) added.\n");

   return SCIP_OKAY;
}

/** Adds the constraints that there must be one parent set to the problem.
 *
 *  @param scip The SCIP instance to add the constraints to.
 *  @param psd The parent set data to use for generating the constraints.
 *  @return SCIP_OKAY if the consraints were all added successfully.
 */
static
SCIP_RETCODE addOneParentSetConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   SCIP_CONS* cons;
   char s[SCIP_MAXSTRLEN];
   int i;
   SCIP_Bool convexconsinitial;
   
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/convexconsinitial", &convexconsinitial) );
   
   for( i = 0; i < psd->n; i++ )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "one_parent_set%d", i);
      SCIP_CALL( SCIPcreateConsSetpart(scip, &cons, s, psd->nParentSets[i], psd->PaVars[i], convexconsinitial, TRUE, TRUE, TRUE, TRUE,
            FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   return SCIP_OKAY;
}

/** Adds the variables and constraints needed for the basic DAG learning problem.
 *
 *  @param scip The SCIP instance to create the problem in.
 *  @param psd The parent set data to base the problem on.
 *  @param scores The scores of each of the parent sets in the problem.
 *  @return SCIP_OKAY if the problem was created successfully.
 */
static
SCIP_RETCODE createBasicModel(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_Real** scores
   )
{
   /* Generate family variables */
   SCIP_CALL( addParentVariables(scip, psd, scores) );

   /* add arrow and edge variables */
   SCIP_CALL( addArrowVariables(scip, psd) );

   /* add any total order  variables */
   SCIP_CALL( addAncestorVariables(scip, psd) );
   
   /* add any total order  variables */
   SCIP_CALL( addTotalorderVariables(scip, psd) );

   /* add any imset variables */
   SCIP_CALL( addImsetVariables(scip, psd) );

   /* add any posind variables */
   SCIP_CALL( addPosindVariables(scip, psd) );

   /* add any pos variables */
   SCIP_CALL( addPosVariables(scip, psd) );

   SCIP_CALL( addNChildrenVars(scip, psd) );

   SCIP_CALL( addPaSizeVariables(scip, psd) );

   /* Add any variables for k best MEC learning */
   SCIP_CALL( addKBestMECVariables(scip, psd) );

   /* Generate constraints */
   SCIP_CALL( addAcyclicityConstraints(scip, psd) );

   /* Add constraint that every node has one parent set */
   SCIP_CALL( addOneParentSetConstraints(scip, psd) );

   /* Set maximization */
   SCIP_CALL_ABORT(SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE));

   SCIPdebugMessage("Model created.\n");

   return SCIP_OKAY;
}



/* Functions for adding additional constraints */
/** Determines if two parent sets differ by just a single given item.
 *
 *  @param scip The SCIP instance the variable sets belong to.
 *  @param psd The parent set data for the problem.
 *  @param i The child variable.
 *  @param ki1 The index of one parent set of i.
 *  @param ki2 The index of another parent set of i.
 *  @param j A variable.
 *
 *  @return True if the parent sets ki1 and ki2 of i are the same except for
 *  the presence of variable j in one of them but not the other.
 */
static
SCIP_Bool differ(
   SCIP* scip,
   ParentSetData* psd,
   int i,
   int ki1,
   int ki2,
   int j
   )
{
   int big;
   int small;
   int l_big;
   int l_small;

   int nParents_diff;

   assert(psd != NULL);

   nParents_diff = psd->nParents[i][ki1] - psd->nParents[i][ki2];

   if( nParents_diff == 1 )
   {
      small = ki2;
      big = ki1;
   }
   else if( nParents_diff == -1 )
   {
      small = ki1;
      big = ki2;
   }
   else
      return FALSE;

   l_small = 0;
   for( l_big = 0;  l_big < psd->nParents[i][big]; ++l_big )
   {
      if( psd->ParentSets[i][big][l_big] == j )
         continue;

      if( l_small < psd->nParents[i][small] && psd->ParentSets[i][big][l_big] == psd->ParentSets[i][small][l_small] )
         l_small++;
      else
         return FALSE;
   }
   return TRUE;
}
/** Adds a constraint enforcing or preventing an immorality constraint.
 *
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data for the problem.
 *  @param i A parent involved in the immorality constraint.
 *  @param j The other parent involved.
 *  @param child The child involved in the constraint.
 *  @param truthvalue TRUE if the immorality must exist, FALSE if it must not exist.
 *
 *  @return SCIP_OK if the constraint was added, or an error code otherwise.
 */
static
SCIP_RETCODE immorality_constraint(
   SCIP* scip,
   ParentSetData* psd,
   int i,
   int j,
   int child,
   SCIP_Bool truthvalue
   )
{
   SCIP_CONS* cons;

   int k;
   int l;

   int found;

   char s[SCIP_MAXSTRLEN];

   SCIP_VAR** vars;
   int nvars = 0;
   SCIP_VAR* edge_i_j;
   SCIP_VAR* negedge_i_j;

   SCIP_Bool infeasible = FALSE;
   SCIP_Bool tightened;

   int tmp;
   int iter;

   assert(scip != NULL);
   assert(psd != NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, psd->nParentSets[child]+psd->nParentSets[i]+psd->nParentSets[j]+1) );

   for( k = 0; k < psd->nParentSets[child]; ++k )
   {
      found = 0;
      for( l = 0; l < psd->nParents[child][k]; ++l )
      {
         if( psd->ParentSets[child][k][l] == i || psd->ParentSets[child][k][l] == j )
            found++;

         if( found == 2 )
         {
            vars[nvars++] = psd->PaVars[child][k];
            break;
         }
      }
   }

   edge_i_j = get_edge(psd,i,j);
   if( truthvalue )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "user_immorality_yes#%d#%d#%d", child, i, j);
      /* for the immorality to occur one of the parent sets for child with both i and j must be chosen */
      SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &cons, s, nvars, vars) );
      /* also the parents must not be 'married'*/
      if( edge_i_j != NULL )
         SCIP_CALL( SCIPtightenVarUb(scip, edge_i_j, 0, TRUE, &infeasible, &tightened) ); 
      if( infeasible )
      {
         SCIPerrorMessage("Infeasibility detected while setting edge variable %s to 0 to enforce immorality %s->%s<-%s.\n", SCIPvarGetName(edge_i_j), 
            psd->nodeNames[i], psd->nodeNames[child], psd->nodeNames[j]);
         SCIPABORT();
         return SCIP_INVALIDDATA;
      }
      
   }
   else
   {
      if( edge_i_j != NULL )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, edge_i_j, &negedge_i_j) );
         vars[nvars++] = negedge_i_j;
      }

      /*tighted by adding in mutually exclusive vars which also imply the edge */

      for( iter = 0; iter < 2; ++iter )
      {
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            found = 0;
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               if( psd->ParentSets[i][k][l] == child || psd->ParentSets[i][k][l] == j )
                  found++;
               
               if( found == 2 )
               {
                  vars[nvars++] = psd->PaVars[i][k];
                  break;
               }
            }
         }
         tmp = i;
         i = j;
         j = tmp;
      }

      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "user_immorality_no#%d#%d#%d",  child, i, j);
      SCIP_CALL( SCIPcreateConsBasicSetpack(scip, &cons, s, nvars, vars) );
   }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIPfreeMemoryArray(scip, &vars);
   return SCIP_OKAY;
}

/** Adds a constraint stating that a given node must be a sink
 *  This is effected by fixing all parent set indicator variables to zero
 *  if the sink is one of the parents
 *  
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data for the problem.
 *  @param sink The node which must be a sink.
 *
 *  @return SCIP_OKAY if the constraint on the edge was added successfully or an error code otherwise.
 */
static
SCIP_RETCODE sink_constraint(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_VAR*** ancestorvars,
   int sink
   )
{
   
   int i;
   int k;
   int l;

   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed;

   for( i = 0; i < psd->n; ++i )
   {
      if( i == sink )
         continue;
    
      for( k = 0;  k < psd->nParentSets[i]; ++k )
          for( l = 0;  l < psd->nParents[i][k]; ++l )
             if( psd->ParentSets[i][k][l] == sink )
             {
                SCIP_CALL( SCIPtightenVarUb(scip, psd->PaVars[i][k], 0, TRUE, &infeasible, &fixed) ); 
                if( infeasible )
                {
                   SCIPerrorMessage("Infeasibility detected while forcing node %d to be a sink.\n",sink);
                   return SCIP_INVALIDDATA;
                }
                break;
             }
   }
   
   return SCIP_OKAY;
}

/** Adds a constraint stating that an edge must be present or absent.
 *  This works by simply fixing the relevant edge or arrow variable
 *  
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data for the problem.
 *  @param i The node from which the edge starts.
 *  @param j The node at which the edge finishes.
 *  @param undirected Whether the edge is undirected or not.
 *  @param truthvalue Whether the edge must or must not appear
 *
 *  @return SCIP_OKAY if the constraint on the edge was added successfully or an error code otherwise.
 */
static
SCIP_RETCODE edge_constraint(
   SCIP* scip,
   ParentSetData* psd,
   int i,
   int j,
   SCIP_Bool undirected,
   SCIP_Bool truthvalue
   )
{
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed;

   SCIP_VAR* var;

   assert(scip != NULL);
   assert(psd != NULL);
   assert(i != j);
   
   if( undirected )
      var = get_edge(psd,i,j);
   else
      var = get_arrow(psd,i,j);

   if( var != NULL )
   {
      if( truthvalue )
         SCIP_CALL( SCIPtightenVarLb(scip, var, 1, TRUE, &infeasible, &fixed) ); 
      else
         SCIP_CALL( SCIPtightenVarUb(scip, var, 0, TRUE, &infeasible, &fixed) ); 
   }

   if( infeasible || (var == NULL && truthvalue) )
   {
      SCIPerrorMessage("Infeasibility detected while setting arrow/edge variable %s to %d.\n", 
         var != NULL ? SCIPvarGetName(var) : "NULL", 
         truthvalue ? 1 : 0);
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }
   else
      return SCIP_OKAY;
}


/** Adds a constraint stating that an ancestor relation must be present or absent.
 *  This works by simply fixing the relevant ancestor variable
 *  
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data for the problem.
 *  @param i The descendant node
 *  @param j The ancestor node
 *  @param truthvalue Whether the relation must or must not appear
 *
 *  @return SCIP_OKAY if the constraint on the edge was added successfully or an error code otherwise.
 */
static
SCIP_RETCODE ancestor_constraint(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_VAR*** ancestorvars,
   int i,
   int j,
   SCIP_Bool truthvalue
   )
{
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed;

   SCIP_VAR* var;

   assert(scip != NULL);
   assert(psd != NULL);
   assert(i != j);

   if( ancestorvars == NULL)
   {
      SCIPerrorMessage("Cannot require presence or absence of ancestor relation without creating ancestor variables.\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }
   
   var = ancestorvars[i][j];
   
   if( truthvalue )
      SCIP_CALL( SCIPtightenVarLb(scip, var, 1, TRUE, &infeasible, &fixed) ); 
   else
      SCIP_CALL( SCIPtightenVarUb(scip, var, 0, TRUE, &infeasible, &fixed) ); 

   if( infeasible  )
   {
      SCIPerrorMessage("Infeasibility detected while setting ancestor variable %s to %d.\n", 
         SCIPvarGetName(var), 
         truthvalue ? 1 : 0);
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }
   else
      return SCIP_OKAY;
}


/** Adds a constraint stating that variable i must come before variable j
 *  in any variable ordering consistent with the learned network
 *  This works by simply fixing the relevant total order variable
 *  An error is generated if the user has not asked for the total order constraint
 *  handler to be used
 *  
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data for the problem.
 *  @param i The node which must be earlier in the order.
 *  @param j The node which must be later in the order.
 *
 *  @return SCIP_OKAY if the constraint on the edge was added successfully or an error code otherwise.
 */
static
SCIP_RETCODE order_constraint(
   SCIP* scip,
   ParentSetData* psd,
   int i,
   int j
   )
{
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed;
   
   SCIP_VAR* var;
   char* endptr;

   char s[SCIP_MAXSTRLEN];

   SCIP_Bool useconslinearordering;

   assert(scip != NULL);
   assert(psd != NULL);
   assert(i != j);

   /* get variable using its name which is a bit of a hack */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "<ord#%d#%d>", j, i);
   
   SCIP_CALL( SCIPparseVarName(scip, s, &var, &endptr) ); 

   if( var == NULL )
   {
      SCIPerrorMessage("Can't set an variable ordering constraint without also creating total order variables!\n Set gobnilp/totalordervars to TRUE\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }
   else
      SCIP_CALL( SCIPtightenVarLb(scip, var, 1, TRUE, &infeasible, &fixed) ); 


   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/useconslinearordering", &useconslinearordering) );

   if( !useconslinearordering )
   {
      SCIPerrorMessage("*** gobnilp/useconslinearordering not set to TRUE, so variable ordering constraint will not work proplerly.***\n");
   }
   
   if( infeasible )
   {
      SCIPerrorMessage("Infeasibility detected while setting order variable %s to TRUE.\n", 
         SCIPvarGetName(var));
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }
   else
      return SCIP_OKAY;
}



/** Parses a string to extract a set from it. */
static
SCIP_RETCODE parseSet(
   SCIP* scip,            /**< The SCIP instance that this belongs to. */
   ParentSetData* psd,    /**< The data set that the set refers to. */
   const char*  str,      /**< string to parse */
   int*  set,             /**< result set */
   int*  n_set,           /**< length of set */
   SCIP_Bool* success     /**< success flag */
)
{
   const char* t;
   char tmp[SCIP_MAXSTRLEN];
   char nodename[SCIP_MAXSTRLEN];
   int k;
   int nodeindex;

   *n_set = 0;
   t = str;
   while( *t != '\0' )
   {
      k = 0;
      while( *t != '\0' && *t != ',' )
         tmp[k++] = *t++;
      tmp[k] = '\0';
      if( sscanf(tmp, "%d", &nodeindex) != 1 )
      {
         if( sscanf(tmp, "%s", nodename) != 1 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Expected variable name. Got: %s\n", tmp);
            *success = FALSE;
            return SCIP_OKAY;
         }
         nodeindex = get_index(nodename, psd);
         if( nodeindex == -1 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Expected variable name. Got: %s\n", tmp);
            *success = FALSE;
            return SCIP_OKAY;
         }
      }

      set[(*n_set)++] = nodeindex;

      if( *t == ',' )
         t++;
   }
   return SCIP_OKAY;
}
/** Adds a conditional independence (ci) constraint.
 *
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data for the problem.
 *  @param a_str A set
 *  @param b_str B set
 *  @param s_str S (separator) set
 *
 *  @return SCIP_OKAY if the constraint was added successfully or an error code otherwise.
 */
static
SCIP_RETCODE ci_constraint(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_VAR*** ancestorvars,
   const char* a_str,
   const char* b_str,
   const char* s_str
   )
{
   SCIP_CONS* cons;

   int a[SCIP_MAXSTRLEN];
   int b[SCIP_MAXSTRLEN];
   int sep[SCIP_MAXSTRLEN];
   int n_a;
   int n_b;
   int n_s;

   SCIP_Bool success;
   char consname[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(psd != NULL);
   assert(a_str != NULL);
   assert(b_str != NULL);
   assert(s_str != NULL);
   
   success = TRUE;
   SCIP_CALL( parseSet(scip, psd, a_str, a, &n_a, &success) );
   if( !success )
      return SCIP_ERROR;
   SCIP_CALL( parseSet(scip, psd, b_str, b, &n_b, &success) );
   if( !success )
      return SCIP_ERROR;
   SCIP_CALL( parseSet(scip, psd, s_str, sep, &n_s, &success) );
   if( !success )
      return SCIP_ERROR;



   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s_|_%s|%s", a_str, b_str, s_str);

   /* generate ci constraint */
   SCIP_CALL(SCIPcreateConsCi(
                scip,
                &cons,
                consname,
                psd,
                ancestorvars,
                a,
                n_a,
                b,
                n_b,
                sep,
                n_s,
                TRUE,
                TRUE,
                TRUE,
                TRUE,
                TRUE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE
             ));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;

}
/** Adds a constraint on the DAG structure.
 *
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param psd The parent set data relating to this problem.
 *  @param line The description of the constraint to add.
 *
 *  @return SCIP_OKAY if the constraint was added or an error otherwise.
 */
static
SCIP_RETCODE process_constraint(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_VAR*** ancestorvars,
   const char* line
   )
{

   int i;
   int j;
   int child;

   char a_str[SCIP_MAXSTRLEN];
   char b_str[SCIP_MAXSTRLEN];
   char s_str[SCIP_MAXSTRLEN];

   if( line[0] == '#' )
      return SCIP_OKAY;

   if( sscanf(line, "%[^~<>-]-%[^~<>-]", a_str, b_str) == 2 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( edge_constraint(scip, psd, i, j, TRUE, TRUE) );
   }
   else if( sscanf(line, "~%[^~<>-]-%[^~<>-]", a_str, b_str) == 2 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( edge_constraint(scip, psd, i, j, TRUE, FALSE) );
   }
   else if( sscanf(line, "%[^~<>-]<-%[^~<>-]", a_str, b_str) == 2 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( edge_constraint(scip, psd, i, j, FALSE, TRUE) );
   }
   else if( sscanf(line, "~%[^~<>-]<-%[^~<>-]", a_str, b_str) == 2 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( edge_constraint(scip, psd, i, j, FALSE, FALSE) );
   }
   else if( sscanf(line, "%[^~<>-]<--%[^~<>-]", a_str, b_str) == 2 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( ancestor_constraint(scip, psd, ancestorvars, i, j, TRUE) );
   }
   else if( sscanf(line, "~%[^~<>-]<--%[^~<>-]", a_str, b_str) == 2 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( ancestor_constraint(scip, psd, ancestorvars, i, j, FALSE) );
   }
   else if( sscanf(line, "%[^~<>-]->%[^~<>-]<-%[^~<>-]", a_str, s_str, b_str) == 3 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      child = get_index(s_str, psd);
      if( i == -1 || j == -1 || child == -1 )
         return SCIP_READERROR;

      SCIP_CALL( immorality_constraint(scip, psd, i, j, child, TRUE) );
   }
   else if( sscanf(line, "~%[^~<>-]->%[^~<>-]<-%[^~<>-]", a_str, s_str, b_str) == 3 )
   {
      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      child = get_index(s_str, psd);
      if( i == -1 || j == -1 || child == -1 )
         return SCIP_READERROR;

      SCIP_CALL( immorality_constraint(scip, psd, i, j, child, FALSE) );
   }
   else if( sscanf(line, "%[^_~<>-]_|_%[^|~<>-]|%[^~<>-]", a_str, b_str, s_str) == 3 )
      ci_constraint(scip, psd, ancestorvars, a_str, b_str, s_str);
   else if( sscanf(line, "%[^_~<>-]_|_%[^~<>-]", a_str, b_str) == 2 )
      ci_constraint(scip, psd, ancestorvars, a_str, b_str, "");
   else if( sscanf(line, "%[^~<>-]<%[^~<>-]", a_str, b_str) == 2 )
   {

      i = get_index(a_str, psd);
      j = get_index(b_str, psd);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( order_constraint(scip, psd, i, j) );
   }
   else if ( sscanf(line, "sink%s", a_str) == 1 )
   {
      i = get_index(a_str, psd);
      if( i == -1 )
         return SCIP_READERROR;
      SCIP_CALL( sink_constraint(scip, psd, ancestorvars, i) );
   }
   else
   {
      SCIPerrorMessage("Not recognised as a DAG constraint: %s\n", line);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}
/** Adds constraints ruling out covered arcs from lower to higher vertices.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addOrderedCoveredArcConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   int i;
   int l;
   int small_i;
   int big_i;
   int small_j;
   int big_j;
   int ki1;
   int ki2;
   int kj1;
   int kj2;
   int j;
   int l2;
   SCIP_Bool ok2;
   SCIP_VAR* arc_tmp[3];
   char s[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;

   SCIP_VAR* arrow_var;
   int nvars;
   
   /*

   rule out covered arcs from higher to lower vertex

   for each i,j,C (j>i,C a set) such that the following variables exist:
   I(i<-C), I(i<-C+j), I(j<-C), I(j<-C+i)

   if I(i<-C+j)=1,I(j<-C)=1
   then the edge i<-j is 'covered' and
   can be replaced by:
   I(i<-C)=1,I(j<-C+i)
   which reverse the edge


   without creating a cycle or changing Markov equivalence
   class, so add constraint:
   I(i<-C+j) + I(j<-C) <= 1
   to rule out first case

   this can be strengthened to
   I(i<-C+j) + I(j<-C) + arrow(j<-i) <= 1
   */

   for( i = 0; i < psd->n; i++ )
   {
      for( ki1 = 0; ki1 < psd->nParentSets[i]; ki1++ )
      {
         for( ki2 = ki1 + 1; ki2 < psd->nParentSets[i]; ki2++ )
         {
            if( psd->nParents[i][ki1] - psd->nParents[i][ki2] == 1 )
            {
               big_i = ki1;
               small_i = ki2;
            }
            else if( psd->nParents[i][ki2] - psd->nParents[i][ki1] == 1 )
            {
               big_i = ki2;
               small_i = ki1;
            }
            else
               /* need to at least find two parent sets for i differing in size by 1 */
               continue;

            for( j = i + 1; j < psd->n; j++ )
            {
               /* ki1 and ki2 can only differ by j */
               if( !differ(scip, psd, i, ki1, ki2, j) )
                  continue;

               for( kj1 = 0; kj1 < psd->nParentSets[j]; kj1++ )
               {
                  for( kj2 = kj1 + 1; kj2 < psd->nParentSets[j]; kj2++ )
                  {
                     if( psd->nParents[j][kj1] - psd->nParents[j][kj2] == 1 )
                     {
                        small_j = kj2;
                        big_j = kj1;
                     }
                     else if( psd->nParents[j][kj2] - psd->nParents[j][kj1] == 1 )
                     {
                        small_j = kj1;
                        big_j = kj2;
                     }
                     else
                        /* need to at least find two parent sets for j differing in size by 1 */
                        continue;

                     if( !differ(scip, psd, j, kj1, kj2, i) )
                        continue;

                     /* small_i and small_j  must be the same, if so big_i is small_j with j added */
                     if( psd->nParents[i][small_i] !=  psd->nParents[j][small_j] )
                        continue;

                     l2 = 0;
                     ok2 = TRUE;
                     for( l = 0; l < psd->nParents[i][small_i]; l++ )
                     {
                        if( psd->ParentSets[i][small_i][l] != psd->ParentSets[j][small_j][l2] )
                        {
                           ok2 = FALSE;
                           break;
                        }
                        else
                           l2++;
                     }

                     if( !ok2 )
                        continue;

                     SCIPdebugMessage("Ruling out having both %s and %s since arc %d<-%d is covered and there exists %s and %s\n",
                                      SCIPvarGetName(psd->PaVars[i][big_i]), SCIPvarGetName(psd->PaVars[j][small_j]), i, j,
                                      SCIPvarGetName(psd->PaVars[i][small_i]), SCIPvarGetName(psd->PaVars[j][big_j]));
                     arc_tmp[0] = psd->PaVars[i][big_i];
                     arc_tmp[1] = psd->PaVars[j][small_j];
                     arrow_var = get_arrow(psd,j,i);
                     if( arrow_var != NULL )
                     {
                        arc_tmp[2] = arrow_var;
                        nvars = 3;
                     }
                     else
                        nvars = 2;
                     SCIPsnprintf(s, SCIP_MAXSTRLEN, "covered_arc#%s#%s", psd->nodeNames[i], psd->nodeNames[j]);
                     SCIP_CALL(SCIPcreateConsSetpack(scip, &cons, s, nvars, arc_tmp,
                                                     FALSE,  /* initial */
                                                     TRUE,
                                                     FALSE, /* enforce */
                                                     FALSE, /* check */
                                                     TRUE, /* propagate */
                                                     FALSE, FALSE, FALSE, FALSE, FALSE));
                     SCIP_CALL( SCIPaddCons(scip, cons) );
                     SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                  }
               }
            }
         }
      }
   }
   return SCIP_OKAY;
}



/** Adds constraints on the number of possible founders.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addFounderConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   int i;
   int k;
   int minfounders;
   int maxfounders;
   SCIP_CONS* cons;

   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/minfounders", &minfounders) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/maxfounders", &maxfounders) );

   if( minfounders == 0 && maxfounders == -1 )
      return SCIP_OKAY;

   SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "founders", 0, NULL, NULL, minfounders,
                                  maxfounders > 0 ? maxfounders : SCIPinfinity(scip),
                                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
         if( psd->nParents[i][k] == 0 )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][k], 1) );
            break;
         }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   return SCIP_OKAY;
}





/** Adds constraints on the number of nodes that can be parents.
 *
 *  This also requires the creation of a number of binary variables.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addParentNumberConstraints(
   SCIP* scip,
   const ParentSetData* psd
   )
{
   char s[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   int i;        /* indexes variables */
   int j;        /* indexes variables */
   int minparents;
   int maxparents;
   SCIP_VAR** ispas;
   int n_ispas;
   SCIP_Real* ones;
   SCIP_VAR* arrow_j_i;

   SCIP_PROBDATA* probdata;
   
   assert(scip != NULL);
   assert(psd != NULL);

   SCIPgetIntParam(scip, "gobnilp/minparents", &minparents);
   SCIPgetIntParam(scip, "gobnilp/maxparents", &maxparents);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->ispa, psd->n) );

   SCIP_CALL( SCIPallocBufferArray(scip, &ispas, psd->n) );


   for( i = 0; i < psd->n; ++i )
   {

      SCIPsnprintf(s, SCIP_MAXSTRLEN, "ispa#%d", i);
      SCIP_CALL(SCIPcreateVar(scip, &(probdata->ispa[i]), s, 0.0, 1.0, 0, SCIP_VARTYPE_BINARY,
                              TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
      SCIP_CALL( SCIPaddVar(scip, probdata->ispa[i]) );
      SCIPdebugMessage("adding variable %s with obj coefficient %f\n", SCIPvarGetName(probdata->ispa[i]), SCIPvarGetObj(probdata->ispa[i]));

      n_ispas=0;
      for( j = 0; j < psd->n; ++j )
      {
         if( j == i )
            continue;

         arrow_j_i = get_arrow(psd,j,i);
         if( arrow_j_i != NULL )
            ispas[n_ispas++] = arrow_j_i;
      }
      
      assert( n_ispas == (psd->n)-1 );
      
      SCIPsnprintf(s, SCIP_MAXSTRLEN, "ispacons#%s", psd->nodeNames[i]);
      SCIP_CALL(SCIPcreateConsOr(scip, &cons, s, probdata->ispa[i], n_ispas, ispas,
                                 TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIPfreeBufferArray(scip, &ispas);

   SCIP_CALL( SCIPallocBufferArray(scip, &ones, psd->n) );
   for( i = 0; i < psd->n; ++i )
      ones[i] = 1.0;

   SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "parentscons", psd->n, probdata->ispa, ones,
                                  minparents, maxparents == -1 ? SCIPinfinity(scip) : maxparents,
                                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIPfreeBufferArray(scip, &ones);

   return SCIP_OKAY;
}
/** Adds constraints on the number of edges in the DAG.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addEdgeNumberConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   SCIP_CONS* cons;
   int minedges, maxedges;
   int i;
   int j;
   SCIP_VAR* edge_i_j;

   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/minedges", &minedges) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/maxedges", &maxedges) );


   if( minedges != 0 || maxedges != -1 )
   {
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "edges", 0, NULL, NULL,
                                     minedges, maxedges > 0 ? maxedges : SCIPinfinity(scip),
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      for( i = 0; i < psd->n; ++i )
      {
         for( j = i+1; j < psd->n; ++j )
         {
            edge_i_j = get_edge(psd,i,j);
            if( edge_i_j != NULL )
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, edge_i_j, 1) );
         }
      }
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}

/** Adds general structural constraints to the DAG.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addGeneralDAGConstraints(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_VAR*** ancestorvars
   )
{
   int status;
   char s[SCIP_MAXSTRLEN];
   char* dagconstraintsfile;
   FILE* dagconstraints;

   SCIPgetStringParam(scip, "gobnilp/dagconstraintsfile", &dagconstraintsfile);
   dagconstraints = fopen(dagconstraintsfile, "r");
   if( dagconstraints == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", dagconstraintsfile);
      return SCIP_NOFILE;
   }

   status = fscanf(dagconstraints, "%[^\n]%*c", s);
   while( status == 1 )
   {
      SCIP_CALL( process_constraint(scip, psd, ancestorvars, s) );
      status = fscanf(dagconstraints, "%[^\n]%*c", s);
   }

   fclose(dagconstraints);
   return SCIP_OKAY;
}

/** Adds constraints ruling out immoralities in the DAG.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
#if 0
static
SCIP_RETCODE addNoImmoralityConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   int i;
   int j;
   int k;
   int l;

   int i0;
   int i1;
   int i2;
   int i3;

   SCIP_Bool feasible;
   SCIP_Bool fixed;

   SCIP_Bool decomposable;
   SCIP_Bool ***store;

   int nodes[4];
   SCIP_VAR** vars;
   SCIP_Real* ones;
   int nvars;
   int mostparentsets;
   
   SCIP_Bool ok;
   int jj;
   SCIP_CONS* cons;
   SCIP_VAR* edge;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/decomposable", &decomposable) );
   
   /* if learning a decomposable model then can arbitrarily fix first variable to be an orphan
      BUT ONLY IF WORKING WITH LIKELIHOOD EQUIVALENCE
   */
   if( decomposable )
   {
      for( k = 0; k < psd->nParentSets[0]; ++k )
         {
            if( psd->nParents[0][k] == 0 )
               SCIP_CALL( SCIPfixVar(scip, psd->PaVars[0][k], 1, &feasible, &fixed) );
            else
               SCIP_CALL( SCIPfixVar(scip, psd->PaVars[0][k], 0, &feasible, &fixed) );
         }
      
   }

   /* create bitmap for faster checking */
   mostparentsets = 0;
   SCIP_CALL( SCIPallocMemoryArray(scip, &store, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      if( psd->nParentSets[i] > mostparentsets)
         mostparentsets = psd->nParentSets[i];

      SCIP_CALL( SCIPallocMemoryArray(scip, &(store[i]), psd->nParentSets[i]) );
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(store[i][k]), psd->n) );
         for( l = 0; l < psd->n; ++l )
            store[i][k][l] = FALSE;
         for( l = 0; l < psd->nParents[i][k]; ++l )
            store[i][k][psd->ParentSets[i][k][l]] = TRUE;
      }
   }
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, 3*mostparentsets ) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ones, 3*mostparentsets ) );
   for( i = 0; i < 3*mostparentsets; ++i )
      ones[i] = 1;

   /* add monotonicity constraints for subsets of size 4 and 3 */
   for( i0 = 0; i0 < psd->n; ++i0 )
   {
      
      nodes[0] = i0;
      for( i1 = i0+1; i1 < psd->n; ++i1 )
      {
         nodes[1] = i1;
         for( i2 = i1+1; i2 < psd->n; ++i2 )
         {
            nodes[2] = i2;

            /* include variables where a node has other two as parents */
            nvars = 0;
            for( j = 0; j < 3; ++j )
            {
               for( k = 0; k < psd->nParentSets[nodes[j]]; ++k )
               {
                  ok = TRUE;
                  for( jj = 0; jj < 3; ++jj )
                  {
                     if( jj != j && !store[nodes[j]][k][nodes[jj]] )
                     {
                        ok = FALSE;
                        break;
                     }
                  }
                  if( ok )
                     vars[nvars++] = psd->PaVars[nodes[j]][k];
               }
            }

            /* include variables where a node has one of other two as parents 
               i.e. the edge variable 
            */
            if( nvars > 0 )
            {
               for( j = 0; j < 2; ++j )
               {
                  for( jj = j+1; jj < 3; ++jj )
                  {
                     SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "todo", nvars, vars, ones, -SCIPinfinity(scip), 0,
                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
                     
                     edge = get_edge(psd,nodes[j],nodes[jj]);
                     if( edge != NULL )
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, edge, -1) );
                     
                     SCIP_CALL( SCIPaddCons(scip, cons) );
                     /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
                     SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                  }
               }
            }


            
            for( i3 = 0; i3 < psd->n; ++i3 )
            {
               if( i3 == i0 || i3 == i1 || i3 == i2 )
                  continue;

               nodes[3] = i3;

               SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "todo", nvars, vars, ones, 0, SCIPinfinity(scip),
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
               
               /* include variables where a node has other three as parents */
               for( j = 0; j < 4; ++j )
               {
                  for( k = 0; k < psd->nParentSets[nodes[j]]; ++k )
                  {
                     ok = TRUE;
                     for( jj = 0; jj < 4; ++jj )
                     {
                        if( jj != j && !store[nodes[j]][k][nodes[jj]] )
                        {
                           ok = FALSE;
                           break;
                        }
                     }
                     if( ok )
                        SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[nodes[j]][k], -1) );
                  }
               }
               
               SCIP_CALL( SCIPaddCons(scip, cons) );
               /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }
      }
   }

   for( i = 0; i < psd->n; ++i )
   {
      for( k = 0; k < psd->nParentSets[i]; ++k )
         SCIPfreeMemoryArray(scip, &(store[i][k]));

      SCIPfreeMemoryArray(scip, &(store[i]));
   }

   SCIPfreeMemoryArray(scip, &store);
   SCIPfreeMemoryArray(scip, &vars);
   SCIPfreeMemoryArray(scip, &ones);

   for( i = 0; i < psd->n; ++i )
   {
   /* add in constraint for single parent sets */
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "todo", 0, NULL, NULL, 0, SCIPinfinity(scip),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      
      for( k = 0; k < psd->nParentSets[i]; ++k )
         if( psd->nParents[i][k] > 1 )
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][k], -1) );

      for( j = 0; j < psd->n; ++j )
      {
         if( j == i)
            continue;

         for( k = 0; k < psd->nParentSets[j]; ++k )
            if( psd->nParents[j][k] == 1 && psd->ParentSets[j][k][0] != i)
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[j][k], 1) );
      }         
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }
   
   return SCIP_OKAY;
}
#endif

/** Adds optional constraints to the problem based on the parameters in the settings file.
 *
 *  @param scip The SCIP instance to add the constraint to.
 *  @param psd The parent set data associated with this problem.
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addAdditionalConstraints(
   SCIP* scip,
   ParentSetData* psd
   )
{
   SCIP_Bool noimmoralities;
   SCIP_Bool decomposable;
   SCIP_Bool orderedcoveredarcs;
   int minparents, maxparents;
   char* dagconstraintsfile;
   SCIP_Bool distinguishedrep;
   SCIP_Bool vanilla;

   SCIP_PROBDATA* probdata;
      
   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   
   /* Constraints on ordered covered arcs */
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/orderedcoveredarcs", &orderedcoveredarcs) );
   if( orderedcoveredarcs )
      SCIP_CALL( addOrderedCoveredArcConstraints(scip, psd) );

   /* Constraints on number of founders */
   SCIP_CALL( addFounderConstraints(scip, psd) );

   /* Constraints on number of parents */
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/minparents", &minparents) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/maxparents", &maxparents) );
   if( minparents > 0 || (maxparents != -1 && maxparents < psd->n) )
      SCIP_CALL( addParentNumberConstraints(scip, psd) );

   /* Constraints on number of edges */
   SCIP_CALL( addEdgeNumberConstraints(scip, psd) );

   /* Constraints on immoralities */
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/noimmoralities", &noimmoralities) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/decomposable", &decomposable) );
   if( decomposable || noimmoralities )
      /* SCIP_CALL( addNoImmoralityConstraints(scip, psd) ); */
   {
      SCIP_CONS* cons;
      SCIP_CALL( SCIPcreateConsBasicChordal(scip, &cons, "chordal", psd) ); 
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }      
   /* Arbitrary DAG Constraints */
   SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/dagconstraintsfile", &dagconstraintsfile) );
   if( strcmp(dagconstraintsfile, "") != 0 )
      SCIP_CALL( addGeneralDAGConstraints(scip, psd, probdata->ancestorvars) );

   /* Pedigree constraints */
   if( PD_inPedigreeMode(scip) )
      SCIP_CALL( PD_addPedigreeSpecificConstraints(scip, psd) );

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/distinguishedrep", &distinguishedrep) );
   if( distinguishedrep )
   {
      SCIP_CONS* cons;
      SCIP_CALL( SCIPincludeConshdlrDistinguishedrep(scip) );
      SCIP_CALL( SCIPcreateConsBasicDistinguishedrep(scip, &cons, "distrep", psd) ); 
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/vanilla", &vanilla) );
   if( vanilla )
   {
      SCIP_CONS* cons;

      SCIP_CALL( SCIPincludeConshdlrVanilla(scip) );
      /* SCIP_CALL( SCIPcreateConsBasicVanilla(scip, &cons, "vanilla", psd, probdata->ancestorvars) ); */
      SCIP_CALL( SCIPcreateConsVanilla(scip, &cons, "vanilla", psd, probdata->ancestorvars,
            FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   return SCIP_OKAY;
}


/** Sets the branching priority of variables. */
static
SCIP_RETCODE setBranchingPriorities(
   SCIP* scip,                /**< The SCIP instance being used. */
   ParentSetData* psd         /**< The data for the problem.     */
   )
{
   int i;
   int j;
   int pavarspriority;
   int edgevarpriority;
   int arrowvarpriority;
   SCIP_VAR* arrow_i_j;
   SCIP_VAR* edge_i_j;

   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/pavarspriority", &pavarspriority) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/edgevarpriority",&edgevarpriority) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/arrowvarpriority",&arrowvarpriority) );

   for( i = 0; i < psd->n; i++ )
      for( j = 0; j < psd->n; j++ )
         if ( j != i )
         {
            arrow_i_j = get_arrow(psd,i,j);
            if( arrow_i_j != NULL )
               SCIP_CALL( SCIPchgVarBranchPriority(scip, arrow_i_j, arrowvarpriority) );
         }
   
   for( i = 0; i < psd->n; i++ )
      for( j = i+1; j < psd->n; j++ )
      {
         edge_i_j = get_edge(psd,i,j);
         if( edge_i_j != NULL )
            SCIP_CALL( SCIPchgVarBranchPriority(scip, edge_i_j, edgevarpriority) );
      }

   for( i = 0; i < psd->n; i++ )
      for( j = 0; j < psd->nParentSets[i]; j++ )
      {
         assert( psd->PaVars[i][j] != NULL );
         SCIP_CALL( SCIPchgVarBranchPriority(scip, psd->PaVars[i][j], pavarspriority) );
      }

   return SCIP_OKAY;
}

/** Sets properties to non-standard values. */
static
SCIP_RETCODE setSpecialProperties(
   SCIP* scip,                /**< The SCIP instance being used. */
   ParentSetData* psd         /**< The data for the problem.     */
   )
{
   int i, j;
   SCIP_VAR* arrow_i_j;
   SCIP_VAR* edge_i_j;

   return SCIP_OKAY;
   
   for (i = 0; i < psd->n; i++)
      for( j = 0; j < psd->nParentSets[i]; j++ )
      {
         assert( psd->PaVars[i][j] != NULL );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip,psd->PaVars[i][j]) );
      }

   for (i = 0; i < psd->n; i++)
      for (j = 0; j < psd->n; j++)
         if ( j != i )
         {
            arrow_i_j = get_arrow(psd,i,j);
            if( arrow_i_j != NULL )
               SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip,arrow_i_j) );
         }

   for (i = 0; i < psd->n; i++)
      for (j = i+1; j < psd->n; j++)
      {
         edge_i_j = get_edge(psd,i,j);
         if( edge_i_j != NULL )
            SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip,edge_i_j) );
      }
   
   return SCIP_OKAY;
}


/* Functions for reading in the problem */
/** Gets the delimiter information needed to read in information from files.
 *
 *  @param scip The SCIP instance being used.
 *  @param num_delims The number of delimiters.
 *  @param delims The delimiters to use.
 *  @param merge_delims Whether consecutive delimiters should be merged into one.
 *  @return SCIP_OKAY if the parameters were successfully read, or an error otherwise.
 */
static
SCIP_RETCODE getDelimiters(
   SCIP* scip,
   int* num_delims,
   char** delims,
   SCIP_Bool* merge_delims
   )
{
   char* delim_string;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/mergedelimiters", merge_delims) );
   SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/delimiter", &delim_string) );
   if( strcmp(delim_string, "whitespace") == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, delims, 2) );
      (*delims)[0] = ' ';
      (*delims)[1] = '\t';
      *num_delims = 2;
   }
   else if( strcmp(delim_string, "tab") == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, delims, 1) );
      (*delims)[0] = '\t';
      *num_delims = 1;
   }
   else if( strlen(delim_string) > 1 )
   {
      SCIPerrorMessage("Value of \"gobnilp/delimiter\" is incorrect.  It must be either a single character, \"whitespace\" or \"tab\".\n");
      return SCIP_ERROR;
   }
   else
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, delims, 1) );
      (*delims)[0] = delim_string[0];
      *num_delims = 1;
   }
   return SCIP_OKAY;
}


static
int countVars(
   FILE* file
   )
{
   int num_vars = 0;
   char this_line[SCIP_MAXSTRLEN];
   SCIP_Bool end = FALSE;

   while( end == FALSE )
   {
      char* retcode = fgets(this_line, SCIP_MAXSTRLEN, file);
      if( retcode == NULL )
         end = TRUE;
      if( strncmp(this_line, "VAR", 3) == 0 )
         num_vars += 1;
   }

   return num_vars;
}

static
int countParentSets(
   FILE* file
   )
{
   int start_pos = ftell(file);
   int num_scores = 0;
   char this_line[SCIP_MAXSTRLEN];
   SCIP_Bool end = FALSE;
   SCIP_Bool next_var_found = FALSE;

   while( end == FALSE && next_var_found == FALSE )
   {
      char* retcode = fgets(this_line, SCIP_MAXSTRLEN, file);
      int line_length = strlen(this_line);
      if( retcode == NULL )
         end = TRUE;
      if( strncmp(this_line, "META", 4) == 0 )
         ; /* Meta Line */
      else if( strncmp(this_line, "#", 1) == 0 )
         ; /* Comment */
      else if( line_length == 1 || line_length == 0 )
         ; /* Blank Line */
      else if( line_length == -1 )
         ; /* End of File */
      else if( strncmp(this_line, "VAR", 3) == 0 )
         next_var_found = TRUE;
      else
         num_scores += 1;
   }

   fseek(file, start_pos, SEEK_SET);
   return num_scores;
}

static
int splitAtWhitespace(
   const char* string,
   char*** parts
   )
{
   int i;

   int part_length = 0;
   char* this_part = malloc((strlen(string) + 1) * sizeof(char));

   int num_parts = 0;
   (*parts) = malloc((strlen(string) + 1) * sizeof(char*));

   for( i = 0; i < (int)strlen(string); i++ )
   {
      char this_char = string[i];
      if( isspace((int)this_char) )
      {
         if( part_length > 0 )
         {
            this_part[part_length] = '\0';
            this_part = realloc(this_part, (part_length + 1) * sizeof(char));
            (*parts)[num_parts] = this_part;
            num_parts += 1;
            this_part = malloc((strlen(string) + 1) * sizeof(char));
            part_length = 0;
         }
      }
      else
      {
         this_part[part_length] = this_char;
         part_length += 1;
      }
   }

   if( part_length > 0 )
   {
      this_part[part_length] = '\0';
      this_part = realloc(this_part, (part_length + 1) * sizeof(char));
      (*parts)[num_parts] = this_part;
      num_parts += 1;
   }
   else
   {
      free(this_part);
   }
   (*parts) = realloc((*parts), num_parts * sizeof(char*));

   return num_parts;
}

static
void extractMetaNameAndValue(
   const char* string,
   char* name,
   char* value
   )
{
   int name_length = 0;
   int value_length = 0;
   int pos = 4;
   char this_char = string[4];

   /* Jump any initial space */
   while( isspace((int)this_char) )
   {
      pos += 1;
      this_char = string[pos];
   }

   /* Collect the name */
   while( !(isspace((int)this_char)) && (this_char != '=') )
   {
      name[name_length] = this_char;
      name_length += 1;
      pos += 1;
      this_char = string[pos];
   }
   name[name_length] = '\0';

   /* Jump any space up to the = */
   while( this_char != '=' )
   {
      pos += 1;
      this_char = string[pos];
   }

   /* Jump the = */
   pos += 1;
   this_char = string[pos];

   /* Jump any space */
   while( isspace((int)this_char) )
   {
      pos += 1;
      this_char = string[pos];
   }

   /* Collect the value */
   while( this_char != '\0' && this_char != '\r' && this_char != '\n' )
   {
      value[value_length] = this_char;
      value_length += 1;
      pos += 1;
      this_char = string[pos];
   }
   value[value_length] = '\0';
}

/** Reads local score BN file (in PSS format)
 *
 *  @param scip The SCIP instance being used.
 *  @param filename The file to read from.
 *  @param psd The parentage data that is read from the file.
 *  @param scores The score of each of the parent set combinations.
 *  @param prop The properties of the data read from the file.
 *  @return SCIP_OKAY if reading was successful, or an appropriate error code otherwise.
 */
static
SCIP_RETCODE readProblemInPSSFormat(
   SCIP* scip,
   const char* filename,
   ParentSetData* psd,
   SCIP_Real*** scores,
   PropertyData* prop
   )
{
   int i;
   int k;
   int l;
   int m;
   FILE *file;
   int n;
   SCIP_Real** Scores;
   int* nParentSets;
   int** nParents;
   int*** ParentSets;
   char** nodeNames;
   char**** tmp_ParentSets;
   int numCandidateParentSets = 0;
   int current_var = -1;
   int current_parent_set = -1;
   SCIP_Bool end = FALSE;

   /* open file */
   if( strcmp(filename, "-") == 0 )
      file = stdin;
   else
      file = fopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", filename);
      return SCIP_NOFILE;
   }

   n = countVars(file);
   rewind(file);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of variables: %d\n", n);
   psd->n = n;
   prop->n = n;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prop->num_properties), prop->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prop->property_names), prop->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prop->property_values), prop->n) );
   for( i = 0; i < prop->n; i++ )
      prop->num_properties[i] = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &Scores, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParentSets, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParents, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ParentSets, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_ParentSets, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeNames, n) );
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nodeNames[i]), SCIP_MAXSTRLEN) );

   *scores = Scores;
   psd->nParentSets = nParentSets;
   psd->nParents = nParents;
   psd->ParentSets = ParentSets;
   psd->nodeNames = nodeNames;

   while( end == FALSE )
   {
      char this_line[SCIP_MAXSTRLEN];
      char* retcode = fgets(this_line, SCIP_MAXSTRLEN, file);
      int line_length = strlen(this_line);
      if( retcode == NULL )
         end = TRUE;

      if( line_length == -1 )
         ; /* End of File */
      else if( line_length == 1 || line_length == 0 )
         ; /* Blank Line */
      else if( strncmp(this_line, "#", 1) == 0 )
         ; /* Comment */
      else if( strncmp(this_line, "META", 4) == 0 )
      {
         /* Meta Line */
         char name[SCIP_MAXSTRLEN];
         char value[SCIP_MAXSTRLEN];
         extractMetaNameAndValue(this_line, name, value);
         if( current_var == -1 )
         {
            if( PR_hasGlobalProperty(scip, prop, name) )
               SCIPwarningMessage(scip, "Repeated global property %s in %s.\n", name, filename);
            PR_setGlobalProperty(scip, prop, name, value);
         }
         else
         {
            if( PR_hasProperty(scip, prop, current_var, name) )
               SCIPwarningMessage(scip, "Repeated property %s for %s in %s.\n", name, nodeNames[current_var], filename);
            PR_setProperty(scip, prop, current_var, name, value);
         }
      }
      else if( strncmp(this_line, "VAR", 3) == 0 )
      {
         /* VAR line */
         current_var += 1;
         current_parent_set = -1;
         sscanf(this_line, "VAR %s", nodeNames[current_var]);
         nParentSets[current_var] = countParentSets(file);
         numCandidateParentSets += nParentSets[current_var];
         SCIP_CALL( SCIPallocMemoryArray(SCIP, &(nParents[current_var]), nParentSets[current_var]) );
         SCIP_CALL( SCIPallocMemoryArray(SCIP, &(Scores[current_var]), nParentSets[current_var]) );
         SCIP_CALL( SCIPallocMemoryArray(SCIP, &(ParentSets[current_var]), nParentSets[current_var]) );
         SCIP_CALL( SCIPallocMemoryArray(SCIP, &(tmp_ParentSets[current_var]), nParentSets[current_var]) );
      }
      else
      {
         /* SCORE line */
         char** line_parts = NULL;
         int num_parts = splitAtWhitespace(this_line, &line_parts);
         current_parent_set += 1;
         Scores[current_var][current_parent_set] = atof(line_parts[0]);
         nParents[current_var][current_parent_set] = num_parts - 1;
         SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[current_var][current_parent_set]), nParents[current_var][current_parent_set]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_ParentSets[current_var][current_parent_set]), nParents[current_var][current_parent_set]) );
         for( i = 1; i < num_parts; i++ )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_ParentSets[current_var][current_parent_set][i - 1]), SCIP_MAXSTRLEN) );
            strcpy(tmp_ParentSets[current_var][current_parent_set][i - 1], line_parts[i]);
         }
         for( i = 0; i < num_parts; i++ )
            free(line_parts[i]);
         free(line_parts);
      }
   }
   fclose(file);

   /* Find the ParentSet numbers now that all names are known */
   for( i = 0; i < n; i++ )
      for( k = 0; k < nParentSets[i]; k++ )
         for( l = 0; l < nParents[i][k]; l++ )
         {
            SCIP_Bool foundMatch = FALSE;
            for( m = 0; m < n && foundMatch == FALSE; m++ )
               if( strcmp(nodeNames[m], tmp_ParentSets[i][k][l]) == 0 )
               {
                  ParentSets[i][k][l] = m;
                  foundMatch = TRUE;
               }
            if( !foundMatch )
            {
               SCIPerrorMessage("Reading failed: unable to identify node %s, a potential parent of node %s.\n", tmp_ParentSets[i][k][l], nodeNames[i]);
               return SCIP_READERROR;
            }
         }

   /* Free the allocated memory that is no longer needed */
   for( i = 0; i < n; i++ )
   {
      for( k = 0; k < nParentSets[i]; k++ )
      {
         for( l = 0; l < nParents[i][k]; l++ )
            SCIPfreeMemoryArray(scip, &(tmp_ParentSets[i][k][l]));
         SCIPfreeMemoryArray(scip, &(tmp_ParentSets[i][k]));
      }
      SCIPfreeMemoryArray(scip, &(tmp_ParentSets[i]));
   }
   SCIPfreeMemoryArray(scip, &tmp_ParentSets);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of candidate parent sets: %d\n", numCandidateParentSets);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "File reading successful\n");
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
   return SCIP_OKAY;
}
/** Reads local score BN file (in Jaakkola format)
 *
 *  Format:
 *  - first line is: number of variables
 *  then local scores for each variable
 *  - first line of section of local scores is "variable number_of_parent_sets"
 *  - other lines are like "score 3 parent1 parent2 parent3" (e.g. when there are 3 parents)
 *
 *  @param scip The SCIP instance being used.
 *  @param filename The file to read from.
 *  @param psd The parentage data that is read from the file.
 *  @param scores The score of each of the parent set combinations.
 *  @param prop The properties of the data read from the file.
 *  @return SCIP_OKAY if reading was successful, or an appropriate error code otherwise.
 */
static
SCIP_RETCODE readProblemInJaakkolaFormat(
   SCIP* scip,
   const char* filename,
   ParentSetData* psd,
   SCIP_Real*** scores,
   PropertyData* prop
)
{
   int i;
   int k;
   int l;
   int m;
   FILE *file;
   int status;
   int n;            /* number of variables */
   SCIP_Real** Scores;
   int* nParentSets;
   int** nParents;
   int*** ParentSets;
   char** nodeNames;
   char**** tmp_ParentSets;
   int numCandidateParentSets = 0;

   /* open file */
   if( strcmp(filename, "-") == 0 )
      file = stdin;
   else
      file = fopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", filename);
      return SCIP_NOFILE;
   }

   /* read number of elements */
   status = fscanf(file, "%d", &n);
   if( ! status )
   {
      SCIPerrorMessage("Reading failed: first line did not state number of variables.\n");
      return SCIP_READERROR;
   }
   assert(0 < n);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of variables: %d\n", n);
   psd->n = n;

   prop->n = psd->n;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prop->num_properties), prop->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prop->property_names), prop->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prop->property_values), prop->n) );
   for( i = 0; i < prop->n; i++ )
   {
      prop->num_properties[i] = 0;
      prop->property_names[i] = NULL;
      prop->property_values[i] = NULL;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &Scores, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParentSets, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParents, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ParentSets, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tmp_ParentSets, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeNames, n) );
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nodeNames[i]), SCIP_MAXSTRLEN) );

   *scores = Scores;
   psd->nParentSets = nParentSets;
   psd->nParents = nParents;
   psd->ParentSets = ParentSets;
   psd->nodeNames = nodeNames;

   for( i = 0; i < n; ++i )
   {
      status = fscanf(file, "%s %d", nodeNames[i], &(nParentSets[i]));
      if( ! status )
      {
         SCIPerrorMessage("Reading failed: did not get number of parents for variable %d.\n", i);
         return SCIP_READERROR;
      }

      if( nParentSets[i] == 0 )
      {
         SCIPerrorMessage("Variable %s has no parent sets. Problem is infeasible.\n", nodeNames[i]);
      }

      /*SCIPmessagePrintInfo("%d %d\n\n", i, nParentSets[i]);*/
      numCandidateParentSets += nParentSets[i];

      SCIP_CALL( SCIPallocMemoryArray(scip, &(nParents[i]), nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(Scores[i]), nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i]), nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_ParentSets[i]), nParentSets[i]) );

      for( k = 0; k < nParentSets[i]; ++k )
      {
         status = fscanf(file, "%lf %d", &(Scores[i][k]), &(nParents[i][k]));
         if( ! status )
         {
            SCIPerrorMessage("Reading failed: did not get size of parent set %d for variable %d.\n", k, i);
            return SCIP_READERROR;
         }

         if( nParents[i][k] < 0 )
         {
            SCIPerrorMessage("Reading failed: size of parent set %d for variable %d is apparently negative!\n", k, i);
            SCIPerrorMessage("Possibly a local score is being interpreted as a parent set size due to an extraneous parent on the preceding line.\n");
            return SCIP_READERROR;
         }

         SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i][k]), nParents[i][k]) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_ParentSets[i][k]), nParents[i][k]) );

         for( l = 0; l < nParents[i][k]; ++l )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(tmp_ParentSets[i][k][l]), SCIP_MAXSTRLEN) );
            status = fscanf(file, "%s", tmp_ParentSets[i][k][l]);
            if( ! status )
            {
               SCIPerrorMessage("Reading failed: did not get parent %d for parent set %d for variable %d.\n", l, k, i);
               return SCIP_READERROR;
            }
         }
      }
   }

   fclose(file);

   /* Find the ParentSet numbers now that all names are known */
   for( i = 0; i < n; i++ )
      for( k = 0; k < nParentSets[i]; k++ )
         for( l = 0; l < nParents[i][k]; l++ )
         {
            SCIP_Bool foundMatch = FALSE;
            for( m = 0; m < n && foundMatch == FALSE; m++ )
               if( strcmp(nodeNames[m], tmp_ParentSets[i][k][l]) == 0 )
               {
                  ParentSets[i][k][l] = m;
                  foundMatch = TRUE;
               }
            if( !foundMatch )
            {
               SCIPerrorMessage("Reading failed: unable to identify node %s, a potential parent of node %s.\n", tmp_ParentSets[i][k][l], nodeNames[i]);
               return SCIP_READERROR;
            }
         }

   /* Free the allocated memory that is no longer needed */
   for( i = 0; i < n; i++ )
   {
      for( k = 0; k < nParentSets[i]; k++ )
      {
         for( l = 0; l < nParents[i][k]; l++ )
            SCIPfreeMemoryArray(scip, &(tmp_ParentSets[i][k][l]));
         SCIPfreeMemoryArray(scip, &(tmp_ParentSets[i][k]));
      }
      SCIPfreeMemoryArray(scip, &(tmp_ParentSets[i]));
   }
   SCIPfreeMemoryArray(scip, &tmp_ParentSets);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of candidate parent sets: %d\n", numCandidateParentSets);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "File reading successful\n");
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
   return SCIP_OKAY;
}

/** Gets a problem name based on the file from which it was read.
 *
 *  This function copied from LOP example written by March Pfetsch.
 *
 *  @return SCIP_OKAY if the name could be found.  An error otherwise.
 */
static
SCIP_RETCODE getProblemName(
   const char* filename,         /**< input filename */
   char*       probname,         /**< output problemname */
   int         maxSize           /**< maximum size of probname */
)
{
   int i = 0;
   int j = 0;
   int l;

   /* first find end of string */
   while( filename[i] != 0 )
      ++i;
   l = i;

   /* go back until '.' or '/' or '\' appears */
   while( (i > 0) && (filename[i] != '.') && (filename[i] != '/') && (filename[i] != '\\') )
      --i;

   /* if we found '.', search for '/' or '\\' */
   if( filename[i] == '.' )
   {
      l = i;
      while( (i > 0) && (filename[i] != '/') && (filename[i] != '\\') )
         --i;
   }

   /* correct counter */
   if( (filename[i] == '/') || (filename[i] == '\\') )
      ++i;

   /* copy name */
   while( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if( j > maxSize - 1 )
         return SCIP_ERROR;
   }
   probname[j] = 0;

   return SCIP_OKAY;
}

/** Reads a problem from a file in a format other than CIP.
 *
 *  @param scip The SCIP instance in to which the problem should be read.
 *  @param filename The file to read the problem from.
 *  @param format The format the problem is given in.
 *  @return SCIP_OKAY if the problem was successfully created, or an error code otherwise.
 */
static
SCIP_RETCODE readProblemInNonCIPFormat(
   SCIP* scip,
   char* frequencyfile,
   const char* filename,
   char* format
   )
{
   int i;
   int k;
   char probname[SCIP_MAXSTRLEN];
   SCIP_Real** scores = NULL;
   ParentSetData* psd = NULL;
   PropertyData* prop = NULL;
   SCIP_Bool onlyscores;

   SCIP_PROBDATA* probdata = NULL;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/onlyscores", &onlyscores) );


   SCIP_CALL( getProblemName(filename, probname, SCIP_MAXSTRLEN) );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "File name:\t\t%s\n", filename);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Problem name:\t\t%s\n", probname);
   if( !onlyscores )
   {
      SCIP_CALL( SCIPallocMemory(scip, &probdata) );
      SCIP_CALL( SCIPcreateProb(scip, probname, probdelorigBN, NULL, NULL, NULL, NULL, NULL, probdata) );
   }

   SCIP_CALL( SCIPallocMemory(scip, &psd) );
   SCIP_CALL( SCIPallocMemory(scip, &prop) );
   PR_initialise(prop);
   psd->PaVars = NULL;
   if( strcmp(format, "dat") == 0 )
   {
      int num_delims;
      char* delims;
      SCIP_Bool merge_delims;
      SCIP_CALL( getDelimiters(scip, &num_delims, &delims, &merge_delims) );
      SCIP_CALL( SC_readProblemInDataFormat(scip, filename, num_delims, delims, merge_delims, psd, &scores, prop) );
      SCIPfreeMemoryArray(scip, &delims);
      if( onlyscores )
      {
         SCIP_CALL( IO_printScoresInJKLFormat(scip, psd, scores) ); 
            for( i = 0; i < psd->n; i++ )
            {
               for( k = 0; k < psd->nParentSets[i]; k++ )
                  SCIPfreeMemoryArray(scip, &(psd->ParentSets[i][k]));
               SCIPfreeMemoryArray(scip, &(scores[i]));
               SCIPfreeMemoryArray(scip, &(psd->nParents[i]));
               SCIPfreeMemoryArray(scip, &(psd->ParentSets[i]));
               SCIPfreeMemoryArray(scip, &(psd->nodeNames[i]));
            }
         SCIPfreeMemoryArray(scip, &scores);
         SCIPfreeMemoryArray(scip, &(psd->nParents));
         SCIPfreeMemoryArray(scip, &(psd->nParentSets));
         SCIPfreeMemoryArray(scip, &(psd->ParentSets));
         SCIPfreeMemoryArray(scip, &(psd->nodeNames));
         SCIPfreeMemoryArray(scip, &psd);
         SCIP_CALL( PR_deallocatePropertyData(scip, &prop) );
         return SCIP_OKAY;
      }
   }
   else if( strcmp(format, "jkl") == 0 )
   {
      SCIP_CALL( readProblemInJaakkolaFormat(scip, filename, psd, &scores, prop) );
   }
   else if( strcmp(format, "pss") == 0 )
   {
      SCIP_CALL( readProblemInPSSFormat(scip, filename, psd, &scores, prop) );
   }
   else if( strcmp(format, "gen") == 0 )
   {
      int num_delims;
      char* delims;
      SCIP_Bool merge_delims;
      SCIP_CALL( getDelimiters(scip, &num_delims, &delims, &merge_delims) );
      SCIP_CALL( PS_readProblemInGenomeFormat(scip, filename, num_delims, delims, merge_delims, frequencyfile, psd, &scores, prop) );
      SCIPfreeMemoryArray(scip, &delims);
   }
   else
   {
      /* This should be picked up before now, but do something sensible if not. */
      SCIPerrorMessage("Trying to read in data in an unknown format - %s\n", format);
      return SCIP_ERROR;
   }

   if( !onlyscores )
   {
      probdata->psd = psd;
      probdata->scores = scores;
      probdata->ancestorvars = NULL;
      probdata->im_vars = NULL;
      probdata->posvars = NULL;
      probdata->totalordervars = NULL;
      probdata->ispa = NULL;
      probdata->isFemale = NULL;
      probdata->pasize = NULL;
      probdata->posind = NULL;
      probdata->imsetvars = NULL;
      probdata->nch = NULL;
   }

   SCIP_CALL( createBasicModel(scip, psd, scores) );
   SCIP_CALL( MD_setParentSetData(scip, psd) );
   SCIP_CALL( MD_setPropertyData(scip, prop) );

   /* for( i = 0; i < psd->n; i++ ) */
   /*    SCIPfreeMemoryArray(scip, &(scores[i])); */
   /* SCIPfreeMemoryArray(scip, &scores); */
   /* SCIP_CALL( PS_deallocateParentSetData(scip, &psd, FALSE) ); */

   SCIP_CALL( PR_deallocatePropertyData(scip, &prop) );

   SCIPdebugMessage("Parent set data structure set up.\n");

   return SCIP_OKAY;
}
/** Reads a problem from a file in cip format.
 *
 *  @param scip The SCIP instance in to which the problem should be read.
 *  @param filename The file to read the problem from.
 *  @return SCIP_OKAY if the problem was successfully created, or an error code otherwise.
 */
static
SCIP_RETCODE readProblemInCIPFormat(
   SCIP* scip,
   const char* filename
   )
{
   SCIP_CALL( SCIPreadProb(scip, filename, "cip") );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "File name:\t\t%s\n", filename);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Problem name:\t\t%s\n", SCIPgetProbName(scip));
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
   return SCIP_OKAY;
}

/** Reads in a problem from a file.
 *
 *  Currently supported file formats are "jak", "cip" and "dat".  The one to use
 *  is determined using the <code>-f</code> cmmand line argument and the
 *  filename extension.  See getInputFormat() for details.
 *
 *  @param scip The SCIP instance to read the problem in to.
 *  @param filename The filename from which to read the problem.
 *  @return SCIP_OKAY if the problem was successfully created, or an error code otherwise.
 */
SCIP_RETCODE BN_readProblem(
   SCIP* scip,
   char* inputformat,
   char* frequencyfile,
   const char* filename
   )
{
   ParentSetData* psd;
   char* format = getInputFormat(scip, inputformat, filename);
   SCIP_Bool onlyscores;

   int i;
   int k;

   SCIP_Bool allowmultiaggregation;
   
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/onlyscores", &onlyscores) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/allowmultiaggregation", &allowmultiaggregation) );

   if( onlyscores && strcmp(format, "dat") != 0 )
   {
      SCIPerrorMessage("Can only use onlyscores parameter when input is data.\n");
      return SCIP_ERROR;
   }      

   if( strcmp(format, "jkl") == 0 || strcmp(format, "dat") == 0 || strcmp(format, "gen") == 0 || strcmp(format, "pss") == 0 )
      SCIP_CALL( readProblemInNonCIPFormat(scip, frequencyfile, filename, format) );
   else if( strcmp(format, "cip") == 0 )
   {
      SCIP_CALL( readProblemInCIPFormat(scip, filename) );
      if( MD_getParentSetData(scip) == NULL )
      {
         SCIPwarningMessage(scip, "No metadata found in input file.\n");
         return SCIP_OKAY;
      }
   }
   else
   {
      /* Shouldn't ever get here if getInputFormat() is working correctly */
      SCIPerrorMessage("This file format is unsupported.\n");
      return SCIP_ERROR;
   }
   SCIPfreeMemoryArray(scip, &format);

   if( onlyscores )
      return SCIP_OKAY;

   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM);

   psd = MD_getParentSetData(scip);

   /* sort the parent sets 
    useful for the distinguished representative constraint, for example 
   */
   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
         SCIPsortInt(psd->ParentSets[i][k],psd->nParents[i][k]);

   SCIP_CALL( addAdditionalConstraints(scip, psd) );
   SCIP_CALL( setBranchingPriorities(scip, psd) );
   if( !allowmultiaggregation )
      SCIP_CALL( setSpecialProperties(scip, psd) );

   SCIPdebugMessage("Additional constraints, branching priorities and special properties all added.\n");	

   return SCIP_OKAY;
}



/* Functions related to finding the n best networks */
/** Gets the number of most likely Bayesian networks that should be found.
 *
 *  @param scip The SCIP inatance used for finding the networks.
 *  @return The number of Bayesian networks that should be found.
 */
int BN_getNumberOfRepeats(
   SCIP* scip
   )
{
   int nbns;
   SCIPgetIntParam(scip, "gobnilp/nbns", &nbns);
   return nbns;
}

/** Adds a constraint that prevents networks in the same MEC as the current best network being found again.
 *
 *  @param scip The SCIP instance being used for the optimisation.
 *  @param run The iteration of the loop that this soilution was found on.
 *
 *  @return SCIP_OKAY if a constraint could be added, or an error otherwise.
 */
static
SCIP_RETCODE BN_addMECNonRepetitionConstraint(
   SCIP* scip,
   int run
   )
{
   int i;
   int j;
   int k;

   char consname[SCIP_MAXSTRLEN];
   int num_included = 0;
   SCIP_VAR** included_vars;
   SCIP_CONS* cons;

   ParentSetData* psd = MD_getParentSetData(scip);
   SCIP_SOL* sol = SCIPgetBestSol(scip);

   SCIP_VAR* edge_i_j;

   SCIP_PROBDATA* probdata;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   
   /* Allocate the memory */
   SCIP_CALL( SCIPallocMemoryArray(scip, &included_vars, (psd->n) * (psd->n) + (psd->n) * (psd->n) * (psd->n)) );

   /* Go through the edge vars */
   for( i = 0; i < psd->n - 1; i++ )
      for( j = i + 1; j < psd->n; j++ )
      {
         edge_i_j = get_edge(psd,i,j);
         if( edge_i_j != NULL )
         {
            if( SCIPgetSolVal(scip, sol, edge_i_j) > 0.5 )
            {
               SCIP_VAR* neg_var;
               SCIP_CALL( SCIPgetNegatedVar(scip, edge_i_j, &neg_var) );
               included_vars[num_included++] = neg_var;
            }
            else
            {
               included_vars[num_included++] = edge_i_j;
            }
         }
      }

   /* Go through the immorality vars */
   for( i = 0; i < psd->n - 1; i++ )
      for( j = i + 1; j < psd->n; j++ )
         for( k = 0; k < psd->n; k++ )
            if( i != k && j != k )
            {
               /* imvars array is global */
               if( SCIPgetSolVal(scip, sol, probdata->im_vars[i][j][k]) > 0.5 )
               {
                  SCIP_VAR* neg_var;
                  SCIP_CALL( SCIPgetNegatedVar(scip, probdata->im_vars[i][j][k], &neg_var) );
                  included_vars[num_included++] = neg_var;
               }
               else
               {
                  included_vars[num_included++] = probdata->im_vars[i][j][k];
               }
            }

   /* Create and add the new constraint */
   SCIP_CALL( SCIPfreeTransform(scip) );

   SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ruleoutMEC#%d", run);

   SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &cons, consname, num_included, included_vars) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* Free the memory */
   SCIPfreeMemoryArray(scip, &included_vars);

   return SCIP_OKAY;
}

/** Adds a constraint that prevents the current best network being found again.
 *
 *  @param scip The SCIP instance being used for the optimisation.
 *  @param run The iteration of the loop that this soilution was found on.
 *
 *  @return SCIP_OKAY if a constraint could be added, or an error otherwise.
 */
SCIP_RETCODE BN_addNonRepetitionConstraint(
   SCIP* scip,
   int run
   )
{
   int n_empty = 0;
   int i;
   int k;
   ParentSetData* psd = MD_getParentSetData(scip);
   SCIP_SOL* sol = SCIPgetBestSol(scip);
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS *cons;
   int* chosen;
   SCIP_Bool kbestMEC;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/kbestMEC", &kbestMEC) );
   if( kbestMEC )
   {
      SCIP_CALL( BN_addMECNonRepetitionConstraint(scip, run) );
      return SCIP_OKAY;
   }
   else
   {

      SCIP_CALL( SCIPallocMemoryArray(scip, &chosen, psd->n) );

      /* record which BN just found before doing 'free transform' */

      for( i = 0; i < psd->n; ++i )
      {
         SCIP_Bool no_parents = TRUE;
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            SCIP_Real val = SCIPgetSolVal(scip, sol, psd->PaVars[i][k]);
            assert(SCIPisIntegral(scip, val));
            if( val > 0.5 )
            {
               chosen[i] = k;
               no_parents = FALSE;
               break;
            }
         }
         if( no_parents )
         {
            n_empty++;
            chosen[i] = -1;
         }
      }

      SCIP_CALL( SCIPfreeTransform(scip) );
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ruleout#%d", run);
      /* maybe change this to set covering constraint */
      /* rather than rely on upgrading */
      /* basically the same as CUTOFF_CONSTRAINT(addBinaryCons) in cons_countsols.c */
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(scip), (psd->n) - 1 - n_empty,
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

      for( i = 0; i < psd->n; ++i )
      {
         if( chosen[i] == -1 )
            for( k = 0; k < psd->nParentSets[i] - 1; ++k )
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][k], -1) );
         else
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][chosen[i]], 1) );
      }
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      SCIPfreeMemoryArray(scip, &chosen);

      return SCIP_OKAY;
   }
}

/* Functions for printing */
/** Prints the current problem.
 *
 *  @param scip The SCIP instance for which the problem should be printed.
 *  @param run The iteration of the main loop that the problem is to be solved on.
 *  @return SCIP_OKAY if printing succeeded or an appropriate error code otherwise.
 */
SCIP_RETCODE BN_printProblem(
   SCIP* scip,
   int run
   )
{
   SCIP_CALL( IO_printProblem(scip, run) );
   return SCIP_OKAY;
}

/** Prints appropriate information about each optimal solution obtained.
 *
 *  @param scip The SCIP instance for which the solution has been found.
 *  @param run The iteration of the main loop that the solution was found on.
 *  @return SCIP_OKAY if printing succeeded or an appropriate error code otherwise.
 */
SCIP_RETCODE BN_doIterativePrint(
   SCIP* scip,
   MA_info* ma_info,
   int run
   )
{
   SCIP_CALL( IO_doIterativePrint(scip, ma_info, MD_getParentSetData(scip), run) );
   return SCIP_OKAY;
}

/** Prints any of the current SCIP or GOBNILP parameters not at their default value.
 *
 *  @param scip The SCIP instance to consult the parameters of.
 *  @return SCIP_OKAY if the parameters were printed correctly, or an error code otherwise.
 */
SCIP_RETCODE BN_printParameters(
   SCIP* scip
   )
{
   SCIP_CALL( IO_printParameters(scip) );
   return SCIP_OKAY;
}

/** Prints a header which describes the GOBNILP and SCIP systems being used.
 *
 *  @param scip The SCIP instance that is being used.
 *  @return SCIP_OKAY if printing succeeded or an error code otherwise.
 */
SCIP_RETCODE BN_printHeader(
   SCIP* scip
   )
{
   SCIP_CALL( IO_printHeader(scip) );
   return SCIP_OKAY;
}

/** Prints the scores for the family variables in Jaakkola format.
 *
 *  @param scip The SCIP instance to print scores for.
 *  @return SCIP_OKAY if printing succeeded or an error otherwise.
 */
SCIP_RETCODE BN_printScores(
   SCIP* scip
   )
{
   SCIP_CALL( IO_printScoresInJKLFormat(scip, MD_getParentSetData(scip), NULL) );
   SCIP_CALL( IO_printScoresInPSSFormat(scip, MD_getParentSetData(scip), MD_getPropertyData(scip)) );
   return SCIP_OKAY;
}

/** Prints solutions found using countsols
 *
 *  @param scip The SCIP instance to print scores for.
 *  @return SCIP_OKAY if printing succeeded or an error otherwise.
 */
SCIP_RETCODE BN_printcountsols(
   SCIP* scip,
   char* filename
   )
{
   SCIP_CALL( IO_printcountsols(scip, MD_getParentSetData(scip), filename) );
   return SCIP_OKAY;
}
