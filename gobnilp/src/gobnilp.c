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
 *  Main entry point for the program.
 *
 *  This file contains only general information about the sequence that operations
 *  should be performed in.  Much of what is in this file should be applicable to
 *  ILP problems in general, not just Bayesian network learning.  Functions that are
 *  more specific to BN learning are found in @link probdata_bn.c @endlink .
 */

/** @mainpage notitle
 *
 * @section whatisthis Should I be reading this documentation?
 * This is the documentation for the @b GOBNILP source code for development purposes.\n
 *
 * For help in using @b GOBNILP, see the PDF manual included in the distribution.\n
 *
 *
 * @section overview Overview of the code layout
 *
 * @subsection bge_matrix bge_matrix.c
 * @c bge_matrix.c contains matrix routines for scoring a Gaussian network
 *
 * @subsection bge_posterior bge_posterior.c
 * @c bge_posterior.c contains functions for computing the posterior matrix for Gaussian networks
 *
 * @subsection bge_score bge_score.c
 * @c bge_score.c contains the main functions for computing BGe scores
 *
 * @subsection bge_vector bge_vector.c
 * @c bge_vector.c has vector routines for Gaussian networks
 *
 * @subsection overview_circuit_cuts circuit_cuts.c
 * @c circuit_cuts.c contains code to look for cluster cuts by searching for cyclic graphs in a given solution (typically the LP relaxation).
 * This code is 'logically' part of @c cons_dagcluster.c since it adds constraints (as cutting planes) which follow from acyclicity. 
 *
 * @subsection overview_cons_chordal cons_chordal.c
 * @c cons_chordal.c implements a constraint handler for constraints which ensure that DAGs have no immoralities and thus are equivalent to
 * chordal undirected graphs (aka decomposable models).
 *
 * @subsection overview_cons_ci cons_ci.c
 * @c cons_ci.c implements a constraint handler for conditional independence constraints. 
 *
 * @subsection overview_cons_dagcluster cons_dagcluster.c
 * Enforcing a constraint to rule out cycles in the Bayesian network is non-trivial.  The
 * approach in @b GOBNILP is to use &lsquo;cluster constraints&rsquo;.  As there are
 * exponentially many of these, they are added during solving as cutting planes.  
 * The main code associated with creating and updating these constraints is contained in
 * @c cons_dagcluster.c .\n
 *
 * @subsection overview_cons_distinguishedrep cons_distinguishedrep.c
 * @c cons_distinguishedrep.c implements a constraint handler for constraints which ensure that only one representative from each
 * equivalence class of Markov equivalent DAGs is feasible.
 *
 * @subsection overview_cons_lop cons_lop.c
 * @c cons_lop.c implements a constraint handler for linear ordering constraints. This constraint handler was written by Marc Pfetsch and 
 * is part of the LOP 
 * example that comes with the SCIP distribution. It is copied across when @b GOBNILP is installed.
 *
 * @subsection overview_cons_partialordering cons_partialordering.c 
 * @c cons_partialordering.c implements a constraint
 * handler for partial order constraints. This constraint handler implements constraints ensuring that DAGs are
 * consistent with a partial order. The partial order can be minimal in which case it is the ancestor relation for the
 * DAG. The constraint handler was written by small modifications to @c cons_lop.c . This file is created when @b
 * GOBNILP is installed by applying a patch to a copy of @c cons_lop.c.
 *
 * @subsection overview_cons_vanilla cons_vanilla.c
 * @c cons_vanilla.c implements a constraint handler for constraints stating a DAG that no 'obviously' suboptimal DAGs are feasible under
 * the assumption that no constraints (other than acyclicity) are imposed on DAGs.
 *
 * @subsection overview_convex4_cuts convex4_cuts.c
 * @c convex4_cuts.c contains code to look for cuts which are (generalisations of) facets of the convex hull of 4-node DAGs. 
 * This code is 'logically' part of @c cons_dagcluster.c since it adds constraints (typically as cutting planes) which follow from acyclicity. 
 *
 * @subsection overview_disp_clearcols disp_clearcols.c
 * @c disp_clearcols.c contains code to generate clearer column headings for @b GOBNILP output
 *
 * @subsection overview_event_splitdag event_splitdag.c
 * @c event_splitdag.c provides a simple event handler to detect when the program should 
 * attempt to split dagcluster constraints into their individual strongly connected components
 *
 * @subsection overview_fractional_circuit_cuts fractional_circuit_cuts.c
 * @c fractional_circuit_cuts.c generalises the approach taken in @c circuit_cuts.c
 * This code is 'logically' part of @c cons_dagcluster.c since it adds constraints (as cutting planes) which follow from acyclicity. 
 *
 * @subsection overview_gobnilp gobnilp.c
 * The main entry point in to the code is in the file @c gobnilp.c .  This
 * file contains a single procedure mostly laying out the general order of operations for
 * a problem in SCIP.  The idea is that this main calling function should able to be
 * largely ignorant of how the underlying implementation of the problem so that
 * representations and functionality can be altered without affecting the general order
 * of operations of the program.
 * Unless you are adding a new plugin, command line parameter or altering the basic
 * workflow of the code, you probably shouldn't be making alterations in this file.
 *
 * @subsection overview_heur_sinks heur_sinks.c
 * @c heur_sinks.c has most of the code necessary to implement a heuristic for finding
 * Bayesian networks.  In general, this file should be self contained, except for the
 * function needed to include the heuristic in the program.  The one excpetion to this is
 * that the heuristic shouldn't assign variables for the pedigree problem.  These should
 * instead be set from @c pedigrees.c in order to maintain the correct modularisation of
 * data.
 *
 * @subsection overview_matroid_cuts matroid_cuts.c
 * @c matroid_cuts.c contains functions for finding 'matroid' cuts (as proposed by Milan Studeny)
 *
 * @subsection overview_metadata metadata.c
 * @c metadata.c contains metadata relating to the problem that is not necessarily part of the ILP that will be solved.
 *
 * @subsection overview_model_averaging model_averaging.c
 * @c model_averaging.c contains all the functionality to perform model averaging over
 * the n best networks found.  As far as possible, the fact that @b GOBNILP is capable
 * of model averaging should be hidden from other parts of the program.
 * There are two exceptions to this.  First, @c gobnilp.c needs to know about model
 * averaging in order that it can make appropriate calls to perform the model averaging.
 * Second, @c output.c needs to be aware that the program can perform model averaging so
 * that the output functions can produce suitable output if a model average is to be
 * printed.
 *
 * @subsection overview_output output.c
 * As the name suggests, @c output.c contains functions related to outputing data.  In
 * fact, with a small number of exceptions, all functions related to file output should
 * appear in this file.  The exceptions to this are those functions which would break
 * the modular nature of the program if put into this file.  At present, the only
 * substantial output function not appearing in this file is that for outputting
 * pedigrees, which instead appears in @c pedigrees.c .
 * Functions appearing in this file for outputting solutions should also be capable of
 * sensible output when called with model averaging data.
 *
 * @subsection overview_data_structures parent_set_data.c
 * @c parent_set_data.c contains code for copying, deleting and altering @c ParentSetData which contains information
 * on candidate parent sets for vertices.
 * The header file @c parent_set_data.h
 * defines several
 * data structures that are used throughout the program.  Most importantly, this is
 * where the main data structure for storing instance data @c ParentSetData
 * is defined.
 *
 * @subsection overview_pedigree_data pedigree_data.c
 * @c pedigree_data.c contains functions for managing pedigree data (stored in a @c PedigreeData instance)
 *
 * @subsection overview_pedigree_scorer pedigree_scorer.c
 * @c pedigree_scorer.c contains functions for creating scores for pedigrees.
 *
 * @subsection overview_pedigrees pedigrees.c
 * @b GOBNILP is capable of learning pedigrees in addition to general purpose Bayesian
 * network learning.  As much of the additional behaviour that this requires is unrelated
 * to normal Bayesian network learning, the code for it is all assembled in one place
 * where it can be altered without affecting the main program.  Any functions related to
 * pedigrees should appear in this file.  Outside @c pedigrees.c , the data and constraints
 * that the pedigree part of the program is using should be treated as a black box.
 * A function @link PD_inPedigreeMode @endlink is defined which allows the rest of the
 * program to determine whether pedigrees are being used and call the appropriate
 * pedigree specific functions if so.
 *
 *
 * @subsection overview_probdata_bn probdata_bn.c
 * The majority of operations that are more specific to Bayesian network learning are in
 * the file @c probdata_bn.c .  For example, the functions that create the
 * variables and constraints needed to learn BNs are found in this file.
 * If adding a relatively small change in functionality to @b GOBNILP , this is probably
 * the first place to consider adding it.
 *
 * @subsection overview_property_data property_data.c
 * @c property_data.c contains functions related to managing @c PropertyData.
 *
 * @subsection overview_scoring scoring.c
 * @c scoring.c computes local scores (currently only BDeu scores) for candidate parent sets from complete discrete data. 
 * Includes code which allows certain constraints (e.g. conditional independence ones) to limit which
 * local scores need be computed. ADtrees (plus some additional caching of log-likelihoods) are used to improve
 * efficiency over a naive approach.
 *
 * @subsection overview_set_packing_cuts set_packing_cuts.c
 * @c set_packing_cuts.c adds initial constraints (not cuts despite the name of the file!) of the sort \f$I(A \leftarrow B,C) +   I(B \leftarrow A,C) +   I(C \leftarrow A,B) \leq 1\f$. These are 'redundant' but improve efficiency. See inequality (12) in Bartlett and Cussens (UAI'13).
 *
 * @subsection overview_solution_info solution_info.c
 * @c solution_info contains code for collecting and storing information on a solution (typically an LP solution)
 * which is then used by eg. several separators
 *
 * @subsection overview_stack stack.c
 * @c stack.c implements a data type representing a stack of integers
 *
 * @subsection overview_subip_cuts subip_cuts.c
 * @c subip_cuts.c constructs and solves a sub-IP for finding cluster cuts. Used by both @c cons_dagcluster.c and @c cons_ci.c
 *
 * @subsection overview_utils utils.c
 * @c utils.c contains a number of functions that may be useful in several other files.
 * At present most of these are shortcuts to several commonly used SCIP functions which have lots
 * of parameters but which are nearly always called with the same parameters for many of
 * these.  The functions in this file just make the rest of the code easier to read by
 * hiding this distraction.
 *
 * @subsection overview_vector vector.c
 * @c vector.c implements a data type representing a resizeable array of integers
 *
 * @subsection overview_vectorlist vectorlist.c
 * @c vectorlist.c implements a data type representing a resizeable array resizeable arrays of integers
 *
 * @subsection overview_versiongit versiongit.h
 * @c versiongit.h just contains some versioning information and probably shouldn't be
 * changed, except to update the version numbers.
 *
 * @subsection overview_summary Summary
 * @c gobnilp.c has the overall structure of the program.  @c probdata_bn.c provides the
 * functions related to creating the Bayesian network problem.  @c parent_set_data.h
 * contain the main shared data structures for the programs.  Pedigree functions are
 * all in @c pedigrees.c , the behaviour of which should be separated from the rest of
 * the program, such that pedigree related constraints can be added or removed without
 * affecting the correctness of the rest of the program.  @c cons_dagcluster.c and @c
 * heur_sinks.c provide SCIP plugins and should be entirely self contained except for
 * calls to pedigree functions.  The program should be unaware of model averaging except
 * for calls from @c gobnilp.c and @c output.c .
 *
 */

#include <string.h>
#include <stdio.h>
#include <scip/cons_countsols.h>
#include "scip/set.h" 
#include "scip/paramset.h"
#include "probdata_bn.h"
#include "utils.h"

#include "scip/struct_set.h"
#include "scip/struct_scip.h" 

/** The default file from which to attempt to read parameters. */
#define DEFAULT_GOBNILP_PARAMS_FILE "gobnilp.set"

#define MAXFILENAMELENGTH 100
#define MAXTMPSTRINGLENGTH 100

/** Information gathered by parsing command line arguments */
typedef struct
{
   char* parameterfile;      /**< The default file format in which the problem is given. */
   char* inputformat;        /**< The input format that the problem file is to be given in. */
   SCIP_Bool exitbeforefirstsolve; /**< Whether the program should just read and generate the problem and not solve it. */
   char* frequencyfile;      /**< The file containing the allele frequencies for pedigree scoring. */
   char* inputfile;      /**< The file containing the data */
   int verblevel;            /**< Verbosity level */
} CLA_info;   

/** Parses the command line arguments.
 *
 *  @param argc The number of command line arguments.
 *  @param argv The command line arguments.
 *  @param keys The keys extracted from the command line arguments.
 *  @param values The values extracted from the command line arguments.
 */
static
void parseArguments(
   int    argc,
   char** argv,
   char*  keys,
   char** values
   )
{
   int i;
   for( i = 1; i < argc - 1; i++ )
      if( argv[i][0] != '-' )
      {
         fprintf(stderr, "ERROR: Each optional argument must be preceded by '-': %s\n", argv[i]);
         exit(-1);
      }
      else
      {
         keys[i - 1] = argv[i][1];
         if( argv[i][2] == '=' )
            values[i - 1] = argv[i] + 3;
         else
            values[i - 1] = argv[i] + 2;
      }
}

/** Sets default values for values which can bet set on the command line */
static
void set_cla(
   CLA_info* cla_info  /**< Command line argument information to be initialised. */
   )
{
  /* initialise values settable via command line to default values */
   cla_info->parameterfile = (char* ) DEFAULT_GOBNILP_PARAMS_FILE;
   cla_info->inputformat = (char* ) "";
   cla_info->exitbeforefirstsolve = FALSE;
   cla_info->frequencyfile = (char* ) "";
   cla_info->verblevel = 3;
}

static
SCIP_Bool bad_status(
   int status
   )
{
   if( status != 1 )
   {
      fprintf(stderr, "Failed to find 'read <<inputfilename>> on standard input'\n");
      return TRUE;
   }
   else
      return FALSE;
}

/** Reads the command line arguments and stores information .
 *
 *  @param cla_info Command line argument information to be updated.
 *  @param argc The number of command line arguments.
 *  @param argv The command line arguments.
 *  @return SCIP_OKAY if reading succeeded or an error otherwise.
 */
static
void readCommandLineArgs(
   CLA_info* cla_info,
   int    argc,
   char** argv
   )
{
   int i;
   char* keys;
   char** values;
   char tmpstr[MAXTMPSTRINGLENGTH];
   char input_filename[MAXFILENAMELENGTH];
   char settings_filename[MAXFILENAMELENGTH];
   int status;
   
   if( argc == 1 )
   {
      fprintf(stderr, "No arguments given, assuming doing 'make test'\n");
      fprintf(stderr, "Looking for 'load <<settings_filename>> on standard input'\n");
      status = scanf("%s", tmpstr);
      if( bad_status(status) )
         exit(-1);
      while( strcmp(tmpstr,"load") != 0 )
      {
         status = scanf("%s", tmpstr);
         if( bad_status(status) )
            exit(-1);
      }
      status = scanf("%s", settings_filename);
      if( bad_status(status) )
         exit(-1);
      cla_info->parameterfile = settings_filename;
      
      fprintf(stderr, "Looking for 'read <<input_filename>> on standard input'\n");
      status = scanf("%s", tmpstr);
      if( bad_status(status) )
         exit(-1);
      while( strcmp(tmpstr,"read") != 0 )
      {
         status = scanf("%s", tmpstr);
         if( bad_status(status) )
            exit(-1);
      }
      status = scanf("%s", input_filename);
      if( bad_status(status) )
         exit(-1);
      cla_info->inputfile = input_filename;

      /* must have maximal verbosity for 'make test' */
      cla_info->verblevel = 5;
      
      return;
   }

   keys = (char *) malloc( (argc-2) * sizeof(char));
   values = (char **) malloc( (argc-2) * sizeof(char *)); 
   parseArguments(argc, argv, keys, values);
   for( i = 0; i < argc - 2; i++ )
   {
      switch(keys[i])
      {
      case 'g':
         /* no need to copy string ('inside' argv) */
         cla_info->parameterfile = values[i];
         break;
      case 'f':
         cla_info->inputformat = values[i];
         break;
      case 'q':
         cla_info->frequencyfile = values[i];
         break;
      case 'x':
         cla_info->exitbeforefirstsolve = TRUE;
         break;
      case 'v':
         cla_info->verblevel = atoi(values[i]);
         break;
      default:
         fprintf(stderr, "Warning: Unrecognised optional argument: %s\n", argv[i + 1]);
         break;
      }
   }
   cla_info->inputfile = argv[argc - 1];
   free(keys);
   free(values);
}


/** Main function which starts the solution of the BN learning problem.
 *  @param argc The number of command line arguments supplied to the program.
 *  @param argv The command line arguments supplied to the program.
 */
int main(int argc, char** argv)
{
   SCIP* scip = NULL;
   const char* filename;
   int i;
   int nbns;
   SCIP_Bool foundAll = FALSE;
   SCIP_Bool onlyscores;
   SCIP_Real objlimit;
   SCIP_Bool countsols;
   SCIP_Bool collect;
   SCIP_Longint sollimit; 
   SCIP_Bool valid;
   char* tmp_countsols_filename;
   char countsols_filename[SCIP_MAXSTRLEN];
   SCIP_Bool gobnilpparams;
   MA_info ma_info;
   CLA_info cla;

   
#if SCIP_VERSION < 400
   printf("This version of GOBNILP requires SCIP 4.0.0 or above.\n");
   exit(1);
#endif

   /* initialise values settable via command line to default values */
   set_cla(&cla);
   
   /* process command line arguments */
   readCommandLineArgs(&cla, argc, argv);
   
   /* Create SCIP instance */
   SCIP_CALL( SCIPcreate(&scip) );

   /* Set verbosity level now to control how header is printed */
   SCIPsetIntParam(scip, "display/verblevel", cla.verblevel);

   /* Print header describing SCIP and GOBNILP versions */
   SCIP_CALL( BN_printHeader(scip) );

   /* Include GOBNILP plugins, e.g. acyclity constraint handler */
   SCIP_CALL( BN_includePlugins(scip) );

   /* Add (most) GOBNILP parameters to SCIP instance */
   SCIP_CALL( BN_addParameters(scip) );

   /* Add GOBNILP parameters controlling model averaging */
   SCIP_CALL( MA_addAveragingParameters(scip, &ma_info) );

   /* Find out whether we should set SCIP parameters to values optimised for GOBNILP */
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/gobnilpparams", &gobnilpparams) );

   /* Conditionally set SCIP parameters to values optimised for GOBNILP */
   if( gobnilpparams )
      SCIP_CALL( BN_setParamaterDefaults(scip) );

   /* Conditionally read in user set parameter values */
   filename = cla.parameterfile;
   if( SCIPfileExists(filename) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Reading parameter file <%s>.\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
   }
   else
   {
      SCIPwarningMessage(scip, "Parameter file <%s> not found - using default settings.\n", filename);
   }

   /* Set which output columns should be displayed */
   SCIP_CALL( BN_suppresscols(scip) );

   /* Conditionally print out parameters not at their default values */
   SCIP_CALL( BN_printParameters(scip) );

   /* Read problem data from file specified as last command line argument */
   /* If onlyscores is set, no problem is created and any computed scores are printed out */
   assert(SCIPgetStage(scip) == SCIP_STAGE_INIT); 	
   SCIP_CALL( BN_readProblem(scip, cla.inputformat, cla.frequencyfile, cla.inputfile) );

   /* if only scores (computed from data) are required then no SCIP problem is created and
    * there is an early exit */
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/onlyscores", &onlyscores) );
   assert(onlyscores || SCIPgetStage(scip) == SCIP_STAGE_PROBLEM); 	
   if( onlyscores )
   {
      SCIP_CALL( SCIPfree(&scip) );
      BMScheckEmptyMemory();
      return 0;
   }

   /* Print out scores if desired. (Default is not to.) */
   SCIP_CALL( BN_printScores(scip) );

   /* Initialise values for model averaging */
   SCIP_CALL( MA_createAverageDataStructure(scip, &ma_info) );

   /* Obtain various parameter values associated with multiple solutions */
   SCIP_CALL( SCIPgetRealParam(scip, "gobnilp/objlimit", &objlimit) );
   if( !SCIPisInfinity(scip, -objlimit) )
      SCIP_CALL( SCIPsetObjlimit(scip,objlimit) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/countsols", &countsols) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/countsols/collect", &collect) );
   SCIP_CALL( SCIPgetLongintParam(scip,"gobnilp/countsols/sollimit", &sollimit) );
   /* make a copy here, since otherwise SCIPsetParamsCountsols somehow messes it up */
   SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/outputfile/countsols", &tmp_countsols_filename) );
   strcpy(countsols_filename,tmp_countsols_filename);

   /* Solve the model (at most n times) */
   nbns = BN_getNumberOfRepeats(scip);
   for( i = 0; i < nbns && !foundAll; i++ )
   {
      SCIP_CALL( BN_printProblem(scip, i) );
      if( i == 0 && cla.exitbeforefirstsolve )
         break;
      SCIPdebugMessage("Start of solving.\n");

      if( countsols )
      {
         SCIP_CALL( SCIPsetParamsCountsols(scip) );
         /* put separating back on! */
         SCIP_CALL( SCIPparamsetSetSeparating(scip->set->paramset, scip->set, scip->messagehdlr, SCIP_PARAMSETTING_DEFAULT, TRUE) ); 
         SCIP_CALL( SCIPsetBoolParam(scip, "constraints/countsols/collect", collect) );
         SCIP_CALL( SCIPsetLongintParam(scip, "constraints/countsols/sollimit", sollimit) );
         SCIP_CALL( SCIPcount(scip) );

         printf("Found this many solutions: %lld\n", SCIPgetNCountedSols(scip, &valid) );
         if( !valid )
            printf("Number of solutions found too great to report correctly.\n");

         if( collect )
         {
            SCIP_CALL( BN_printcountsols(scip, countsols_filename) ); 
            printf("Solutions written to %s.\n", countsols_filename);
         }
      }
      else
      {
         /* Normal solving: only one BN required */
         SCIP_CALL( SCIPsolve(scip) );
      }

      if( !countsols )
      {
         if( SCIPgetBestSol(scip) == NULL )
         {
            if( i == 0 )
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No solutions possible.\n");
            else
               SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No further solutions possible.\n");
            foundAll = TRUE;
         }
         else
         {
            SCIP_CALL( MA_updateAverageDataStructure(scip, &ma_info, SCIPgetBestSol(scip)) );
            SCIP_CALL( BN_doIterativePrint(scip, &ma_info, i) );
            SCIP_CALL( SCIPfreeTransform(scip) );
            /* SCIP_CALL( BN_addNonRepetitionConstraint(scip, i) ); */
         }
      }
   }

   /* Deallocate any memory being used */
   SCIP_CALL( MA_destroyAverageDataStructure(scip, &ma_info) );
   SCIP_CALL( SCIPfreeProb(scip) );
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   return 0;
}
