/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2016 James Cussens, Mark Bartlett        *
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
 *  Defines all the functionality needed to find a heuristic solution to the Bayesian network problem.
 *
 *  The heurstic is based on repeatedly choosing the best parent set for a node to add to the network in
 *  a greedy manner, subject to the constraint that adding the new edges doesn't result in a cycle being
 *  formed.
 */
/*#define SCIP_DEBUG*/
#include <string.h>

#include "heur_sinks.h"
#include "scip/scip.h"
#include "pedigrees.h"
#include "parent_set_data.h"
#include "metadata.h"
#include "utils.h"
#include "scip/pub_var.h"

#define HEUR_NAME             "sinks"                          /**< The name of the heuristic. */
#define HEUR_DESC             "primal heuristic template"      /**< A description of the heuristic. */
#define HEUR_DISPCHAR         'k'                              /**< The character to display in the solver when the heuristic is used. */
#define HEUR_PRIORITY         10                               /**< The calling priority of the heuristic. */
#define HEUR_FREQ             1                                /**< The heurstic is called only in the root */
#define HEUR_FREQOFS          0                                /**< The heurstic is called from the first node onwards. */
#define HEUR_MAXDEPTH         -1                               /**< Set no depth limit for calling the heuristic. */
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPNODE   /**< Call the heuristic after each LP solve during the cut-and-price loop. */
#define HEUR_USESSUBSCIP      FALSE                            /**< The heuristic doesn't use a secondary SCIP instance. */


#define EPSILON              1e-9                              /**< How much improvement a candidate parent set to add must have before we chose it. */

#define max(A,B) ((A) > (B) ? (A) : (B))

#define DEFAULT_MAXDIVEDEPTH    100    /**< maximum dive depth when probing */
#define DEFAULT_PROBING         FALSE  /**< whether to use probing */
#define DEFAULT_ASSUMENOPOSOBJ  TRUE   /**< whether to assume that no problem variable has a positive objective coefficient */
#define DEFAULT_SEED            0      /**< the seed to use for random values. */
#define DEFAULT_NRUNS           1      /**< how often to run the heuristic on each LP solution */
#define DEFAULT_PALIM          -1      /**< only DAGS with parent sets of size at most this will be considered (-1 is unlimited) */

/** There is no copy method. */
#define heurCopySinks NULL
/** There is no deinitialization method. */
#define heurExitSinks NULL
/** There is no initialization method. */
#define heurInitsolSinks NULL
/** There is no deinitialization method. */
#define heurExitsolSinks NULL

/** Data needed by the sinks primal heuristic. */
struct SCIP_HeurData
{
   SCIP_Bool    printsols;         /**< whether to print out candidate primal solutions */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generation */
   FILE*        file;              /**< where to print candidate primal solutions, NULL for stdout */
   SCIP_Bool    probing;           /**< whether to use probing */
   int          maxdivedepth;      /**< maximum dive depth when probing */
   SCIP_Bool    assumenoposobj;    /**< whether to assume that no problem variable has a positive objective coefficient */
   int          nruns;             /**< how often to run the heuristic on each LP solution */
   int          palim;             /**< only DAGS with parent sets of size at most this will be considered (-1 is unlimited) */
   SCIP_Bool    distinguishedrep;  /**< whether only distinguished reps of MECs will be feasible */
   int          seed;              /**< initial value for random seed */
};

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSinks)
{
   /*lint --e{715}*/
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);
   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSinks)
{
   /*lint --e{715}*/
   char* filesols;
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   SCIP_CALL( SCIPgetStringParam(scip, "heuristics/sinks/filesols", &filesols) );
   if( strcmp(filesols, "") != 0 )
      heurdata->file = fopen(filesols, "w");
   else
      heurdata->file = NULL;
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/distinguishedrep", &(heurdata->distinguishedrep)) );
   SCIPheurSetData(heur, heurdata);
   return SCIP_OKAY;
}

/** execution method of primal heuristic */

/* sink finding algorithm: for each selected parentset I(v<-W)
   there is a loss of 1 - x*(v<-W) where x*(v<-W) is LP value
   of I(v<-W). Algorithm greedily attempts to select parentsets representing a DAG which minimises
   the sum of such losses
*/

static
SCIP_DECL_HEUREXEC(heurExecSinks)
{
   const ParentSetData* psd;                     /** parent set data (or NULL if not set) */

   SCIP_Real** lpsolvals;                        /** needed for storing the LP solution
                                                     asking for values from NULL using SCIPgetSolVals 
                                                     (at the wrong point) will render the pseudo-solution */

   SCIP_SOL* incumbent  = SCIPgetBestSol(scip);  /** incumbent solution (NULL if none found so far ) */
   SCIP_Real z_incumbent = -SCIPinfinity(scip);  /** objective value of incumbent solution */
   SCIP_Real z;                                  /** running value of solution under construction */


   int* bestparents = NULL;                      /** bestparents[j] is the index for the best ( highest scoring )
                                                     remaining parentset for node j */
   int* not_sinks = NULL;                        /** not_sinks[i] is the ith remaining non-sink */
   SCIP_Bool* sink = NULL;                       /** sink[i] = TRUE if i has been chosen as a sink */
   SCIP_Real* loss_lb = NULL;                    /** loss_lb[j] is a lower bound on the 'local' loss
                                                     incurred by chosing a parent set for variable j */

   int i;                                        /** indexes not_sinks array */
   int j;                                        /** indexes a BN vertex */ 
   int k;                                        /** indexes parent sets */  
   int l;                                        /** indexes parents in a parent set */
   int sinks_to_choose;                          /** number of sinks remaining to be found */
   SCIP_Bool success = FALSE;                    /** TRUE is a feasible solution good enough to keep is constructed */
   
   SCIP_SOL* sol = NULL;                         /** solution being constructed (if not probing) */
   SCIP_SOL* bestsol = NULL;                     /** best solution so far (over multiple runs) */
   SCIP_Real bestsolz = -SCIPinfinity(scip);     /** objective value of best solution so far (over multiple runs) */

   SCIP_HEURDATA* heurdata;                      /** heuristic data */
   
   SCIP_Bool rejected;                           /** TRUE if a vertex rejected as a sink because it is the parent of a vertex which is not yet a sink */
   int nfixed ;                                  /** number of family variables fixed to 1 */
   int* fixedc;                                  /** child vertices of family variables fixed to 1 */
   int* fixedp;                                  /** indices of parent sets of family variables fixed to 1 */
   int child;                                    /** an element of fixedc */
   int parent_set;                               /** an element of fixedp */

   int bestvariable;                             /** best vertex (so far) to choose as next sink */
   SCIP_VAR* bestpavar;                          /** best family variable (so far) to set to 1 since child is best vertex to choose as next sink */
   SCIP_Real bestloss;                           /** loss associated with choosing bestvariable as next sink */
   int bestindex;                                /** index of bestvariable in not_sinks array */

   SCIP_VAR* pavar;                              /** a family variable */
   SCIP_Real j_loss;                             /** loss associated with choosing a vertex 'j' as the next sink */
   SCIP_Real improvement;                        /** improvement in loss if currently considered vertex were chosen as next sink rather than current best sink */

   SCIP_Bool current_bestparents_for_j_allowed;  /** whether the current best parent set of vertex 'j' are allowed, now that a new sink has been chosen */

   SCIP_Bool cutoff;                             /** if probing, records whether probing node can be cutoff */
   SCIP_VAR** cands;                             /** if probing, stores pseudo branching candidates */
   int ncands;                                   /** if probing, stores number of pseudo branching candidates */
   SCIP_VAR* var;                                /** if probing, first pseudo branching candidate variable */
   int divedepth;                                /** if probing, records depth of dive */

   int run;                                      /** current run of the heuristic */

   int parent_to_check;                          /** a parent to check not to already be a sink */

   SCIP_VAR* arrow_edgevar;
   int h;

   SCIP_SOL* dr_sol = NULL;
   
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   psd = MD_getParentSetData(scip);
   if( psd == NULL )
   {
      SCIPdebugMessage("Couldn't find parent set data structure.\n");
      return SCIP_OKAY;
   }
   
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &heurdata->randnumgen, 
         SCIPinitializeRandomSeed(scip, heurdata->seed), TRUE) );

   
   /* only run when not in the root if distinguished rep constraint in play */
   if( !heurdata->distinguishedrep && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   /* try to create a distinguished representative for a MEC from the LP solution */
   if( heurdata->distinguishedrep )
   {
      SCIP_Bool ok = TRUE;
      SCIP_VAR** vars = SCIPgetVars(scip);
      
      for( i = 0; i < SCIPgetNVars(scip); i++ )
         if( !SCIPisFeasIntegral(scip,SCIPgetSolVal(scip,NULL,vars[i])) )
         {
            ok = FALSE;
            break;
         }

      if( ok )
      {
         /* create initial solution of all variables set to zero */
         SCIP_CALL( SCIPcreateSol(scip, &dr_sol, heur) );

         /* setting 4th argument to TRUE means that this function call
            attempts to construct a DAG which is a distinguished representative
            of a MEC
         */
         if( is_dr_feasible(scip, psd, NULL, TRUE, dr_sol) )
         {
            SCIP_CALL(SCIPtrySol(scip, dr_sol,
                  FALSE, /* don't want violations to be displayed */
                  FALSE, /* don't check violations */
                  TRUE,  /* do check bounds */
                  FALSE, /* no need to check integrality */
                  TRUE,  /* do check LP rows */
                  &success));
         }

         SCIP_CALL( SCIPfreeSol(scip, &dr_sol) );
         
         if( success )
         {
            *result = SCIP_FOUNDSOL;
            return SCIP_OKAY;
         }
      }
   }
      
   if( heurdata->assumenoposobj )
   {
      if( incumbent != NULL )
         z_incumbent = SCIPgetSolOrigObj(scip, incumbent);
      
      SCIPdebugMessage("Need to beat %f\n", z_incumbent);
   }
   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bestparents, psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &not_sinks, psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sink, psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &loss_lb, psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &lpsolvals, psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &fixedc, psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &fixedp, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(lpsolvals[i]), psd->nParentSets[i]) );
      /* at this point SCIPgetSolVals(scip,NULL,..) gets LP values */
      SCIP_CALL( SCIPgetSolVals(scip, NULL, psd->nParentSets[i], psd->PaVars[i], lpsolvals[i]) );
   }

   
   for( run = 0; run < heurdata->nruns; ++run )
   {
   
      /* Initialise */
      z = 0;
      sinks_to_choose = psd->n;
      for( i = 0; i < sinks_to_choose; ++i )
      {
         /* assume parent sets are ordered with best score first, so: */
         bestparents[i] = 0;

         /* initially no vertex has been chosen as a sink */
         not_sinks[i] = i;
         sink[i] = FALSE;

         loss_lb[i] = 0.0;
      }

      /* shuffle not_sinks array */
      SCIPrandomPermuteIntArray(heurdata->randnumgen,not_sinks,0,sinks_to_choose);

      /* create initial solution of all variables set to zero */
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

      if( heurdata->probing )
      {
         if( SCIPinProbing(scip) )
            SCIP_CALL( SCIPendProbing(scip) );
         SCIP_CALL( SCIPstartProbing(scip) );
      }

      /* rule out best scoring parents which are already fixed to zero */
      /* or those which exceed palim */
      nfixed = 0;
      for( j = 0; j < psd->n; ++j )
      {
         while( bestparents[j] < psd->nParentSets[j] && (
            SCIPvarGetUbLocal(psd->PaVars[j][bestparents[j]]) < 0.5 
            || (heurdata->palim != -1 && psd->nParents[j][bestparents[j]] > heurdata->palim)) )
         {
            loss_lb[j] += lpsolvals[j][bestparents[j]];
            bestparents[j]++;
         }
         
         /* if somehow have already ruled out all parent sets
          then abort
         */
         if( bestparents[j] == psd->nParentSets[j] )
            goto ABORT;
            
         /* record whether already fixed to 1 */
         if( SCIPvarGetLbLocal(psd->PaVars[j][bestparents[j]]) > 0.5 )
         {
            fixedc[nfixed] = j;
            fixedp[nfixed++] = bestparents[j];
         }
      }
   
      while( sinks_to_choose )
      {
         /* look for a new sink */
         bestvariable = -1;
         bestpavar = NULL;
         bestloss = SCIPinfinity(scip);
         bestindex = -1;
      
         for( i = 0; i < sinks_to_choose; ++i )
         {
         
            /* consider variable j as a new sink */
            j = not_sinks[i];
            assert(j > -1);
            assert(j < psd->n);
            assert(!sink[j]);
            assert(bestparents[j] > -1);
            assert(bestparents[j] < psd->nParentSets[j]);

            pavar = psd->PaVars[j][bestparents[j]];
            SCIPdebugMessage("considering variable %s as a sink.\n", SCIPvarGetName(pavar));
         
            /* loss for choosing j is how much bigger actual loss (1 - lpsolvals[j][bestparents[j]])
               is than the best we can hope for ( by this stage ) for j (ie loss_lb[j])
            */

            /* can't be a sink if it is the parent of something fixed to 1 which is not yet a sink .. */
            /* dubious use of 'k' */
            rejected = FALSE;
            for( k = 0; k < nfixed; ++k )
            {
               child = fixedc[k];
            
               if( !sink[child] )
               {
                  parent_set = fixedp[k];
                  for( l = 0; l < psd->nParents[child][parent_set]; ++l )
                  {
                     if( psd->ParentSets[child][parent_set][l] == j )
                     {
                        rejected = TRUE;
                        break;
                     }
                  }
                  if( rejected )
                     break;
               }
            }

            if( rejected )
            {
               SCIPdebugMessage("variable %s rejected as sink since parent of node %s.\n", SCIPvarGetName(pavar), psd->nodeNames[child]);
               continue;
            }

            /* assert(SCIPisGE(scip, (1 - lpsolvals[j][bestparents[j]]), loss_lb[j])); */
            j_loss = ( (1 - lpsolvals[j][bestparents[j]]) - loss_lb[j] );

            /* how much better is the loss for j compared to current best choice? */
            improvement = bestloss - j_loss;
            if( improvement > EPSILON )
            {
               SCIPdebugMessage("variable %s is best sink so far.\n", SCIPvarGetName(pavar));
               bestindex = i;
               bestloss = j_loss;
               bestvariable = j;
               bestpavar = pavar;

            }
            else
            {
               SCIPdebugMessage("variable %s rejected as sink.\n", SCIPvarGetName(pavar));
            }
         }

         /* if all rejected, since perhaps sent some weird LP sol, then give up */
         if( bestpavar == NULL )
            goto ABORT;

         assert(SCIPvarGetUbLocal(bestpavar) > 0.5);
         assert(SCIPisLE(scip, bestloss, 1.0));
         assert(bestindex > -1);
         assert(bestindex < psd->n);
         assert(bestvariable > -1);
         assert(bestvariable < psd->n);

         SCIPdebugMessage("Fixing variable <%s>[%g,%g] = 1\n",
            SCIPvarGetName(bestpavar), SCIPvarGetLbLocal(bestpavar), SCIPvarGetUbLocal(bestpavar));

         if( heurdata->probing )
         {
            if( SCIPvarIsActive(bestpavar) &&
               ( SCIPvarGetStatus(SCIPvarGetTransVar(bestpavar)) != SCIP_VARSTATUS_MULTAGGR  || SCIPvarGetMultaggrNVars(SCIPvarGetTransVar(bestpavar)) == 1 ))
               SCIP_CALL( SCIPfixVarProbing(scip, bestpavar, 1.0) );
            SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );

            if( cutoff )
            {
               SCIPdebugMessage("Could not fix variable %s to 1.\n", SCIPvarGetName(bestpavar));
               /* skip this one and try again */
               loss_lb[bestvariable] += lpsolvals[bestvariable][bestparents[bestvariable]];
               bestparents[bestvariable]++;
               while( bestparents[bestvariable] < psd->nParentSets[bestvariable] && SCIPvarGetUbLocal(psd->PaVars[bestvariable][bestparents[bestvariable]]) < 0.5 )
               {
                  loss_lb[bestvariable] += lpsolvals[bestvariable][bestparents[bestvariable]];
                  bestparents[bestvariable]++;
               }
               /* if no choices left for bestvariable then abort */
               if( !(bestparents[bestvariable] < psd->nParentSets[bestvariable]) )
                  goto ABORT;
            }
         }
         else
         {
            /* set relevant family, arrow and edge variables to 1, abort if any are fixed to 0 */ 
            assert( bestpavar != NULL );
            if( SCIPvarIsActive(bestpavar) && SCIPvarGetStatus(SCIPvarGetTransVar(bestpavar)) != SCIP_VARSTATUS_MULTAGGR  )
            {
               SCIP_CALL( SCIPsetSolVal(scip, sol, bestpavar, 1.0) );
            }
            k = bestparents[bestvariable];
            for( l = 0; l < psd->nParents[bestvariable][k]; ++l )
            {
               for(h = 0; h < 2; h++)
               {
                  arrow_edgevar = get_arrowedge(psd,bestvariable,psd->ParentSets[bestvariable][k][l],(h==0));
                  if( arrow_edgevar == NULL || SCIPvarGetUbLocal(arrow_edgevar) < 0.5 ) 
                     goto ABORT;
                  if( SCIPvarIsActive(arrow_edgevar) &&
                     SCIPvarGetStatus(arrow_edgevar) != SCIP_VARSTATUS_MULTAGGR &&
                     SCIPvarGetStatus(SCIPvarGetTransVar(arrow_edgevar)) != SCIP_VARSTATUS_MULTAGGR &&
                     SCIPvarGetStatus(SCIPvarGetTransVar(arrow_edgevar)) != SCIP_VARSTATUS_AGGREGATED )
                  {
                     SCIP_CALL( SCIPsetSolVal(scip, sol, arrow_edgevar, 1.0) );
                  }
               }
            }
         }
      
         if( heurdata->assumenoposobj )
         {
            z += SCIPvarGetObj(bestpavar);
            if( z < z_incumbent )
            {
               SCIPdebugMessage("Aborting heuristic since objective value of partially constructed solution is %g and current incumbent has objective of %g\n", z, z_incumbent);
               goto ABORT;
            }
         }

         sink[bestvariable] = TRUE;
         SCIPdebugMessage("variable %s chosen as a sink.\n", SCIPvarGetName(bestpavar));

         /* best variable now a sink */
         /* remove from list of not_sinks by overwriting its entry with last entry in this array */
         /* decrement sinks_to_choose as well */
         sinks_to_choose--;
         not_sinks[bestindex] = not_sinks[sinks_to_choose];

         /* update allowed parent sets */
         /* and update upper bound on distance moved by any future rounding */
         for( i = 0; i < sinks_to_choose; ++i )
         {
         
            j = not_sinks[i];
            do
            {
               /* if local upper bound on 'bestparents' less than 1 then
                  can never be chosen and thus not best
                  otherwise have to check that choice of sink has
                  not ruled it out
               */
            
               if( SCIPvarGetUbLocal(psd->PaVars[j][bestparents[j]]) < 0.5 )
                  current_bestparents_for_j_allowed = FALSE;
               else
               {
                  current_bestparents_for_j_allowed = TRUE;
                  /* check each parent of current best parent set for j */
                  for( l = 0; l < psd->nParents[j][bestparents[j]]; ++l )
                  {
                     parent_to_check = psd->ParentSets[j][bestparents[j]][l];
                     /* if parent_to_check is already a sink it cannot be a parent of j (which is currently a non-sink) */
                     if( sink[parent_to_check] )
                     {
                        /* current best parent set no longer allowed */
                        current_bestparents_for_j_allowed = FALSE;
                        break;
                     }
                  }
               }

               if( !current_bestparents_for_j_allowed )
               {
                  if( (!heurdata->probing && SCIPvarGetLbLocal(psd->PaVars[j][bestparents[j]]) > 0.5) || (heurdata->probing && cutoff) )
                  {
                     SCIPdebugMessage("Could not rule out parent set %s, aborting heuristic.\n", SCIPvarGetName(psd->PaVars[j][bestparents[j]]));
                     goto ABORT;
                  }

                  /* if not probing, check that ruled out parents for j have value 0 */
                  if( !heurdata->probing && SCIPgetSolVal(scip, sol, psd->PaVars[j][bestparents[j]]) == 1 ) 
                  {
                     /* ruled out parent set already selected (!) due to setting arrow/edge variables */
                     /* have to abort */
                     goto ABORT;
                  }
                  /* update the least we can expect to lose by rounding some other, future parent set to 1 */
                  /* for future choices we just care about how much *more* is lost */
                  loss_lb[j] += lpsolvals[j][bestparents[j]];
                  /* update which is the best parent set for j */
                  bestparents[j]++;
                  while( bestparents[j] < psd->nParentSets[j] && 
                     ( SCIPvarGetUbLocal(psd->PaVars[j][bestparents[j]]) < 0.5  || 
                        ( heurdata->palim != -1 && psd->nParents[j][bestparents[j]] > heurdata->palim) ) )
                  {
                     SCIPdebugMessage("This variable %s set to zero.\n", SCIPvarGetName(psd->PaVars[j][bestparents[j]]));
                     loss_lb[j] += lpsolvals[j][bestparents[j]];
                     bestparents[j]++;
                  }
                  /* if no choices left for j then abort */
                  if( !(bestparents[j] < psd->nParentSets[j]) )
                  {
                     SCIPdebugMessage("No parent set choices left for node %s, aborting heuristic.\n", psd->nodeNames[j]);
                     goto ABORT;
                  }
               }
            
            }
            while( !current_bestparents_for_j_allowed );
         }   /* all allowed parent set now updated */
      }  /* main loop finished, heuristic solution now constructed */

      if( bestsol == NULL || z > bestsolz )
      {
         SCIPcreateSolCopy(scip,&bestsol,sol);
         bestsolz = z;
      }
      SCIP_CALL( SCIPfreeSol(scip, &sol) );
      sol = NULL;

   }

   if( heurdata->probing )
   {
      /* mop up any unfixed variables */
      
      SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );
      
      /* similar code to heur_fixandinfer.c: author Tobias Achterberg */
      
      SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, &ncands, NULL) );
      
      divedepth = 0;
      
      while( !cutoff && ncands > 0 && divedepth < heurdata->maxdivedepth && !SCIPisStopped(scip) )
      {
         /* fix first unfixed variable to its current pseudo solution value */
         var = cands[0];
         SCIP_CALL( SCIPnewProbingNode(scip) );
         divedepth++;
         SCIP_CALL( SCIPfixVarProbing(scip, var, SCIPgetVarSol(scip, var)) );
         SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );
         if( !cutoff )
            SCIP_CALL( SCIPgetPseudoBranchCands(scip, &cands, &ncands, NULL) );
      }
      
      if( !cutoff )
      {
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPprintSol(scip, NULL, NULL, FALSE) );
#endif
         SCIP_CALL( SCIPtryCurrentSol(scip, heur, FALSE, FALSE, FALSE, TRUE, &success) );
      }
   }
   else  /* not probing */
   {
      if( heurdata->printsols )
      {
         FILE* output = heurdata->file == NULL ? stdout : heurdata->file;
         fprintf(output, "START sink heuristic solution.\n");
         SCIP_CALL( SCIPprintSol(scip, bestsol, heurdata->file, FALSE) );
         fprintf(output, "END sink heuristic solution.\n");
      }
#ifdef SCIP_DEBUG
      /* trysolFree ( if using it ) clears sol so print it now */
      SCIP_CALL( SCIPprintSol(scip,  bestsol, NULL, FALSE) );
#endif
      
      SCIP_CALL(SCIPtrySol(scip, bestsol,
            FALSE, /* don't want violations to be displayed */
            FALSE, /* don't check violations */
            TRUE,  /* do check bounds */
            FALSE, /* no need to check integrality */
            TRUE,  /* do check LP rows */
            &success));

      if( !success && heurdata->distinguishedrep )
      {
         /* create initial solution of all variables set to zero */
         SCIP_CALL( SCIPcreateSol(scip, &dr_sol, heur) );
         if( is_dr_feasible(scip, psd, bestsol, TRUE, dr_sol) )
            SCIP_CALL(SCIPtrySolFree(scip, &dr_sol,
                  FALSE, /* don't want violations to be displayed */
                  FALSE, /* don't check violations */
                  TRUE,  /* do check bounds */
                  FALSE, /* no need to check integrality */
                  TRUE,  /* do check LP rows */
                  &success));
         /* if( !success ) */
         /* { */
         /*    SCIP_CALL( SCIPprintSol(scip,  bestsol, NULL, FALSE) ); */
         /*    printf("*********\n"); */
         /*    SCIP_CALL( SCIPprintSol(scip,  dr_sol, NULL, FALSE) ); */
         /* } */
      }
   }
   
ABORT:
   if( success )
   {
      *result = SCIP_FOUNDSOL;
      SCIPdebugMessage("ok\n");
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMessage("not ok\n");
   }
   
   SCIPfreeBlockMemoryArray(scip, &bestparents, psd->n);
   SCIPfreeBlockMemoryArray(scip, &not_sinks, psd->n);
   SCIPfreeBlockMemoryArray(scip, &sink, psd->n);
   SCIPfreeBlockMemoryArray(scip, &loss_lb, psd->n);
   SCIPfreeBlockMemoryArray(scip, &fixedc, psd->n);
   SCIPfreeBlockMemoryArray(scip, &fixedp, psd->n);

   for( i = 0; i < psd->n; ++i )
      SCIPfreeBlockMemoryArray(scip, &(lpsolvals[i]), psd->nParentSets[i]);
   SCIPfreeBlockMemoryArray(scip, &lpsolvals, psd->n);
   
   if( heurdata->probing )
      SCIP_CALL( SCIPendProbing(scip) );
   else
   {
      if( sol != NULL )
         SCIP_CALL( SCIPfreeSol(scip, &sol) );
      if( bestsol != NULL )
         SCIP_CALL( SCIPfreeSol(scip, &bestsol) );
   }
   
   return SCIP_OKAY;
}

/** creates the sinks primal heuristic and includes it in SCIP.
 *  @return SCIP_OKAy if the operation suceeded, or an appropriate error message otherwise.
 */
SCIP_RETCODE HS_includePrimal(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_HEURDATA* heurdata;

   /* create sinks primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL(SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                             HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
                             heurCopySinks,
                             heurFreeSinks, heurInitSinks, heurExitSinks,
                             heurInitsolSinks, heurExitsolSinks, heurExecSinks,
                             heurdata));

   /* add sinks primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   SCIP_CALL(SCIPaddBoolParam(
                scip,
                "heuristics/sinks/printsols",
                "whether to print *every* BN found by sink heuristic (in SCIP solution format)",
                &heurdata->printsols,
                FALSE,
                FALSE,
                NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(
                scip,
                "heuristics/sinks/probing",
                "whether to use probing",
                &heurdata->probing,
                FALSE,
                DEFAULT_PROBING,
                NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(
                scip,
                "heuristics/sinks/maxdivedepth",
                "maximum dive depth when using probing",
                &heurdata->maxdivedepth,
                FALSE,
                DEFAULT_MAXDIVEDEPTH, 0, INT_MAX,
                NULL, NULL));


   SCIP_CALL(SCIPaddIntParam(
                scip,
                "heuristics/sinks/seed",
                "initial value for random seed",
                &heurdata->seed,
                FALSE,
                DEFAULT_SEED, 0, INT_MAX,
                NULL, NULL));


   SCIP_CALL(SCIPaddIntParam(
                scip,
                "heuristics/sinks/nruns",
                "how many times to run the heuristic on each LP solution",
                &heurdata->nruns,
                FALSE,
                DEFAULT_NRUNS, 0, INT_MAX,
                NULL, NULL));


   SCIP_CALL(SCIPaddBoolParam(
                scip,
                "heuristics/sinks/assumenoposobj",
                "whether to assume that no problem variable has a positive objective coefficient",
                &heurdata->assumenoposobj,
                FALSE,
                DEFAULT_ASSUMENOPOSOBJ,
                NULL, NULL));


   SCIP_CALL(SCIPaddStringParam(
                scip,
                "heuristics/sinks/filesols",
                "where to print solutions found by sink heuristic",
                NULL,
                FALSE,
                "",
                NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(
                scip,
                "heuristics/sinks/palim",
                "only DAGS with parent sets of size at most this will be considered (-1 is unlimited)",
                &heurdata->palim,
                FALSE,
                DEFAULT_PALIM, -1, INT_MAX,
                NULL, NULL));


   return SCIP_OKAY;
}
