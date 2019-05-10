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

/** @file cons_dagcluster.c
 *  @brief  constraint handler for acyclicity constraints
 *  @author James Cussens
 *  @author Mark Bartlett
 *  @author Tom Pottage
 *
 *  Implements a constraint for preventing cycles in the graph.
 *
 *  The constraint states that for any group of k nodes, at least one node must have no parents in
 *  that cluster.  This generalises such that there must also be at least two nodes with at most one parent,
 *  at least three with at most two parents and so on.
 *
 *  As there are exponentially many of these constraints, it is implemented as a series of cutting planes.
 *
 *  In the price-and-cut loop, cluster cutting planes are sought first (due to high SEPAPRIORITY), then other separators
 *  e.g. Gomory may kick in. Once the price-and-cut loop is finished, cluster cutting planes are looked for again
 *  (due to high ENFOPRIORITY). This may succeed since eg Gomory cuts will have 'moved' the LP solution.
 *  If one is found then LP solving is re-invoked (p35 Achterberg) which may lead to yet more cluster
 *  cutting planes.
 *
 *  Doxygen documentation for the locally defined structs @c SCIP_ConsData and @c SCIP_ConshdlrData is unfortunately not available
 *  since structs with the same name are defined in other constraint handlers and Doxygen cannot handle the name clash. You will have to consult
 *  the source code to view this documentation.
 */

/* This file was created by editing the constraint handler template file in SCIP */
/*#define SCIP_DEBUG*/
#include "cons_dagcluster.h"
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "scip/scipdefplugins.h"
#include "parent_set_data.h"
#include "solution_info.h"
#include "circuit_cuts.h"
#include "fractional_circuit_cuts.h"
#include "convex4_cuts.h"
#include "set_packing_cuts.h"
#include "subip_cuts.h"
#include "matroid_cuts.h"
#include "utils.h"

#define DEFAULT_KMAX 1
#define DEFAULT_KMAX_ROOT 1
#define DEFAULT_CLUSTERCUTS_ADDTOPOOL TRUE
#define DEFAULT_PROPAGATESPC TRUE
#define DEFAULT_USEINCUMBENTCONS FALSE
#define DEFAULT_SEPASOL_USECONVEX4 FALSE
#define DEFAULT_SUBIPALWAYS TRUE
#define DEFAULT_SUBIPEVER TRUE
#define DEFAULT_SUBIPDEPTHLIMIT -1
#define DEFAULT_FRACTIONALALWAYS FALSE
#define DEFAULT_FRACTIONALEVER FALSE
#define DEFAULT_FORCECUTS TRUE
#define DEFAULT_EXTRAPROPS FALSE
#define DEFAULT_MATROIDPAIRCUTS FALSE
#define DEFAULT_MATROIDPAIRMAXCUTS 20
#define DEFAULT_MATROIDPAIRCUTSLIMIT 6
#define DEFAULT_MATROIDCUTS FALSE
#define DEFAULT_CIRCUIT_SIZE_LIM 3
#define DEFAULT_GROUND_SET_SIZE_LIM 8
#define DEFAULT_ADDSPC TRUE


#define DEFAULT_CONSSEPALPSUBIPTIMELIMIT   10 /*1e+20*/
#define DEFAULT_CONSSEPASOLSUBIPTIMELIMIT  10 /*1e+20*/

#define DEFAULT_CONSSEPALPSUBIPGAPLIMIT  0
#define DEFAULT_CONSSEPASOLSUBIPGAPLIMIT 0

#define DEFAULT_CONSSEPALPSUBIPABSGAPLIMIT  0
#define DEFAULT_CONSSEPASOLSUBIPABSGAPLIMIT 0

#define DEFAULT_ROOTGAPSEPALIMIT  0.001

#define CHECK_WITH_ARROWS TRUE

#define EPSILON 0.0001
#define min(A,B) ((A) > (B) ? (B) : (A))
#define max(A,B) ((A) > (B) ? (A) : (B))

/* for SCIP 3.0.0 */
#define consGetVarsDagcluster NULL
#define consGetNVarsDagcluster NULL

/* Constraint handler properties */
#define CONSHDLR_NAME            "dagcluster"          /**< Name of the constraint handler. */
#define CONSHDLR_DESC            "DAG cluster-based acyclicity constraint handler" /**< Description of the constraint handler. */
#define CONSHDLR_SEPAPRIORITY      100000000           /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       90                /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY      -900000            /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ                10            /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ                 1            /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ              100            /**< frequency for using all instead of only the useful constraints in separation, propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS             -1           /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA             FALSE           /**< should separation method be delayed, if other separators found cuts?, */
#define CONSHDLR_DELAYPROP             FALSE           /**< should propagation method be delayed, if other propagators found reductions?, */
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_ALWAYS
#define CONSHDLR_NEEDSCONS              TRUE           /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler, */


/* Data structures */
/** constraint data for dagcluster constraints
    Note that 'local', 'global' and 'sol' info
    is constraint (not constraint handler) specific
    since all depend on which variables are included in the constraint
*/
struct SCIP_ConsData
{
   ParentSetData* psd;
   SCIP_Bool*** store;                /**< store[i][k][j] = TRUE if j is in kth parent set for i */
   CircuitCutsStorage* ccs;
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int kmax;                        /**< maximum value of k when searching for k-cluster inequalities (outside the root) */
   int kmaxroot;                    /**< maximum value of k when searching for k-cluster inequalities (in the root) */
   SCIP_Bool clustercuts_addtopool; /**< whether to add cluster cuts to the global cut pool */
   SCIP_Bool propagatespc;          /**< whether to propagate added set packing constraints */
   SCIP_Bool nosep_dagcluster;      /**< whether the dagcluster separator failed to find an efficacious cut when last called on LP solution */
   SCIP_Bool nosepasol_dagcluster;  /**< whether the dagcluster separator failed to find an efficacious cut when last called on an arbitrary primal solution */
   SCIP_Bool useincumbentcons;      /**< whether to (additionally) look for cluster cuts which are tight for the current incumbent */
   SCIP_Bool sepasol_useconvex4;    /**< whether to search for convex4 cuts when separating an arbitrary primal solution */
   SCIP_Bool subipalways;           /**< whether to always use the subIP to look for cluster cuts */
   SCIP_Bool subipever;             /**< whether to ever use the subIP to look for cluster cuts */
   int subipdepthlimit;             /**< maximum depth for using the subIP to search for cluster cuts (-1 for no limit) */
   SCIP_Bool fractionalalways;      /**< whether to always use fractional cycles to look for cluster cuts */
   SCIP_Bool fractionalever;        /**< whether to ever use fractional cycles to look for cluster cuts */
   SCIP_Bool forcecuts;             /**< whether to force cuts to be added */
   SCIP_Bool extraprops;            /**< whether to do extra (slow) propagations */

   int sepalp_calls;                /**< number of calls to separate the LP solution */
   int writemip_minsepalp_calls;    /**< only write out pre-separation MIP relaxation after reaching this number of calls to LP separator */
   int writemip_maxsepalp_calls;    /**< only write out pre-separation MIP relaxation if not beyond this number of calls to LP separator */

   SCIP_Real conssepalpsubiptimelimit;
   SCIP_Real conssepasolsubiptimelimit;
   SCIP_Real consenfolpsubiptimelimit;

   SCIP_Real conssepalpsubipgaplimit;
   SCIP_Real conssepasolsubipgaplimit;
   SCIP_Real consenfolpsubipgaplimit;

   SCIP_Real conssepalpsubipabsgaplimit;
   SCIP_Real conssepasolsubipabsgaplimit;
   SCIP_Real consenfolpsubipabsgaplimit;

   SCIP_Real rootgapsepalimit;

   SCIP_Bool matroidpaircuts;
   int matroidpaircutslimit;
   int matroidpairmaxcuts;
   
   SCIP_Bool matroidcuts;
   int circuit_size_lim;
   int ground_set_size_lim;

   SCIP_Bool addspc;
   
};

/** Checks whether a solution is feasible for a given constraint.
 *  Uses arrow variables only to check for acyclicity
 *  Solution is infeasible only if there is a cycle involving arrows with positive
 *  values in the solution
 * 
 *  @param scip The SCIP instance conatining the constraint and solution.
 *  @param cons The acyclicity constraint to check.
 *  @param The solution to check for feasibility.  If NULL is given, then the
 *         current solution will be used.
 *  @return TRUE if the solution is feasible for this constraint, or FALSE otherwise.
 */
static
SCIP_Bool isFeasible_arrows(
   SCIP* scip,
   SCIP_CONS* cons,
   SCIP_SOL* sol
   )
{

   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   const int n = consdata->psd->n;
   SCIP_Bool** ancestor;            /* ancestor[i][j] = TRUE iff j is an ancestor of i in current solution */
   SCIP_Bool* orphan;               /* orphan[i] = TRUE iff there is no j for which ancestor[i][j] is TRUE */

   int i;
   int j;
   int kk;

   SCIP_VAR* arrow_i_j;
   SCIP_Bool result = TRUE;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &orphan, n) );   
   SCIP_CALL( SCIPallocBufferArray(scip, &ancestor, n) );

   for( i = 0; i < n; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(ancestor[i]), n) );
      orphan[i] = TRUE;
      for( j = 0; j < n; ++j )
      {
         if( i == j)
            continue;
         
         arrow_i_j = get_arrow(consdata->psd,i,j);
         if( arrow_i_j != NULL && SCIPisPositive(scip, SCIPgetSolVal(scip, sol, arrow_i_j)) )
         {
            ancestor[i][j] = TRUE;
            orphan[i] = FALSE;
         }
         else
            ancestor[i][j] = FALSE;
      }
   }

   /* no need to check for 2-cycles since constraint arising from edge variables 
      takes care of that 
   */
   
   /* compute transitive closure (with early exit if a cycle detected) */

   /* Cormen et al:
      ancestor[i][j][k] = TRUE iff there exists
      a path from i to j and any intermediate vertices
      are in {0, ..., k-1}
      But drop k index and just update ancestor[i][j]
   */
   
   for( kk = 0; kk < n; ++kk )
   {
      if( orphan[kk] ) 
         continue; 
      
      for( i = 0; i < n; ++i )
      {
         if( i == kk || !ancestor[i][kk] ) 
            continue; 
         
         for( j = 0; j < n; ++j )
         {
            /* only care about cycles so don't bother
               finding descendants of orphans, 
               hence the "orphan[j]" on next line
            */
            if( orphan[j] || j == kk || j == i )
               continue;
            if( ancestor[kk][j] )
            {
               ancestor[i][j] = TRUE;
               if( ancestor[j][i] )
               {
                  result = FALSE;
                  goto TERMINATE;
               }
            }
         }
      }
   }

TERMINATE:
   for( i = n-1; i > -1; --i )
      SCIPfreeBufferArray(scip, &(ancestor[i]));
   SCIPfreeBufferArray(scip, &ancestor);
   SCIPfreeBufferArray(scip, &orphan);

   return result;
}


/** Checks whether a solution is feasible for a given constraint.
 *
 *  Uses family variables only to check for acyclicity
 *  The graph to check is is constructed by finding, for each vertex, the first
 *  parent set with positive value in the solution.
 *  This function gives the same answer as isFeasible_arrows if there is no more
 *  than one positive parent set for each vertex
 *
 *  @param scip The SCIP instance conatining the constraint and solution.
 *  @param cons The acyclicity constraint to check.
 *  @param The solution to check for feasibility.  If NULL is given, then the
 *         current solution will be used.
 *  @return TRUE if the solution is feasible for this constraint, or FALSE otherwise.
 */
static
SCIP_Bool isFeasible_fvs(
   SCIP* scip,
   SCIP_CONS* cons,
   SCIP_SOL* sol
   )
{
   int i;
   int ii;
   int k;
   int l;
   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   int n = consdata->psd->n;
   int** parents_of;
   int* num_parents_of;
   SCIP_Bool* dealt_with;
   int* todo;
   SCIP_Bool made_progress = TRUE;
   int n_todo = 0;
   int n_new_todo;
   int n_new_parents;
   int* new_parents;
   int* new_todo;
   int* tmp;

   SCIP_Bool result = TRUE;
   ParentSetData* psd = NULL;
   int parent;
   
   psd = consdata->psd;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &num_parents_of, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &parents_of, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &new_parents, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dealt_with, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &todo, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &new_todo, n) );

   /* find which parent set is positive */
   for( i = 0; i < n; i++ )
   {
      /* no parents found for i initially */
      dealt_with[i] = TRUE;
      parents_of[i] = NULL;
      for( k = 0; k < psd->nParentSets[i]; k++ )
         if( SCIPisPositive(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k])) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(parents_of[i]), psd->nParents[i][k]) );
            /* need to copy parent set, since this will be altered in this function */
            num_parents_of[i] = 0;
            for( l = 0; l < psd->nParents[i][k]; l++ )
            {
               parent = psd->ParentSets[i][k][l];
               assert( parent > -1 );
               if( parent < n )
               {
                  /* be robust to the possibility of 'extra' vertices appearing as parent,
                     treat them as source vertices */
                  parents_of[i][num_parents_of[i]++] = parent;
               }
            }
            if( num_parents_of[i] > 0 )
            {
               dealt_with[i] = FALSE;
               todo[n_todo++] = i;
            }
            break;
         }
   }

   /* this is basically topological sort */
   /* specifically: repeatedly remove those parents from parent sets
      which have been 'dealt_with'. A vertex is 'dealt_with' if all
      its parents are 'dealt_with'.
   */
   while( made_progress )
   {
      made_progress = FALSE;
      n_new_todo = 0;
      for( ii = 0; ii < n_todo; ++ii )
      {
         i = todo[ii];
         n_new_parents = 0;
         for( l = 0; l < num_parents_of[i]; ++l )
            if( !dealt_with[parents_of[i][l]] )
               new_parents[n_new_parents++] = parents_of[i][l];
         if( n_new_parents == 0 )
         {
            dealt_with[i] = TRUE;
            made_progress = TRUE;
         }
         else
         {
            if( n_new_parents < num_parents_of[i] )
            {
               /* have to copy (possibly reduced) parent set over */
               for( l = 0; l < n_new_parents; ++l )
                  parents_of[i][l] = new_parents[l];
               num_parents_of[i] = n_new_parents;
            }
            new_todo[n_new_todo++] = i;
         }
      }
      tmp = todo;
      todo = new_todo;
      new_todo = tmp;
      n_todo = n_new_todo;
   }

   if( n_todo == 0 )
      result = TRUE;
   else
      result = FALSE;

   for( i = n - 1; i > -1; i-- )
      if( parents_of[i] != NULL )
         SCIPfreeBufferArray(scip, &(parents_of[i]));
   SCIPfreeBufferArray(scip, &new_todo);
   SCIPfreeBufferArray(scip, &todo);
   SCIPfreeBufferArray(scip, &dealt_with);
   SCIPfreeBufferArray(scip, &new_parents);
   SCIPfreeBufferArray(scip, &parents_of);
   SCIPfreeBufferArray(scip, &num_parents_of);

   return result;
}

/** Checks whether a solution is feasible for a given set of constraints.
 *
 *  A solution is feasible for a constraint if
 *  - it has no more than 1 parent set for each node involved in the constraint, and
 *  - there are no cycles formed by the nodes involved in the constraint.
 *
 *  @param scip The SCIP instance conatining the constraints and solution.
 *  @param conss The acyclicity constraints to check.
 *  @param nconss The number of acyclicity constraints to check.
 *  @param sol The solution to check for feasibility.  If NULL is given, then the
 *         current solution will be used.
 *  @return TRUE if the solution is feasible for all the constraints, or FALSE otherwise.
 */
static
SCIP_Bool areFeasible(
   SCIP* scip,
   SCIP_CONS** conss,
   int nconss,
   SCIP_SOL* sol
   )
{
   int c;
   if( CHECK_WITH_ARROWS )
   {
      for( c = 0; c < nconss; c++ )
         if( !isFeasible_arrows(scip, conss[c], sol) )
            return FALSE;
   }
   else
   {
      for( c = 0; c < nconss; c++ )
         if( !isFeasible_fvs(scip, conss[c], sol) )
            return FALSE;
   }
   return TRUE;
}

/** Creates the data for a constraint.
 *  @param scip The SCIP instance to which the constraint belongs.
 *  @param consdata The location to store the new constraint data.
 *  @param psd The parent set data on which the constraint is based.
 *  @return SCIP_OKAY if successful, or an appropriate error otherwise.
 */
static
SCIP_RETCODE createConsData(
   SCIP* scip,
   SCIP_CONSDATA** consdata,
   ParentSetData* psd
   )
{
   int i;
   int k;
   int l;

   assert( scip != NULL);
   assert( consdata != NULL);
   assert( psd != NULL);
   
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( PS_copyParentSetData(scip, psd, &((*consdata)->psd)) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store), (*consdata)->psd->n) );
   for( i = 0; i < (*consdata)->psd->n; ++i )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store[i]), (*consdata)->psd->nParentSets[i]) );
   for( i = 0; i < (*consdata)->psd->n; ++i )
      for( k = 0; k < (*consdata)->psd->nParentSets[i]; ++k )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store[i][k]), (*consdata)->psd->n) );
         for( l = 0; l < (*consdata)->psd->n; ++l )
            (*consdata)->store[i][k][l] = FALSE;
         for( l = 0; l < (*consdata)->psd->nParents[i][k]; ++l )
            (*consdata)->store[i][k][(*consdata)->psd->ParentSets[i][k][l]] = TRUE;
      }

   SCIP_CALL( SCIPallocBlockMemory(scip, &((*consdata)->ccs)) );
   CC_initialise(scip, (*consdata)->ccs, (*consdata)->psd);
   /* following function not currently used */
   /* FC_initialise(scip, (*consdata)->ccs, (*consdata)->psd); */

   return SCIP_OKAY;
}

#if 0
// Correctly creates a transformed consdata struct from the original psd structure
static
SCIP_RETCODE createTransConsData(
   SCIP* scip,
   SCIP_CONSDATA** consdata,
   ParentSetData* psd
   )
{
   int i;
   int k;
   int l;

   assert( scip != NULL);
   assert( consdata != NULL);
   assert( psd != NULL);
   
   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( PS_transformParentSetData(scip, psd, &((*consdata)->psd)) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store), (*consdata)->psd->n) );
   for( i = 0; i < (*consdata)->psd->n; ++i )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store[i]), (*consdata)->psd->nParentSets[i]) );
   for( i = 0; i < (*consdata)->psd->n; ++i )
      for( k = 0; k < (*consdata)->psd->nParentSets[i]; ++k )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store[i][k]), (*consdata)->psd->n) );
         for( l = 0; l < (*consdata)->psd->n; ++l )
            (*consdata)->store[i][k][l] = FALSE;
         for( l = 0; l < (*consdata)->psd->nParents[i][k]; ++l )
            (*consdata)->store[i][k][(*consdata)->psd->ParentSets[i][k][l]] = TRUE;
      }

   SCIP_CALL( SCIPallocBlockMemory(scip, &((*consdata)->ccs)) );
   CC_initialise(scip, (*consdata)->ccs, (*consdata)->psd);
   /* following function not currently used */
   /* FC_initialise(scip, (*consdata)->ccs, (*consdata)->psd); */

   return SCIP_OKAY;
}
#endif

/** Try to split each of the locally active constraints into multiple constraints locally
 *  by looking for strongly connected components in the data.
 *
 *  @param scip The SCIP instance that the constraints belong to.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error otherwise.
 */
SCIP_RETCODE DC_tryToSplit(
   SCIP* scip
   )
{
   int i;
   int j;
   SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   int num_active_cons = SCIPconshdlrGetNActiveConss(conshdlr);
   SCIP_CONS** original_conss = SCIPconshdlrGetConss(conshdlr);
   SCIP_NODE* current_node = SCIPgetCurrentNode(scip);

   /* I don't know if this is actually needed or not. */
   /* We add and 'delete' cons from the node while */
   /* iterating over conss, so I make a copy to iterate */
   /* over just to be sure. */
   SCIP_CONS** conss;
   SCIP_CALL( SCIPallocMemoryArray(scip, &conss, num_active_cons) );
   for( i = 0; i < num_active_cons; i++ )
      conss[i] = original_conss[i];

   /*{
      SCIP_NODE* parent = SCIPnodeGetParent(current_node);
      long long int parent_num = -1;
      if (parent != NULL)
         parent_num = SCIPnodeGetNumber(parent);
      printf("*** Trying to split %lld(%lld) ***\n", SCIPnodeGetNumber(current_node), parent_num);
      printf("%d active dagcluster constraints at this node\n", num_active_cons);
   }*/

   /* Look at each of the enabled constraints in turn */
   for( i = 0; i < num_active_cons; i++ )
   {
      SCIP_CONS* cons = conss[i];
      if( SCIPconsIsEnabled(cons) )
      {
         int num_components;
         ParentSetData** components;
         SCIP_CONSDATA* consdata = SCIPconsGetData(cons);

         /* Find the strongly connected components in this constraint's data */
         SCIP_CALL( PS_splitToComponents(scip, consdata->psd, &num_components, &components) );

         /*{
            printf("Found %d components\n", num_components);
            for (j = 0; j < num_components; j++)
               printf("%d ", components[j]->n);
            printf("\n");
         }*/

         /* If there is more than one SCC, */
         /* add a new local constraint for any SCC with more than two nodes in it */
         /* and deactivate the current constraint locally */
         if( num_components > 1 )
         {
            SCIP_CONS* new_cons;
            for( j = 0; j < num_components; j++ )
            {
               if( components[j]->n > 2 )
               {
                  SCIP_CALL( DC_createCons(scip, &new_cons, "DagCluster", components[j], TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPaddConsLocal(scip, new_cons, current_node) );
                  SCIP_CALL( SCIPreleaseCons(scip, &new_cons) );
               }
            }
            SCIP_CALL( SCIPdelConsLocal(scip, cons) );
         }

         /* Clean up */
         for( j = 0; j < num_components; j++ )
            SCIP_CALL( PS_deallocateParentSetData(scip, &(components[j]), FALSE) );
         SCIPfreeMemoryArray(scip, &components);
      }
   }

   SCIPfreeMemoryArray(scip, &conss);

   return SCIP_OKAY;
}






static
SCIP_Bool shouldTryFractionalCircuitCuts(
   SCIP_CONSHDLRDATA* conshdlrdata,
   int nGen
   )
{
   if( conshdlrdata->fractionalalways )
      return TRUE;
   if( conshdlrdata->fractionalever && nGen == 0 )
      return TRUE;
   return FALSE;
}

static
SCIP_Bool shouldTrySubIPCuts(
   SCIP_CONSHDLRDATA* conshdlrdata,
   int nGen,
   int depth
   )
{
   SCIP_Bool depth_less_than_limit = (depth <= conshdlrdata->subipdepthlimit);
   SCIP_Bool depth_limit_is_off = (conshdlrdata->subipdepthlimit == -1);
   SCIP_Bool depth_is_acceptable = depth_less_than_limit || depth_limit_is_off;

   if( depth_is_acceptable )
   {
      if( conshdlrdata->subipalways )
         return TRUE;
      if( conshdlrdata->subipever && nGen == 0 )
         return TRUE;
   }
   return FALSE;
}

static
SCIP_Bool shouldTryMatroidCuts(
   SCIP_CONSHDLRDATA* conshdlrdata
   )
{
   return conshdlrdata->matroidcuts;
}

static
SCIP_RETCODE DagClusterSeparate(
   SCIP*           scip,               /**< SCIP data structure */
   SCIP_CONSDATA*  consdata,           /**< constraint data */
   SCIP_SOL*       sol,                /**< solution to be separated */
   int*            nGen,               /**< *nGen is number of cutting planes added ( even non-efficacious ones are added ) */
   int             k_lb,               /**< lowerbound on 'k' values for k-cluster searching, always positive */
   int             k_ub,               /**< upperbound on 'k' values for k-cluster searching */
   SCIP_CONSHDLR*  conshdlr,           /**< constraint handler */
   SCIP_Bool       addtopool,          /**< whether to add any found cut to the global cut pool */
   SCIP_Bool       forcecuts,           /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,   /**< to return whether an efficacious cutting plane was found */
   SCIP_Real       limits_time,             /**< limit on how long to spend sub-IP solving */
   SCIP_Real       limits_gap,              /**< limit on size of gap in sub-IP */
   SCIP_Real       limits_absgap,           /**< limit on size of the absolute gap in sub-IP */
   SCIP_Bool       incumbent_cons,          /**< whether to consider only cutting planes on which the incumbent lies */
   SolutionInfo*   solinfo,
   SCIP_Bool*      cutoff                   /**< output: pointer to store whether we detected a cutoff */
   )
{
   int num_gen;
   SCIP_Bool found_efficacious;

   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   (*nGen) = 0;
   (*found_efficacious_ptr) = FALSE;
   *cutoff = FALSE;

   if( SCIPgetDepth(scip) == 0 && SCIPgetGap(scip) < conshdlrdata->rootgapsepalimit )
      return SCIP_OKAY;

   num_gen = 0;
   found_efficacious = FALSE;
   SCIP_CALL( CC_findCuts(scip, conshdlr, consdata->psd, sol, consdata->ccs, &num_gen, forcecuts, &found_efficacious, cutoff) );
   if( *cutoff )
      return SCIP_OKAY;
   (*nGen) += num_gen;
   (*found_efficacious_ptr) |= found_efficacious;

   if( shouldTryFractionalCircuitCuts(conshdlrdata, *nGen) )
   {
      num_gen = 0;
      found_efficacious = FALSE;
      SCIP_CALL( FC_findCuts(scip, conshdlr, consdata->psd, sol, consdata->ccs, &num_gen, forcecuts, &found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
      (*nGen) += num_gen;
      (*found_efficacious_ptr) |= found_efficacious;
   }

   if( shouldTryMatroidCuts(conshdlrdata) )
   {
      num_gen = 0;
      found_efficacious = FALSE;

      SCIP_CALL( Matroid_findCuts(scip, consdata->psd, solinfo, sol, &num_gen, conshdlr, addtopool, forcecuts, &found_efficacious, limits_time, limits_gap, limits_absgap, conshdlrdata->circuit_size_lim, conshdlrdata->ground_set_size_lim, NULL, 0, NULL, 0, cutoff) );

      if( *cutoff )
         return SCIP_OKAY;
      (*nGen) += num_gen;
      (*found_efficacious_ptr) |= found_efficacious;
   }

   if( shouldTrySubIPCuts(conshdlrdata, *nGen, SCIPgetDepth(scip)) )
   {
      num_gen = 0;
      found_efficacious = FALSE;

      SCIP_CALL( IP_findCuts(scip, consdata->psd, solinfo, sol, &num_gen, k_lb, k_ub, conshdlr, addtopool, forcecuts, &found_efficacious, limits_time, limits_gap, limits_absgap, incumbent_cons, NULL, 0, NULL, 0, FALSE, conshdlrdata->matroidpaircuts, conshdlrdata->matroidpaircutslimit, conshdlrdata->matroidpairmaxcuts, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
      (*nGen) += num_gen;
      (*found_efficacious_ptr) |= found_efficacious;
   }

   return SCIP_OKAY;
}

/* Callback methods of constraint handler */
/** copy method for constraint handler plugins (called when SCIP copies plugins) */

static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyDagcluster)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( DC_includeConshdlr(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyDagcluster)
{
   SCIP_CONSDATA* sourcedata;
   const char* consname;
   ParentSetData* temppsd;

   assert(scip != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);
   assert(consname != NULL);

   /* Get source data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( PS_deepCopyParentSetData(scip, sourcescip, sourcedata->psd, &temppsd, consmap, varmap, global, valid) );

   if(*valid)
   {
      /* create the copied constraint into the targetSCIP */
      SCIPdebugMessage(" Create and capture copied constraint \n") ;
      SCIP_CALL( DC_createCons (scip, cons, consname, temppsd, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }
   /* Clean up: the cloned parent set data */
   SCIP_CALL( PS_deallocateParentSetData (scip , &temppsd, FALSE) );
   SCIPdebugMessage("In copying method : Memory freed \n");

   // {
   //    int nVars = SCIPgetNVars(scip);
   //    int n = SCIPgetNVars(sourcescip);
   //    assert( nVars == n);
   // }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitDagcluster NULL
/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreDagcluster NULL
/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolDagcluster NULL
/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#define consExitsolDagcluster NULL
/** presolving method of constraint handler */
#define consPresolDagcluster NULL
/** constraint activation notification method of constraint handler */
#define consActiveDagcluster NULL
/** constraint deactivation notification method of constraint handler */
#define consDeactiveDagcluster NULL
/** constraint enabling notification method of constraint handler */
#define consEnableDagcluster NULL
/** constraint disabling notification method of constraint handler */
#define consDisableDagcluster NULL
/** variable deletion of constraint handler */
#define consDelvarsDagcluster NULL


/* Initialisation Callbacks */
/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreDagcluster)
{

   int c;

   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;

   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !(conshdlrdata->addspc) )
      return SCIP_OKAY;

   SCIPdebugMessage("Adding set packing constraints in dagcluster constraint handler initialisation.\n");

   for( c = 0; c < nconss; ++c )
   {

      cons = conss[c];
      assert(cons != NULL);

      consdata = SCIPconsGetData(cons);

      /* now add in 'SPC' constraints */

      SCIP_CALL( SP_add_spc_constraints(scip, consdata->psd, consdata->store, conshdlrdata->propagatespc) ); 

   }

   SCIPdebugMessage("Added all set packing constraints in dagcluster constraint handler initialisation.\n");

   return SCIP_OKAY;
}
/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpDagcluster)
{
   int c;

   SCIP_CONSDATA* consdata;

   SolutionInfo solinfo;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* loop through all constraints */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);
      SCIPdebugMessage("adding initial rows for dagcluster constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* now add in cutting planes derived from low dimension convex hull */

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* use a dummy solution (final arg=TRUE) to facilitate generating C4 constraints */
      SCIP_CALL( SI_setsolinfo(scip, &solinfo, consdata->psd, NULL, TRUE, TRUE) );
      SCIP_CALL( C4_add_initconvexhull_constraints(scip, conshdlr, &solinfo, consdata->psd, consdata->store) );
      SI_freesolinfo(&solinfo, consdata->psd->n);
   }

   return SCIP_OKAY;
}


/* Transformation Callbacks */
/** Transformation method of constraint handler */
#if 0
static
SCIP_DECL_CONSTRANS(consTransDagcluster)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( createTransConsData(scip, &targetdata, sourcedata->psd) );

   SCIP_CALL(SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                            SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                            SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                            SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
                            SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)));

   return SCIP_OKAY;
}
#else
#define consTransDagcluster NULL
#endif

#define consInitDagcluster NULL
/** Transformation initialisation method of constraint handler */
/* static SCIP_DECL_CONSINIT(consInitDagcluster) { */
/*    /\* All we need to do is transform each of the PaVars to their transformed equivalent. *\/ */
/*    int i, j, k; */

/*    for (i = 0; i < nconss; i++) { */
/*       SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]); */
/*       for (j = 0; j < consdata->psd->n; j++) */
/*          for (k = 0; k < consdata->psd->nParentSets[j]; k++) */
/*             consdata->psd->PaVars[j][k] = SCIPvarGetTransVar(consdata->psd->PaVars[j][k]); */
/*    } */
/*    return SCIP_OKAY; */

/* } */


/* Finalisation Callbacks */
/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDagcluster)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}
/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteDagcluster)
{

   int i;
   int k;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPdebugMessage("deleting DAG cluster constraint <%s>.\n", SCIPconsGetName(cons));

   for( i = 0; i < (*consdata)->psd->n; ++i )
   {
      for( k = 0;  k < (*consdata)->psd->nParentSets[i]; ++k )
         SCIPfreeMemoryArray(scip, &((*consdata)->store[i][k]));
      SCIPfreeMemoryArray(scip, &((*consdata)->store[i]));
   }
   SCIPfreeMemoryArray(scip, &((*consdata)->store));

   SCIP_CALL( CC_finalise(scip, (*consdata)->ccs, (*consdata)->psd->n) );
   /* SCIP_CALL( FC_finalise(scip, (*consdata)->ccs, (*consdata)->psd->n) ); */
   SCIPfreeBlockMemory(scip, &((*consdata)->ccs));

   SCIP_CALL( PS_deallocateParentSetData(scip, &((*consdata)->psd), FALSE) );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/* CIP Related Callbacks */
/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintDagcluster)
{
   SCIP_CONSDATA* consdata = SCIPconsGetData(cons);
   PS_writeToFile(scip, file, consdata->psd);
   return SCIP_OKAY;
}
/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseDagcluster)
{
   ParentSetData* psd;
   SCIP_CALL( PS_parse(scip, (char*)str, &psd) );
   SCIP_CALL( DC_createCons(scip, cons, name, psd, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   SCIP_CALL( PS_deallocateParentSetData(scip, &psd, FALSE) );

   *success = TRUE;
   return SCIP_OKAY;
}


/* Feasibility Checking Callbacks */
/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDagcluster)
{
   if( areFeasible(scip, conss, nconss, NULL) )
   {
      *result = SCIP_FEASIBLE;
      SCIPdebugMessage("pseudo solution satifies the DAG cluster constraints.\n");
   }
   else
   {
      *result = SCIP_INFEASIBLE;
      SCIPdebugMessage("pseudo solution DOES NOT satify the DAG cluster constraints.\n");
   }
   return SCIP_OKAY;
}
/** feasibility check method of constraint handler for integral solutions **/
static
SCIP_DECL_CONSCHECK(consCheckDagcluster)
{
   if( areFeasible(scip, conss, nconss, sol) )
   {
      *result = SCIP_FEASIBLE;
      SCIPdebugMessage("primal solution satifies the DAG cluster constraints.\n");
   }
   else
   {
      *result = SCIP_INFEASIBLE;
      SCIPdebugMessage("primal solution DOES NOT satify the DAG cluster constraints.\n");
   }
   return SCIP_OKAY;
}


/* Separation Callbacks */
/** separation method of constraint handler for DAG  solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpDagcluster)
{
   int c;
   int nGen = 0;
   int totalnGen = 0;

   char mip_name[SCIP_MAXSTRLEN];

   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;

   SCIP_Bool found_efficacious = FALSE;
   SCIP_Bool found_any_efficacious = FALSE;

   SolutionInfo solinfo;
   SCIP_Bool cutoff = FALSE;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->sepalp_calls++;

   if( conshdlrdata->sepalp_calls <= conshdlrdata->writemip_maxsepalp_calls  && conshdlrdata->sepalp_calls >= conshdlrdata->writemip_minsepalp_calls )
   {
      (void) SCIPsnprintf(mip_name, SCIP_MAXSTRLEN, "gobnilp_%d.lp", conshdlrdata->sepalp_calls);
      SCIP_CALL( SCIPwriteMIP(scip, mip_name, FALSE, TRUE, FALSE) );
   }
   *result = SCIP_DIDNOTRUN;

   /* to avoid cycling */
   if( conshdlrdata->nosep_dagcluster )
   {
      SCIPdebugMessage("Skipping dag cluster separator since it failed to find any efficacious cuts when last called (next time won't be skipped) .\n");
      conshdlrdata->nosep_dagcluster = FALSE;
      return SCIP_OKAY;
   }

   for( c = 0; c < nconss; ++c )
   {


      cons = conss[c];
      assert(cons != NULL);
      SCIPdebugMessage("separating this LP solution for dag cluster constraint <%s>.\n", SCIPconsGetName(cons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintSol(scip, NULL, NULL, FALSE) );
#endif

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;

      /* store useful information on NULL = the LP solution */
      SCIP_CALL( SI_setsolinfo(scip, &solinfo, consdata->psd, NULL, TRUE, FALSE) );

      if( conshdlrdata->useincumbentcons )
      {
         SCIP_CALL( DagClusterSeparate(scip, consdata, NULL, &nGen, 1, 1, conshdlr, conshdlrdata->clustercuts_addtopool, conshdlrdata->forcecuts, &found_efficacious, conshdlrdata->conssepalpsubiptimelimit, conshdlrdata->conssepalpsubipgaplimit, conshdlrdata->conssepalpsubipabsgaplimit, TRUE, &solinfo, &cutoff) );
         if( cutoff )
         {
            *result = SCIP_CUTOFF;
            SI_freesolinfo(&solinfo, consdata->psd->n);
            return SCIP_OKAY;
         }
         totalnGen += nGen;
         found_any_efficacious = found_any_efficacious || found_efficacious;
      }

      SCIP_CALL( DagClusterSeparate(scip, consdata, NULL, &nGen, 1, 1, conshdlr, conshdlrdata->clustercuts_addtopool, conshdlrdata->forcecuts, &found_efficacious, conshdlrdata->conssepalpsubiptimelimit, conshdlrdata->conssepalpsubipgaplimit, conshdlrdata->conssepalpsubipabsgaplimit, FALSE, &solinfo, &cutoff) );
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         SI_freesolinfo(&solinfo, consdata->psd->n);
         return SCIP_OKAY;
      }

      totalnGen += nGen;
      found_any_efficacious = found_any_efficacious || found_efficacious;

      SCIP_CALL( C4_add_convexhull_constraints(scip, &solinfo, consdata->psd, consdata->store, conshdlr, NULL, &nGen, conshdlrdata->forcecuts, &found_efficacious, &cutoff) );
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         SI_freesolinfo(&solinfo, consdata->psd->n);
         return SCIP_OKAY;
      }

      totalnGen += nGen;
      found_any_efficacious = found_any_efficacious || found_efficacious;

      if( !found_efficacious )
      {
         /* look (possibly) for k-cluster constraints */

         if( SCIPgetDepth(scip) == 0 )
         {
            if( conshdlrdata->kmaxroot > 1 )
            {
               /*printf("separating in root with %d\n",conshdlrdata->kmaxroot);*/
               SCIP_CALL( DagClusterSeparate(scip, consdata, NULL, &nGen,  2, conshdlrdata->kmaxroot, conshdlr, conshdlrdata->clustercuts_addtopool, conshdlrdata->forcecuts, &found_efficacious, conshdlrdata->conssepalpsubiptimelimit, conshdlrdata->conssepalpsubipgaplimit, conshdlrdata->conssepalpsubipabsgaplimit, FALSE, &solinfo, &cutoff) );
               if( cutoff )
               {
                  *result = SCIP_CUTOFF;
                  SI_freesolinfo(&solinfo, consdata->psd->n);
                  return SCIP_OKAY;
               }

               totalnGen += nGen;
               found_any_efficacious = found_any_efficacious || found_efficacious;
            }
         }
         else
         {
            if( conshdlrdata->kmax > 1 )
            {
               /*printf("separating outside of root with %d\n",conshdlrdata->kmax);*/
               SCIP_CALL( DagClusterSeparate(scip, consdata, NULL, &nGen,  2, conshdlrdata->kmax, conshdlr, conshdlrdata->clustercuts_addtopool, conshdlrdata->forcecuts, &found_efficacious, conshdlrdata->conssepalpsubiptimelimit, conshdlrdata->conssepalpsubipgaplimit, conshdlrdata->conssepalpsubipabsgaplimit, FALSE, &solinfo, &cutoff) );
               if( cutoff )
               {
                  *result = SCIP_CUTOFF;
                  SI_freesolinfo(&solinfo, consdata->psd->n);
                  return SCIP_OKAY;
               }

               totalnGen += nGen;
               found_any_efficacious = found_any_efficacious || found_efficacious;
            }
         }
      }

      /* else if ( SCIPgetDepth(scip) == 0 ) */
      /*     SCIP_CALL(  SCIPwriteMIP(scip,"foo",FALSE,TRUE)  ); */

      SI_freesolinfo(&solinfo, consdata->psd->n);
      
   }

   if( found_any_efficacious )
      *result = SCIP_SEPARATED;
   else
   {
      conshdlrdata->nosep_dagcluster = TRUE;
      SCIPdebugMessage("Could not find efficacious cut in separation method.\n");
   }

   SCIPdebugMessage("separated %d cuts in separation method.\n", totalnGen);

   return SCIP_OKAY;
}
/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolDagcluster)
{
   int c;
   int totalnGen = 0;
   int nGen = 0;

   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;

   SCIP_Bool found_efficacious = FALSE;
   SCIP_Bool found_any_efficacious = FALSE;

   SolutionInfo solinfo;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;

   /* to avoid cycling */
   if( conshdlrdata->nosepasol_dagcluster )
   {
      SCIPdebugMessage("Skipping dag cluster separator for arbitrary primal solutions since it failed to find any efficacious cuts when last called (next time won't be skipped) .\n");
      conshdlrdata->nosepasol_dagcluster = FALSE;
      return SCIP_OKAY;
   }

   /* loop through all constraints */
   for( c = 0; c < nconss; ++c )
   {

      cons = conss[c];
      assert(cons != NULL);
      SCIPdebugMessage("separating this solution for an arbitrary primal solution for the DAG cluster constraint <%s>.\n", SCIPconsGetName(cons));
#ifdef SCIP_DEBUG
      SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      if( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;

      /* store useful information on primal solution */
      SCIP_CALL( SI_setsolinfo(scip, &solinfo, consdata->psd, sol, TRUE, FALSE) );

      SCIP_CALL( DagClusterSeparate(scip, consdata, sol, &nGen, 1, 1, conshdlr, conshdlrdata->clustercuts_addtopool, conshdlrdata->forcecuts, &found_efficacious, conshdlrdata->conssepasolsubiptimelimit, conshdlrdata->conssepasolsubipgaplimit, conshdlrdata->conssepasolsubipabsgaplimit, FALSE, &solinfo, &cutoff) );
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         SI_freesolinfo(&solinfo, consdata->psd->n);
         return SCIP_OKAY;
      }

      totalnGen += nGen;
      found_any_efficacious = found_any_efficacious || found_efficacious;

      if( conshdlrdata->sepasol_useconvex4 )
      {
         SCIP_CALL( C4_add_convexhull_constraints(scip, &solinfo, consdata->psd, consdata->store, conshdlr, sol, &nGen, conshdlrdata->forcecuts, &found_efficacious, &cutoff) );
         if( cutoff )
         {
            *result = SCIP_CUTOFF;
            SI_freesolinfo(&solinfo, consdata->psd->n);
            return SCIP_OKAY;
         }

         totalnGen += nGen;
         found_any_efficacious = found_any_efficacious || found_efficacious;
      }

      SI_freesolinfo(&solinfo, consdata->psd->n);

   }
   if( found_any_efficacious )
      *result = SCIP_SEPARATED;
   else
   {
      conshdlrdata->nosepasol_dagcluster = TRUE;
      SCIPdebugMessage("Could not find efficacious cut in separation method for arbitrary primal solutions.\n");
   }
   SCIPdebugMessage("separated %d cuts in arbitrary primal solution method.\n", totalnGen);

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDagcluster)
{
   int c;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;

   SCIP_Bool ancestralbranching;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/ancestralbranching", &ancestralbranching) );
   
   if( ancestralbranching )
   {
      for( c = 0; c < nconss; ++c )
      {
         int n;
         ParentSetData* psd;
         int i;
         int k;
         int l;
         int* selected;
         SCIP_Bool* ancestral;
         SCIP_Bool made_progress;
         SCIP_Bool all_ancestral;

         cons = conss[c];
         assert(cons != NULL);
         consdata = SCIPconsGetData(cons);
         assert(consdata != NULL);
         psd = consdata->psd;
         n = psd->n;

         SCIP_CALL( SCIPallocBufferArray(scip, &ancestral, n) );
         SCIP_CALL( SCIPallocBufferArray(scip, &selected, n) );

         /* find any selected parent sets and mark vertices with selected empty parent sets as ancestral */
         made_progress = FALSE;
         for( i = 0; i < n; ++i )
         {
            ancestral[i] = FALSE;
            selected[i] = -1;
            for( k = 0 ; k < psd->nParentSets[i]; ++k )
            {
               if( SCIPvarGetLbLocal(psd->PaVars[i][k]) > 0.5 )
               {
                  if( psd->nParents[i][k] == 0 )
                  {
                     ancestral[i] = TRUE;
                     made_progress = TRUE;
                  }
                  selected[i] = k;
                  break;
               }
            }
         }
      
         while( made_progress )
         {
            made_progress = FALSE;
            for( i = 0; i < n; ++i )
            {
            
               k = selected[i];
               if( k == -1 || ancestral[i] )
                  continue;

               all_ancestral = TRUE;
               for(l = 0; l < psd->nParents[i][k]; ++l )
               {
                  if( !ancestral[psd->ParentSets[i][k][l]] )
                  {
                     all_ancestral = FALSE;
                     break;
                  }
               }

               if( all_ancestral )
               {
                  ancestral[i] = TRUE;
                  made_progress = TRUE;
               }
            }
         }
   
         /* update branching priorities of (unfixed) family variables */
         for( i = 0; i < n; ++i )
         {
            if( selected[i] != -1 )
               continue;
         
            for( k = 0 ; k < psd->nParentSets[i]; ++k )
            {
               all_ancestral = TRUE;
               for(l = 0; l < psd->nParents[i][k]; ++l )
               {
                  if( !ancestral[psd->ParentSets[i][k][l]] )
                  {
                     all_ancestral = FALSE;
                     break;
                  }
               }
               SCIP_CALL( SCIPchgVarBranchPriority(scip, psd->PaVars[i][k], all_ancestral ? 100 : -100 ) );
            }
         }
         SCIPfreeBufferArray(scip, &selected);
         SCIPfreeBufferArray(scip, &ancestral);
      }
   }

   if( areFeasible(scip, conss, nconss, NULL) )
   {
      *result = SCIP_FEASIBLE;
      SCIPdebugMessage("LP solution satifies the DAG cluster constraints.\n");
   }
   else
   {
      *result = SCIP_INFEASIBLE;
      SCIPdebugMessage("LP solution DOES NOT satify the DAG cluster constraints.\n");
   }
   return SCIP_OKAY;
}

/* Other Callbacks */
/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropDagcluster)
{
   int c;
   int nGen = 0;

   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   int i;
   int j;
   int kk;

   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   int n;                           /* number of nodes involved in a given dagcluster constraint */
   SCIP_Bool** ancestor;            /* ancestor[i][j] = TRUE iff j is an ancestor of i in current node */
   SCIP_Bool** pot_ancestor = NULL;        /* pot_ancestor[i][j] = TRUE iff j is a potential ancestor of i in current node */
   SCIP_Bool* orphan;               /* orphan[i] = TRUE iff there is no j for which ancestor[i][j] is TRUE */
   SCIP_Bool* allowed;              /* temporary storage */


   SCIP_VAR* arrow_i_j;
   SCIP_VAR* arrow_j_i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss != NULL);
   assert(result != NULL);


   *result = SCIP_DIDNOTRUN;

   /* loop through all constraints */
   for( c = 0; c < nconss; ++c )
   {

      cons = conss[c];
      assert(cons != NULL);
      SCIPdebugMessage("propagating dagcluster constraint <%s>.\n", SCIPconsGetName(cons));

      if( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
      consdata = SCIPconsGetData(cons);

      assert(consdata != NULL);
      assert(consdata->psd != NULL);
      assert(consdata->psd->arrow != NULL);

      n = consdata->psd->n;

      SCIP_CALL( SCIPallocBufferArray(scip, &orphan, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &allowed, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ancestor, n) );
      for( i = 0; i < n; ++i )
         SCIP_CALL( SCIPallocBufferArray(scip, &(ancestor[i]), n) );

      for( i = 0; i < n; ++i )
      {
         orphan[i] = TRUE;
         for( j = 0; j < n; ++j )
         {
            if( i == j)
               continue;

            arrow_i_j = get_arrow(consdata->psd,i,j);
            if( arrow_i_j != NULL && SCIPvarGetLbLocal(arrow_i_j) > 0.5 )
            {
               ancestor[i][j] = TRUE;
               orphan[i] = FALSE;
            }
            else
               ancestor[i][j] = FALSE;
         }
      }

      /* compute transitive closure */

      /* Cormen et al:
         ancestor[i][j][k] = TRUE iff there exists
         a path from i to j and any intermediate vertices
         are in {0, ..., k-1}
         But drop k index and just update ancestor[i][j]
      */

      for( kk = 0; kk < n; ++kk )
      {
         if( orphan[kk] ) 
            continue; 
         
         for( i = 0; i < n; ++i )
         {
            if( i == kk || orphan[i] ) 
               continue; 
            
            for( j = 0; j < n; ++j )
            {
               if( j == kk || j == i )
                  continue;
               if( ancestor[i][kk] && ancestor[kk][j] )
                  ancestor[i][j] = TRUE;
            }
         }
      }
      
      /* propagations using ancestor */
      for( i = 0 ; i < n ; ++i )
      {
         if( orphan[i] )
            continue;

         for( j = 0; j < n; ++j )
            if( i != j && ancestor[i][j] )
            {
               arrow_j_i = get_arrow(consdata->psd,j,i);
               if( arrow_j_i != NULL )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip,arrow_j_i,0,TRUE,&infeasible,&tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible.\n");
                     *result = SCIP_CUTOFF;
                     goto TIDYUP;
                  }
                  if( tightened )
                  {
                     SCIPdebugMessage("Ruling out: %s (to prevent a cycle) \n", SCIPvarGetName(arrow_j_i));
                     ++nGen;
                  }
               }
            }
      }


   TIDYUP:
      for( i = n-1; i > -1; --i )
         SCIPfreeBufferArray(scip, &(ancestor[i]));
      SCIPfreeBufferArray(scip, &ancestor);
      SCIPfreeBufferArray(scip, &allowed);
      SCIPfreeBufferArray(scip, &orphan);
      if( pot_ancestor != NULL )
      {
         for( i = 0; i < n; ++i )
            SCIPfreeMemoryArray(scip, &(pot_ancestor[i]));
         SCIPfreeMemoryArray(scip, &pot_ancestor);
      }

      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if( nGen > 0 )
      *result = SCIP_REDUCEDDOM;
   SCIPdebugMessage("propagated %d domains.\n", nGen);

   return SCIP_OKAY;
}
/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropDagcluster)
{
   /*SCIP_CONSDATA* consdata;*/

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(infervar != NULL);
   assert(bdchgidx != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Propagation resolution of constraint <%s>.\n", SCIPconsGetName(cons));
   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDagcluster)
{

   SCIP_CONSDATA* consdata;
   int i;
   int j;
   SCIP_VAR* arrow_i_j;
      

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->psd->arrow != NULL);

   /* only lock arrow variables.
      an acyclicity constraint can only be violated by adding, not removing, an arrow
      hence "nlocksneg, nlockspos" below
   */
   
   for( i = 0; i < consdata->psd->n; ++i )
      for( j = 0; j < consdata->psd->n; ++j )
         if( i != j )
         {
            arrow_i_j = get_arrow(consdata->psd, i, j);
            if( arrow_i_j != NULL )
               SCIP_CALL( SCIPaddVarLocks(scip, arrow_i_j, nlocksneg, nlockspos) );
         }
   
   return SCIP_OKAY;
}




/** creates the handler for dagcluster constraints and includes it in SCIP */
SCIP_RETCODE DC_includeConshdlr(SCIP* scip)
{

   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* include constraint handler */

   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
      CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
      consEnfolpDagcluster, consEnfopsDagcluster, consCheckDagcluster, consLockDagcluster,
      conshdlrdata) );
   assert(conshdlr != NULL);
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreDagcluster) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpDagcluster) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyDagcluster, consCopyDagcluster) );
   /*SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransDagcluster) );*/
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeDagcluster) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteDagcluster) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintDagcluster) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseDagcluster) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpDagcluster, consSepasolDagcluster,
         CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropDagcluster, CONSHDLR_PROPFREQ,
         CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropDagcluster) );

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/"CONSHDLR_NAME"/kmax",
                             "maximum k to try for k-cluster cutting planes",
                             &conshdlrdata->kmax, FALSE, DEFAULT_KMAX, 1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/"CONSHDLR_NAME"/kmaxroot",
                             "maximum k to try for k-cluster cutting planes in the root",
                             &conshdlrdata->kmaxroot, FALSE, DEFAULT_KMAX_ROOT, 1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/clustercuts_addtopool",
                              "whether to add cluster cuts to the global cut pool",
                              &conshdlrdata->clustercuts_addtopool, FALSE, DEFAULT_CLUSTERCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/propagatespc",
                              "whether to propagate added set packing constraints",
                              &conshdlrdata->propagatespc, FALSE, DEFAULT_PROPAGATESPC, NULL, NULL));

   conshdlrdata->nosep_dagcluster = FALSE;
   conshdlrdata->nosepasol_dagcluster = FALSE;
   conshdlrdata->sepalp_calls = 0;

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/"CONSHDLR_NAME"/writemip_minsepalp_calls",
                             "only write out pre-separation MIP relaxation after reaching this number of calls to separator",
                             &conshdlrdata->writemip_minsepalp_calls, FALSE, 1, 1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/"CONSHDLR_NAME"/writemip_maxsepalp_calls",
                             "only write out pre-separation MIP relaxation if not beyond this number of calls to separator",
                             &conshdlrdata->writemip_maxsepalp_calls, FALSE, 0, 0, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/"CONSHDLR_NAME"/subipdepthlimit",
                             "do not use the subIP to search for cluster cuts if deeper than this (-1 for no limit)",
                             &conshdlrdata->subipdepthlimit, FALSE, DEFAULT_SUBIPDEPTHLIMIT, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/useincumbentcons",
                              "whether to look for cluster cuts which are tight for the current incumbent",
                              &conshdlrdata->useincumbentcons, FALSE, DEFAULT_USEINCUMBENTCONS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/subipalways",
                              "whether to always use the subIP to look for cluster cuts",
                              &conshdlrdata->subipalways, FALSE, DEFAULT_SUBIPALWAYS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/subipever",
                              "whether to ever use the subIP to look for cluster cuts",
                              &conshdlrdata->subipever, FALSE, DEFAULT_SUBIPEVER, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/fractionalalways",
                              "whether to always use fractional cycles to look for cluster cuts",
                              &conshdlrdata->fractionalalways, FALSE, DEFAULT_FRACTIONALALWAYS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/fractionalever",
                              "whether to ever use fractional cycles to look for cluster cuts",
                              &conshdlrdata->fractionalever, FALSE, DEFAULT_FRACTIONALEVER, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/sepasol_useconvex4",
                              "whether to search for convex4 cuts when separating an arbitrary primal solution",
                              &conshdlrdata->sepasol_useconvex4, FALSE, DEFAULT_SEPASOL_USECONVEX4, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/conssepalpsubiptimelimit",
                              "time limit on sub IP for finding cluster cuts in LP separation",
                              &conshdlrdata->conssepalpsubiptimelimit, FALSE, DEFAULT_CONSSEPALPSUBIPTIMELIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/conssepasolsubiptimelimit",
                              "time limit on sub IP for finding cluster cuts in arbitrary primal solution separation",
                              &conshdlrdata->conssepasolsubiptimelimit, FALSE, DEFAULT_CONSSEPASOLSUBIPTIMELIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/conssepalpsubipgaplimit",
                              "gap limit on sub IP for finding cluster cuts in LP separation",
                              &conshdlrdata->conssepalpsubipgaplimit, FALSE, DEFAULT_CONSSEPALPSUBIPGAPLIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/conssepasolsubipgaplimit",
                              "gap limit on sub IP for finding cluster cuts in arbitrary primal solution separation",
                              &conshdlrdata->conssepasolsubipgaplimit, FALSE, DEFAULT_CONSSEPASOLSUBIPGAPLIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/conssepalpsubipabsgaplimit",
                              "absgap limit on sub IP for finding cluster cuts in LP separation",
                              &conshdlrdata->conssepalpsubipabsgaplimit, FALSE, DEFAULT_CONSSEPALPSUBIPABSGAPLIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/conssepasolsubipabsgaplimit",
                              "absgap limit on sub IP for finding cluster cuts in arbitrary primal solution separation",
         &conshdlrdata->conssepasolsubipabsgaplimit, FALSE, DEFAULT_CONSSEPASOLSUBIPABSGAPLIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/forcecuts",
                              "whether to force all cuts to be added",
                              &conshdlrdata->forcecuts, TRUE, DEFAULT_FORCECUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/"CONSHDLR_NAME"/extraprops",
                              "whether to do extra (slow) propagations",
                              &conshdlrdata->extraprops, TRUE, DEFAULT_EXTRAPROPS, NULL, NULL));

   SCIP_CALL(SCIPaddRealParam(scip,
                              "constraints/"CONSHDLR_NAME"/rootgapsepalimit",
                              "no cuts generated in root if gap is below this value",
                              &conshdlrdata->rootgapsepalimit, FALSE, DEFAULT_ROOTGAPSEPALIMIT, 0, SCIP_REAL_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/matroidpaircuts",
         "whether to search for matroid pair cuts",
         &conshdlrdata->matroidpaircuts, TRUE, DEFAULT_MATROIDPAIRCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/matroidpaircutslimit",
         "limit on size of clusters to allow searching for matroid pair cuts",
         &conshdlrdata->matroidpaircutslimit, FALSE, DEFAULT_MATROIDPAIRCUTSLIMIT, 2, INT_MAX, NULL, NULL));

      SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/matroidpairmaxcuts",
         "limit on number of matroid pair cuts which can be added for a give cluster",
         &conshdlrdata->matroidpairmaxcuts, FALSE, DEFAULT_MATROIDPAIRMAXCUTS, 0, INT_MAX, NULL, NULL));

   
   SCIP_CALL(SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/matroidcuts",
         "whether to search for matroid cuts",
         &conshdlrdata->matroidcuts, TRUE, DEFAULT_MATROIDCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/circuit_size_lim",
         "limit on size of circuits when searching for matroid cuts",
         &conshdlrdata->circuit_size_lim, FALSE, DEFAULT_CIRCUIT_SIZE_LIM, 2, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/ground_set_size_lim",
         "limit on size of ground set when searching for matroid cuts",
         &conshdlrdata->ground_set_size_lim, FALSE, DEFAULT_GROUND_SET_SIZE_LIM, 2, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/addspc",
         "whether to add 'set packing' constraints",
         &conshdlrdata->addspc, TRUE, DEFAULT_ADDSPC, NULL, NULL));

   
   SCIP_CALL( CC_addParams(scip) );
   SCIP_CALL( FC_addParams(scip) );
   SCIP_CALL( C4_addParams(scip) );

   return SCIP_OKAY;
}

/** creates and captures a dagcluster constraint */
SCIP_RETCODE DC_createCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                     *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
)
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the dagcluster constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("dagcluster constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* Initialise constraint data */
   SCIP_CALL( createConsData(scip, &consdata, psd) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

