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
 *  Implements a method for finding and adding cluster cuts.
 *
 *  The technique involves finding all elementary circuits in the current
 *  (possibly partially fractional) solution of the linear relaxation, and
 *  adding cuts to rule out each of these.
 */

#include "circuit_cuts.h"
#include "utils.h"

#define DEFAULT_MAX_CYCLES 1000
#define DEFAULT_ADD_FRACTIONAL_CUTS TRUE
#define DEFAULT_ADD_CLUSTER_CUTS TRUE
#define DEFAULT_ADD_CYCLE_CUTS FALSE
#define DEFAULT_ADD_CLUSTER_CUTS_TO_POOL FALSE
#define DEFAULT_ADD_CYCLE_CUTS_TO_POOL FALSE
#define DEFAULT_MAX_CYCLE_LENGTH 6

/** Adds parameters to those recognised by SCIP.
 *  @param scip The SCIP instance to add parameters to.
 *  @return SCIP_OKAY if the operation succeded or an error othewrwise.
 */
SCIP_RETCODE CC_addParams(
   SCIP* scip                  /**< SCIP data structure */
   )
{

   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/circuit_cuts/add_cluster_cuts",
         "whether cluster cuts should be added based on cycles in the current solution",
         DEFAULT_ADD_CLUSTER_CUTS
         ));

   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/circuit_cuts/add_cluster_cuts_to_pool",
         "whether cluster cuts should be added to the pool based on cycles in the current solution",
         DEFAULT_ADD_CLUSTER_CUTS_TO_POOL
         ));

   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/circuit_cuts/add_cycle_cuts",
         "whether cycle cuts should be added based on cycles in the current solution",
         DEFAULT_ADD_CYCLE_CUTS
         ));

   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/circuit_cuts/add_cycle_cuts_to_pool",
         "whether cycle cuts should be added to the pool based on cycles in the current solution",
         DEFAULT_ADD_CYCLE_CUTS_TO_POOL
         ));

   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/circuit_cuts/use_fractional",
         "whether cuts should be added based on fractional cycles in the current solution",
         DEFAULT_ADD_FRACTIONAL_CUTS
         ));

   SCIP_CALL(SCIPaddIntParam(scip,
         "gobnilp/circuit_cuts/max_length",
         "the maximum length of cycle to search for",
         NULL, FALSE, DEFAULT_MAX_CYCLE_LENGTH, 0, INT_MAX, NULL, NULL));

   return SCIP_OKAY;
}

/** Sets up data structures once, to avoid creating them everytime the separation routine is called.
 *
 *  @param scip The SCIP instance the separation is to take place in.
 *  @param data The problem data that the separation will occur on.
 *
 *  @return SCIP_OKAY if the initialisation was successful or an appropriate error code otherwise.
 */
SCIP_RETCODE CC_initialise(
   SCIP* scip, 
   CircuitCutsStorage* ccs, 
   ParentSetData* data
   )
{
   int i;
   int num_vars;

   assert( scip != NULL );
   assert( ccs != NULL );
   assert( data != NULL );
   
   /* Read the parameters */
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/circuit_cuts/add_cluster_cuts",          &(ccs->add_cluster_cuts)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/circuit_cuts/add_cluster_cuts_to_pool",  &(ccs->add_cluster_cuts_to_pool)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/circuit_cuts/add_cycle_cuts",            &(ccs->add_cycle_cuts)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/circuit_cuts/add_cycle_cuts_to_pool",    &(ccs->add_cycle_cuts_to_pool)) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/circuit_cuts/use_fractional",            &(ccs->add_fractional_cuts)) );
   SCIP_CALL( SCIPgetIntParam(scip,  "gobnilp/circuit_cuts/max_length",                &(ccs->max_cycle_length)) );

   /* this one hard coded for some reason */
   ccs->max_cycles = DEFAULT_MAX_CYCLES;
   
   /* Create adjacency matrix */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->adj_matrix), data->n) );
   for( i = 0; i < data->n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->adj_matrix[i]), data->n) );

   /* Create a stack for storing nodes on */
   ccs->s = StackCreate(data->n);

   /* Create the structure containing the found SCC */
   ccs->components = VectorListCreate(data->n);

   /* Create the structure containing the found cycles */
   ccs->cycles = VectorListCreate(ccs->max_cycles);

   /* Create the data structures for the SCC algorithm */
   ccs->index_array = VectorCreate(data->n);
   ccs->lowlink = VectorCreate(data->n);

   /* Create the b array from the cycle finding algorithm */
   ccs->b = VectorListCreate(data->n);
   for( i = 0; i < data->n; i++ )
      VectorListAppend(ccs->b, VectorCreate(data->n));

   /* Create arrays for the variables to include in the cuts */
   num_vars = 0;
   for( i = 0; i < data->n; i++ )
      num_vars += data->nParentSets[i];
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->included_cluster), num_vars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->excluded_cluster), num_vars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->included_cycle), num_vars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->excluded_cycle), num_vars) );

   /* Create space for the blocked list */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(ccs->blocked), data->n) );

   return SCIP_OKAY;
}
/** Frees all memory allocated by @c CC_initialise.
 *
 *  @return SCIP_OKAY if the finalisation was successful or an appropriate error code otherwise.
 */
SCIP_RETCODE CC_finalise(
   SCIP* scip,                   /**< SCIP data structure */
   CircuitCutsStorage* ccs,      /**< Circuit storage to free */
   int n                         /**< Number of DAG nodes */
   )
{
   int i;

   /* Free adjacency matrix */

   for( i = 0; i < n; i++ )
      SCIPfreeMemoryArray(scip, &(ccs->adj_matrix[i]));
   SCIPfreeMemoryArray(scip, &(ccs->adj_matrix));

   /* Delete the node stack */
   StackDelete(&(ccs->s));

   /* Free the blocked list */
   SCIPfreeMemoryArray(scip, &(ccs->blocked));

   /* Free the list of SCC */
   VectorListDelete(&(ccs->components));

   /* Free the list of found cycles */
   VectorListDelete(&(ccs->cycles));

   /* Free the b array */
   VectorListDelete(&(ccs->b));

   /* Free the arrays from the SCC finding algorithm */
   VectorDelete(&(ccs->index_array));
   VectorDelete(&(ccs->lowlink));

   /* Free the cut creating arrays */
   SCIPfreeMemoryArray(scip, &(ccs->included_cluster));
   SCIPfreeMemoryArray(scip, &(ccs->excluded_cluster));
   SCIPfreeMemoryArray(scip, &(ccs->included_cycle));
   SCIPfreeMemoryArray(scip, &(ccs->excluded_cycle));

   return SCIP_OKAY;
}

/** Return whether a value is strictly less than 1
 *  Uses SCIPepsilon to guard against rounding errors
 *  @return TRUE is the input value is strictly less than 1, else FALSE
 */
 
static 
SCIP_Bool isLessThanOne(
   SCIP* scip,              /**< SCIP data structure */ 
   float value              /**< value to test */
   )
{
   return value < (1 - SCIPepsilon(scip));
}


/** Finds Strongly Connected Components in the current solution.
 *
 *  @param current_node The node to consider next.
 *  @param current_index The current depth index.
 *  @return SCIP_OKAY if the procedure worked correctly, or an error otherwsie.
 */
static 
SCIP_RETCODE findStronglyConnectedComponents(
   int current_node, 
   int* current_index, 
   int n, 
   Vector* index_array, 
   Vector* lowlink, 
   Stack* s, 
   SCIP_Real** adj_matrix, 
   VectorList* components
   )
{
   int i;

   index_array->items[current_node] = (*current_index);
   lowlink->items[current_node] = (*current_index);
   (*current_index) = (*current_index) + 1;
   StackPush(s, current_node);

   for( i = 0; i < n; i++ )
      if( adj_matrix[current_node][i] != 0 )
      {
         if( index_array->items[i] == -1 )
         {
            findStronglyConnectedComponents(i, current_index, n, index_array, lowlink, s, adj_matrix, components);
            if( lowlink->items[i] < lowlink->items[current_node] )
               lowlink->items[current_node] = lowlink->items[i];
         }
         else
         {
            if( StackContains(s, i) )
               if( index_array->items[i] < lowlink->items[current_node] )
                  lowlink->items[current_node] = index_array->items[i];
         }
      }

   if( index_array->items[current_node] == lowlink->items[current_node] )
   {
      Vector* new_component = VectorCreate(StackSize(s));
      i = StackPop(s);
      VectorAppend(new_component, i);
      while( i != current_node )
      {
         i = StackPop(s);
         VectorAppend(new_component, i);
      }
      if( new_component->size == 1 )
         VectorDelete(&new_component);
      else
         VectorListAppend(components, new_component);
   }

   return SCIP_OKAY;
}
/** Unblocks a node during the cycle finding algorithm.
 *
 *  @param u The node to unblock.
 *  @param nodes The nodes in the current SCC.
 *  @return SCIP_OKAY if the unblocking succeeded.
 */
static 
SCIP_RETCODE unblock(
   int u, 
   Vector* nodes, 
   VectorList* b, 
   SCIP_Bool* blocked
   )
{
   int i;
   Vector* b_u = b->items[nodes->items[u]];
   blocked[nodes->items[u]] = FALSE;
   for( i = 0; i < b_u->size; i++ )
   {
      if( blocked[nodes->items[b_u->items[i]]] )
         unblock(b_u->items[i], nodes, b, blocked);
   }
   VectorClear(b_u);
   return SCIP_OKAY;
}

/** Find the elementary cycles of a SCC.
 *
 *  @param scip SCIP data structure 
 *  @param nodes The nodes in the SCC.
 *  @param current_index The index of the current node.
 *  @param start_index The index of the start node for this SCC.
 *  @return True if a cycle was found.
 */
static SCIP_Bool circuit(
   SCIP* scip, 
   Vector* nodes, 
   int current_index, 
   int start_index, 
   VectorList* b, 
   SCIP_Bool* blocked, 
   Stack* s, 
   SCIP_Real** adj_matrix, 
   VectorList* cycles, 
   SCIP_Real weight,
   int max_cycles,
   int max_cycle_length
   )
{
   int i;
   int j;
   int current_value = nodes->items[current_index];
   SCIP_Bool return_value = FALSE;

   if( cycles->size == max_cycles )
   {
      /* Enough cuts already found */
      return FALSE;
   }
   else if( StackSize(s) == max_cycle_length )
   {
      /* Maximum cycle length has been reached */
      /* Returns TRUE to prevent algorithm inferring that no cycle from start_index can go through current_index */
      return TRUE;
   }
   else if( !isLessThanOne(scip, weight) )
   {
      /* "Weight" of path has become too high to be useful cut */
      /* Returns TRUE to prevent algorithm inferring that no cycle from start_index can go through current_index */
      return TRUE;
   }
   else
   {
      /* Normal iteration of algorithm */
      StackPush(s, current_value);
      blocked[current_value] = TRUE;

      for( i = start_index; i < nodes->size; i++ )
      {
         SCIP_Real edge_weight = adj_matrix[current_value][nodes->items[i]];
         if( edge_weight != 0 )
         {
            if( i == start_index && cycles->size < max_cycles )
            {
               if( isLessThanOne(scip, weight + (1 - edge_weight)) )
               {
                  Vector* cycle = VectorCreate(StackSize(s) + 1);
                  for( j = 0; j < StackSize(s); j++ )
                     VectorAppend(cycle, s->items[j]);
                  VectorAppend(cycle, nodes->items[start_index]);
                  VectorListAppend(cycles, cycle);
               }
               return_value = TRUE;
            }
            else if( !blocked[nodes->items[i]] )
               if( circuit(scip, nodes, i, start_index, b, blocked, s, adj_matrix, cycles, weight + (1 - edge_weight),
                     max_cycles, max_cycle_length) )
                  return_value = TRUE;
         }
      }

      if( return_value )
         unblock(current_index, nodes, b, blocked);
      else
         for( i = start_index; i < nodes->size; i++ )
            if( adj_matrix[current_value][nodes->items[i]] != 0 )
               if( !VectorContains(b->items[nodes->items[i]], current_index) )
                  VectorAppend(b->items[nodes->items[i]], current_index);

      StackPop(s);
      return return_value;
   }
}
/** Whether a parent set should be included in the current cut.
 *  @param cycle The nodes in the current cycle.
 *  @param previous_node The previous node in the cycle.
 *  @param child The current node in the cycle.
 *  @param parent_set The index number of the current parent set.
 *  @param in_cluster Whether the parent set should be included in the current cluster cut.
 *  @param in_cycle Whether the parent set should be included in the current cycle cut.
 *  @return SCIP_OKAY if the operation succeded or an error othewrwise.
 */
static 
SCIP_RETCODE shouldIncludeInCuts(
   Vector* cycle, 
   ParentSetData* psd, 
   int previous_node, 
   int child, 
   int parent_set, 
   SCIP_Bool* in_cluster, 
   SCIP_Bool* in_cycle
   )
{
   int i;
   (*in_cycle) = FALSE;
   (*in_cluster) = FALSE;
   for( i = 0; i < psd->nParents[child][parent_set]; i++ )
   {
      if( previous_node ==  psd->ParentSets[child][parent_set][i] )
      {
         (*in_cycle) = TRUE;
         (*in_cluster) = TRUE;
         return SCIP_OKAY;
      }
      else if( VectorContains(cycle, psd->ParentSets[child][parent_set][i]) )
      {
         (*in_cluster) = TRUE;
      }
   }
   return SCIP_OKAY;
}

/** Adds cuts to rule out each of the found cycles.
 *  @param scip The SCIP instance.
 *  @param conshdlr The constraint handler to add cuts for.
 *  @param sol The current solution.
 *  @param forcecuts Whether to force cuts to be added
 *  @param found_efficacious_ptr Whether at least one efficacious cut was added.
 *  @return SCIP_OKAY if the operation succeded or an error othewrwise.
 */

static SCIP_RETCODE addCuts(
   SCIP* scip, 
   SCIP_CONSHDLR* conshdlr, 
   ParentSetData* psd, 
   SCIP_SOL* sol, 
   SCIP_VAR** included_cluster, 
   SCIP_VAR** excluded_cluster, 
   SCIP_VAR** included_cycle, 
   SCIP_VAR** excluded_cycle, 
   VectorList* cycles, 
   SCIP_Bool forcecuts,
   SCIP_Bool* found_efficacious_ptr,
   SCIP_Bool* cutoff,
   CircuitCutsStorage* ccs
   )
{
   int i;
   int j;
   int k;
   SCIP_ROW* cluster_cut;
   SCIP_ROW* cycle_cut;

   for( i = 0; i < cycles->size; i++ )
   {
      Vector* this_cycle = cycles->items[i];
      int cycle_rhs = this_cycle->size - 2;
      int cluster_rhs = this_cycle->size - 2;

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cluster_cut, conshdlr, "clustercut", -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cycle_cut, conshdlr, "cyclecut", -SCIPinfinity(scip), SCIPinfinity(scip), FALSE, FALSE, TRUE) );

      for( j = 1; j < this_cycle->size; j++ )
      {
         int this_node = this_cycle->items[j];
         int num_included_cluster = 0, num_excluded_cluster = 0;
         int num_included_cycle = 0,   num_excluded_cycle = 0;

         for( k = 0; k < psd->nParentSets[this_node]; k++ )
         {
            SCIP_Bool in_cluster, in_cycle;
            SCIP_VAR* this_var = psd->PaVars[this_node][k];
            SCIP_CALL( shouldIncludeInCuts(this_cycle, psd, this_cycle->items[j - 1], this_node, k, &in_cluster, &in_cycle) );
            if( in_cluster )
               included_cluster[num_included_cluster++] = this_var;
            else
               excluded_cluster[num_excluded_cluster++] = this_var;
            if( in_cycle )
               included_cycle[num_included_cycle++] = this_var;
            else
               excluded_cycle[num_excluded_cycle++] = this_var;
         }

         if( num_excluded_cluster < num_included_cluster )
         {
            SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cluster_cut, num_excluded_cluster, excluded_cluster, -1.0) );
            cluster_rhs -= 1;
         }
         else
         {
            SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cluster_cut, num_included_cluster, included_cluster, 1.0) );
         }
         if( num_excluded_cycle < num_included_cycle )
         {
            SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cycle_cut, num_excluded_cycle, excluded_cycle, -1.0) );
            cycle_rhs -= 1;
         }
         else
         {
            SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cycle_cut, num_included_cycle, included_cycle, 1.0) );
         }

      }

      SCIP_CALL( SCIPchgRowRhs(scip, cluster_cut, cluster_rhs) );
      SCIP_CALL( SCIPchgRowRhs(scip, cycle_cut, cycle_rhs) );

      if( ccs->add_cluster_cuts )
      {

         SCIP_CALL( SCIPaddRow(scip, cluster_cut, forcecuts, cutoff) );
         if( *cutoff )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &cluster_cut) );
            SCIP_CALL( SCIPreleaseRow(scip, &cycle_cut) );
            return SCIP_OKAY;
         }
         if( SCIPisCutEfficacious(scip, sol, cluster_cut) )
            *found_efficacious_ptr = TRUE;
      }
      if( ccs->add_cycle_cuts )
      {

         SCIP_CALL( SCIPaddRow(scip, cycle_cut, forcecuts, cutoff) );
         if( *cutoff )
         {
            SCIP_CALL( SCIPreleaseRow(scip, &cluster_cut) );
            SCIP_CALL( SCIPreleaseRow(scip, &cycle_cut) );
            return SCIP_OKAY;
         }
      }
      if( ccs->add_cluster_cuts_to_pool )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cluster_cut) );
      }
      if( ccs->add_cycle_cuts_to_pool )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, cycle_cut) );
      }


      /*VectorPrint(this_cycle); */
      /*/*SCIP_CALL(  SCIPprintRow(scip, cluster_cut, NULL)  ); */
      /*printf("Cluster-cut: act=%f, rhs=%f, vars=%d, norm=%f, feas=%f, eff=%f\n", */
      /*SCIPgetRowLPActivity(scip, cluster_cut), SCIProwGetRhs(cluster_cut), SCIProwGetNNonz(cluster_cut), */
      /*SCIProwGetNorm(cluster_cut), SCIPgetRowLPFeasibility(scip, cluster_cut), */
      /*SCIPgetCutEfficacy(scip, NULL, cluster_cut)); */

      /*/*SCIP_CALL(  SCIPprintRow(scip, cycle_cut, NULL)  ); */
      /*printf("Cycle cut: act=%f, rhs=%f, vars=%d, norm=%f, feas=%f, eff=%f\n", */
      /*SCIPgetRowLPActivity(scip, cycle_cut), SCIProwGetRhs(cycle_cut), SCIProwGetNNonz(cycle_cut), */
      /*SCIProwGetNorm(cycle_cut), SCIPgetRowLPFeasibility(scip, cycle_cut), */
      /*SCIPgetCutEfficacy(scip, NULL, cycle_cut)); */
      /*printf("\n"); */

      SCIP_CALL( SCIPreleaseRow(scip, &cluster_cut) );
      SCIP_CALL( SCIPreleaseRow(scip, &cycle_cut) );
   }

   return SCIP_OKAY;
}

/** Find cuts to rule out cycles in the current solution.
 *
 *  @param scip The SCIP instance the cuts are to be found for.
 *  @param conshdlr The constraint handler responsible for the cuts.
 *  @param sol The current relaxed solution.
 *  @param nGen The number of cuts that were added.
 *  @param forcecuts Whether to force cuts to be added
 *  @param found_efficacious_ptr A pointer to return whether any efficacious cuts were found.
 *  @return SCIP_OKAY if the algorithm terminated correctly, or an appropriate error otherwise.
 */
SCIP_RETCODE CC_findCuts(
   SCIP* scip, 
   SCIP_CONSHDLR* conshdlr, 
   ParentSetData* psd, 
   SCIP_SOL* sol, 
   CircuitCutsStorage* ccs, 
   int* nGen,
   SCIP_Bool forcecuts,
   SCIP_Bool* found_efficacious_ptr,
   SCIP_Bool* cutoff
   )
{

   if( ccs->add_cluster_cuts || ccs->add_cycle_cuts || ccs->add_cluster_cuts_to_pool || ccs->add_cycle_cuts_to_pool )
   {
      int i;
      int j;
      int k;

      /* Initialise */
      for( i = 0; i < psd->n; i++ )
      {
         ccs->index_array->items[i] = -1;
         ccs->lowlink->items[i] = -1;
         for( j = 0; j < psd->n; j++ )
            ccs->adj_matrix[i][j] = 0;
      }
      VectorListClear(ccs->components);
      VectorListClear(ccs->cycles);

      /* Set adjacency matrix entries based on the current solution */
      for( i = 0; i < psd->n; i++ )
         for( j = 0; j < psd->nParentSets[i]; j++ )
            for( k = 0; k < psd->nParents[i][j]; k++ )
               ccs->adj_matrix[psd->ParentSets[i][j][k]][i] += SCIPgetSolVal(scip, sol, psd->PaVars[i][j]);

      /* Round the adjacency matrix entries if they are within tolerance levels or we are not using fractional cuts*/
      for( i = 0; i < psd->n; i++ )
         for( j = 0; j < psd->n; j++ )
            if( SCIPisIntegral(scip, ccs->adj_matrix[i][j]) )
            {
               if( ccs->adj_matrix[i][j] < 0.5 )
                  ccs->adj_matrix[i][j] = 0;
               else
                  ccs->adj_matrix[i][j] = 1;
            }
            else if( ccs->add_fractional_cuts == FALSE )
            {
               ccs->adj_matrix[i][j] = 0;
            }

      /* Find the strongly connected components */
      j = 0;
      for( i = 0; i < psd->n; i++ )
         if( ccs->index_array->items[i] == -1 )
            SCIP_CALL( findStronglyConnectedComponents(i, &j, psd->n, ccs->index_array, ccs->lowlink, ccs->s, ccs->adj_matrix, ccs->components) );

      /* Sort the components (no real need for this) */
      for( i = 0; i < ccs->components->size; i++ )
         SCIPsortInt(ccs->components->items[i]->items, ccs->components->items[i]->size);

      /* For each SCC, find the cycles in the component */
      for( i = 0; i < ccs->components->size; i++ )
      {
         Vector* this_component = ccs->components->items[i];
         for( j = 0; j < this_component->size; j++ )
         {
            for( k = 0; k < this_component->size; k++ )
            {
               VectorClear(ccs->b->items[this_component->items[k]]);
               ccs->blocked[this_component->items[k]] = FALSE;
            }
            circuit(scip, this_component, j, j, ccs->b, ccs->blocked, ccs->s, ccs->adj_matrix, ccs->cycles, 0,
               ccs->max_cycles, ccs->max_cycle_length);
         }
      }

      /* Do something with the cycles */
      addCuts(scip, conshdlr, psd, sol, ccs->included_cluster, ccs->excluded_cluster, ccs->included_cycle, ccs->excluded_cycle, ccs->cycles, forcecuts, found_efficacious_ptr, cutoff, ccs);

      /* Return */
      (*nGen) = ccs->cycles->size;

   }
   else
   {

      (*nGen) = 0;

   }

   return SCIP_OKAY;
}

