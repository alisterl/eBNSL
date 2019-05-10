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
 *  The technique involves finding cycles in the current fractional solution
 *  of the linear relaxation, and adding cuts to rule out some of these.
 */

#include "fractional_circuit_cuts.h"
#include "utils.h"

#define DEFAULT_ADD_CUTS TRUE
#define DEFAULT_ADD_CUTS_TO_POOL FALSE

/** A collection of storage items for finding fractional cuts */
typedef struct
{
   SCIP* myscip;
   int num_nodes;
   int** successors;
   int*  num_successors;
   SCIP_Real* distance_to;
   int* predecessor;
   SCIP_Bool* in_queue;
   SCIP_Real** adj_matrix;
} FC_info;

SCIP_RETCODE FC_addParams(
   SCIP* scip          /**< SCIP data structure */
   )
{

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/fractional_circuit_cuts/add_cluster_cuts",
                             "whether cluster cuts should be added based on fractional cycles in the current solution",
                             DEFAULT_ADD_CUTS
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/fractional_circuit_cuts/add_cluster_cuts_to_pool",
                             "whether cluster cuts should be added to the pool based on fractional cycles in the current solution",
                             DEFAULT_ADD_CUTS_TO_POOL
                            ));

   return SCIP_OKAY;
}

/* SCIP_RETCODE FC_initialise( */
/*    SCIP* scip,                  /\**< SCIP data structure *\/ */
/*    CircuitCutsStorage* ccs,  */
/*    ParentSetData* data */
/*    ) */
/* { */
/*    SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/fractional_circuit_cuts/add_cluster_cuts",          &ADD_CUTS) ); */
/*    SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/fractional_circuit_cuts/add_cluster_cuts_to_pool",  &ADD_CUTS_TO_POOL) ); */
/*    /\* Add any memory allocation here. *\/ */
/*    return SCIP_OKAY; */
/* } */

/* SCIP_RETCODE FC_finalise( */
/*    SCIP* scip,               /\**< SCIP data structure *\/ */
/*    CircuitCutsStorage* ccs,  */
/*    int n */
/*    ) */
/* { */
/*    /\* Add any memory deallocation here. *\/ */
/*    return SCIP_OKAY; */
/* } */


static 
void clearAdjMatrix(
   FC_info* fc_info    /**< Fractional cuts info/storage */ 
   )
{
   int i;
   int j;
   for( i = 0; i < fc_info->num_nodes; i++ )
      for( j = 0; j < fc_info->num_nodes; j++ )
         (fc_info->adj_matrix)[i][j] = 1;
}

static 
void setAdjMatrix(
   FC_info* fc_info,         /**< Fractional cuts info/storage */ 
   SCIP_SOL* sol, 
   ParentSetData* psd
   )
{
   int i;
   int j;
   int k;

   for( i = 0; i < psd->n; i++ )
      for( j = 0; j < psd->nParentSets[i]; j++ )
         for( k = 0; k < psd->nParents[i][j]; k++ )
            (fc_info->adj_matrix)[psd->ParentSets[i][j][k]][i] -= SCIPgetSolVal(fc_info->myscip, sol, psd->PaVars[i][j]);
}

static 
void roundAdjMatrix(
   FC_info* fc_info   /**< Fractional cuts info/storage */ 
   )
{
   int i;
   int j;

   for( i = 0; i < fc_info->num_nodes; i++ )
      for( j = 0; j < fc_info->num_nodes; j++ )
         if( SCIPisIntegral(fc_info->myscip, (fc_info->adj_matrix)[i][j]) )
         {
            if( (fc_info->adj_matrix)[i][j] < 0.5 )
               (fc_info->adj_matrix)[i][j] = 0;
            else
               (fc_info->adj_matrix)[i][j] = 1;
         }
}

static 
void initAdjMatrix(
   FC_info* fc_info, /**< Fractional cuts info/storage */ 
   SCIP_SOL* sol, 
   ParentSetData* psd
   )
{
   clearAdjMatrix(fc_info);
   setAdjMatrix(fc_info, sol, psd);
   roundAdjMatrix(fc_info);
}

static 
void clearNumSuccessors(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      (fc_info->num_successors)[i] = 0;
}

static 
void countNumSuccessors(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   int j;
   for( i = 0; i < fc_info->num_nodes; i++ )
      for( j = 0; j < fc_info->num_nodes; j++ )
         if( (fc_info->adj_matrix)[i][j] != 1 )
            (fc_info->num_successors)[i] += 1;
}

static 
SCIP_RETCODE allocateSuccessorsMemory(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      SCIP_CALL( SCIPallocMemoryArray(fc_info->myscip, &((fc_info->successors)[i]), (fc_info->num_successors)[i]) );
   return SCIP_OKAY;
}

static
void setSuccessors(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   int j;
   int next;

   for( i = 0; i < fc_info->num_nodes; i++ )
   {
      next = 0;
      for( j = 0; j < fc_info->num_nodes; j++ )
         if( (fc_info->adj_matrix)[i][j] != 1 )
         {
            (fc_info->successors)[i][next] = j;
            next += 1;
         }
   }
}

static 
SCIP_RETCODE findSuccessors(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   clearNumSuccessors(fc_info);
   countNumSuccessors(fc_info);
   SCIP_CALL( allocateSuccessorsMemory(fc_info) );
   setSuccessors(fc_info);
   return SCIP_OKAY;
}

static 
SCIP_RETCODE allocDataStructures(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   SCIP_CALL( SCIPallocMemoryArray(fc_info->myscip, &(fc_info->distance_to), fc_info->num_nodes) );
   SCIP_CALL( SCIPallocMemoryArray(fc_info->myscip, &(fc_info->predecessor), fc_info->num_nodes) );
   SCIP_CALL( SCIPallocMemoryArray(fc_info->myscip, &(fc_info->in_queue), fc_info->num_nodes) );
   SCIP_CALL( SCIPallocMemoryArray(fc_info->myscip, &(fc_info->num_successors), fc_info->num_nodes) );
   SCIP_CALL( SCIPallocMemoryArray(fc_info->myscip, &(fc_info->successors), fc_info->num_nodes) );
   return SCIP_OKAY;
}

static 
void freeDataStructures(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      SCIPfreeMemoryArray(fc_info->myscip, &((fc_info->successors)[i]));

   SCIPfreeMemoryArray(fc_info->myscip, &(fc_info->distance_to));
   SCIPfreeMemoryArray(fc_info->myscip, &(fc_info->predecessor));
   SCIPfreeMemoryArray(fc_info->myscip, &(fc_info->in_queue));
   SCIPfreeMemoryArray(fc_info->myscip, &(fc_info->num_successors));
   SCIPfreeMemoryArray(fc_info->myscip, &(fc_info->successors));
}

static 
void setDistanceTo(
   FC_info* fc_info,  /**< Fractional cuts info/storage */ 
   int node, 
   SCIP_Real value
   )
{
   (fc_info->distance_to)[node] = value;
}

static 
void updateDistance(
   FC_info* fc_info,  /**< Fractional cuts info/storage */ 
   int from_node, 
   int to_node
   )
{
   SCIP_Real new_distance = (fc_info->distance_to)[from_node] + (fc_info->adj_matrix)[from_node][to_node];
   if( new_distance < (fc_info->distance_to)[to_node] )
   {
      (fc_info->distance_to)[to_node] = new_distance;
      (fc_info->predecessor)[to_node] = from_node;
   }
}

static 
void updateNeighboursDistances(
   FC_info* fc_info,  /**< Fractional cuts info/storage */ 
   int from_node
   )
{
   int i;
   int* neighbours = (fc_info->successors)[from_node];
   int num_neighbours = (fc_info->num_successors)[from_node];
   for( i = 0; i < num_neighbours; i++ )
      updateDistance(fc_info, from_node, neighbours[i]);
}

static 
void clearPredecessorStructure(
   FC_info* fc_info  /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      (fc_info->predecessor)[i] = -1;
}

static 
void clearQueueStructure(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      (fc_info->in_queue)[i] = TRUE;
}

static 
void clearDistanceStructure(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      (fc_info->distance_to)[i] = SCIPinfinity(fc_info->myscip);
}

static 
void clearDataStructures(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   clearPredecessorStructure(fc_info);
   clearQueueStructure(fc_info);
   clearDistanceStructure(fc_info);
}

static 
void initialiseForStartNode(
   FC_info* fc_info, /**< Fractional cuts info/storage */ 
   int start
   )
{
   clearDataStructures(fc_info);
   setDistanceTo(fc_info, start, 0);
}

static 
SCIP_Bool queueEmpty(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   for( i = 0; i < fc_info->num_nodes; i++ )
      if( (fc_info->in_queue)[i] )
         return FALSE;
   return TRUE;
}

static 
int getClosestQueuedNode(
   FC_info* fc_info /**< Fractional cuts info/storage */ 
   )
{
   int i;
   SCIP_Real min_distance = SCIPinfinity(fc_info->myscip);
   int min_node = -1;
   for( i = 0; i < fc_info->num_nodes; i++ )
      if( (fc_info->in_queue)[i] && (fc_info->distance_to)[i] < min_distance )
      {
         min_distance = (fc_info->distance_to)[i];
         min_node = i;
      }
   return min_node;
}

static 
void removeFromQueue(
   FC_info* fc_info,      /**< Fractional cuts info/storage */ 
   int node
   )
{
   (fc_info->in_queue)[node] = FALSE;
}


static 
SCIP_Bool isLessThanOne(
   SCIP* scip,                /**< SCIP data structure */ 
   float value
   )
{
   return value < (1 - SCIPepsilon(scip));
}

static 
SCIP_Bool runDijkstra(
   FC_info* fc_info,        /**< Fractional cuts info/storage */ 
   SCIP* scip,              /**< SCIP data structure */
   int start
   )
{
   initialiseForStartNode(fc_info, start);
   updateNeighboursDistances(fc_info, start);
   setDistanceTo(fc_info, start, SCIPinfinity(fc_info->myscip));
   while( !queueEmpty(fc_info) )
   {
      int closest_node = getClosestQueuedNode(fc_info);
      if( closest_node == -1 )
      {
         /* Ran out of reachable nodes without finding a loop. */
         return FALSE;
      }
      else if( !isLessThanOne(scip, (fc_info->distance_to)[closest_node]) )
      {
         /* Any loop we find would be too long. */
         return FALSE;
      }
      else if( closest_node == start )
      {
         /* We have a loop and it is not too long. */
         return TRUE;
      }
      else
      {
         /* Do a normal step of the algorithm. */
         removeFromQueue(fc_info, closest_node);
         updateNeighboursDistances(fc_info, closest_node);
      }
   }
   return (isLessThanOne(scip, (fc_info->distance_to)[start]));
}

static 
int findCycleLength(
   FC_info* fc_info, /**< Fractional cuts info/storage */ 
   int start
   )
{
   int length = 1;
   int current_node;
   for( current_node = (fc_info->predecessor)[start]; current_node != start; current_node = (fc_info->predecessor)[current_node] )
      length += 1;
   return length;
}

static 
SCIP_RETCODE extractNodesInCycle(
   FC_info* fc_info, /**< Fractional cuts info/storage */ 
   int start, 
   int** cycle, 
   int* cycle_length
   )
{
   int current_node;
   int i;

   (*cycle_length) = findCycleLength(fc_info, start);
   SCIP_CALL( SCIPallocMemoryArray(scip, cycle, (*cycle_length)) );

   current_node = start;
   for( i = 0; i < (*cycle_length); i++ )
   {
      (*cycle)[i] = current_node;
      current_node = (fc_info->predecessor)[current_node];
   }

   return SCIP_OKAY;
}

static 
SCIP_RETCODE addClusterCut(
   SCIP* scip,                    /**< SCIP data structure */
   SCIP_ROW* cut, 
   SCIP_SOL* sol, 
   SCIP_Bool* was_added,
   SCIP_Bool forcecuts,
   SCIP_Bool* is_efficacious,
   SCIP_Bool* cutoff
   )
{
   SCIP_CALL( SCIPaddRow(scip, cut, forcecuts, cutoff) );
   (*is_efficacious) = SCIPisCutEfficacious(scip, sol, cut);
   (*was_added) = TRUE;
   return SCIP_OKAY;
}

static 
SCIP_Bool shouldIncludeVariableInCut(
   int* parent_set, 
   int parent_set_size, 
   int* nodes, 
   int num
   )
{
   int i;
   int j;
   for( i = 0; i < parent_set_size; i++ )
      for( j = 0; j < num; j++ )
         if( parent_set[i] == nodes[j] )
            return TRUE;
   return FALSE;
}

static 
SCIP_RETCODE createClusterCut(
   SCIP* scip,                            /**< SCIP data structure */
   SCIP_CONSHDLR* conshdlr, 
   SCIP_VAR** include, 
   int num_include, 
   SCIP_VAR** exclude, 
   int num_exclude, 
   int num, 
   int num_nodes_excluded, 
   SCIP_ROW** cut
   )
{
   const char* name = "clustercut";
   SCIP_Real lhs = -SCIPinfinity(scip);
   SCIP_Real rhs = num - num_nodes_excluded - 1;
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, cut, conshdlr, name, lhs, rhs, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, (*cut), num_include, include, 1) );
   SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, (*cut), num_exclude, exclude, -1) );
   return SCIP_OKAY;
}

static 
void partitionVariablesToInclude(
   ParentSetData* psd, 
   int node, 
   int* nodes, 
   int num, 
   SCIP_VAR** included, 
   int* num_vars_included, 
   SCIP_VAR** excluded, 
   int* num_vars_excluded
   )
{
   int i;

   for( i = 0; i < psd->nParentSets[node]; i++ )
   {
      SCIP_VAR* this_var = psd->PaVars[node][i];
      SCIP_Bool include = shouldIncludeVariableInCut(psd->ParentSets[node][i], psd->nParents[node][i], nodes, num);
      if( include )
      {
         included[(*num_vars_included)] = this_var;
         (*num_vars_included) = (*num_vars_included) + 1;
      }
      else
      {
         excluded[(*num_vars_excluded)] = this_var;
         (*num_vars_excluded) = (*num_vars_excluded) + 1;
      }
   }
}

static 
SCIP_RETCODE addClusterCutOverNodes(
   SCIP* scip,                            /**< SCIP data structure */
   SCIP_CONSHDLR* conshdlr, 
   SCIP_SOL* sol, 
   SCIP_VAR** included, 
   SCIP_VAR** excluded, 
   ParentSetData* psd, 
   int* nodes, 
   int num, 
   SCIP_Bool* was_added, 
   SCIP_Bool forcecuts, 
   SCIP_Bool* is_efficacious, 
   SCIP_Bool* cutoff
   )
{
   int i;
   int num_vars_included = 0;
   int num_vars_excluded = 0;
   int previous_num_included, previous_num_excluded;
   int num_included_added, num_excluded_added;
   int num_nodes_excluded = 0;
   SCIP_ROW* cut;

   SCIP_Bool addcuts;
   SCIP_Bool addcuts_to_pool;
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/fractional_circuit_cuts/add_cluster_cuts", &addcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/fractional_circuit_cuts/add_cluster_cuts_to_pool", &addcuts_to_pool) );

   for( i = 0; i < num; i++ )
   {
      previous_num_included = num_vars_included;
      previous_num_excluded = num_vars_excluded;
      partitionVariablesToInclude(psd, nodes[i], nodes, num, included, &num_vars_included, excluded, &num_vars_excluded);
      num_included_added = num_vars_included - previous_num_included;
      num_excluded_added = num_vars_excluded - previous_num_excluded;

      if( num_included_added > num_excluded_added )
      {
         num_vars_included = previous_num_included;
         num_nodes_excluded += 1;
      }
      else
      {
         num_vars_excluded = previous_num_excluded;
      }
   }

   SCIP_CALL( createClusterCut(scip, conshdlr, included, num_vars_included, excluded, num_vars_excluded, num, num_nodes_excluded, &cut) );
   if( addcuts )
      SCIP_CALL( addClusterCut(scip, cut, sol, was_added, forcecuts, is_efficacious, cutoff) );
   if( addcuts_to_pool )
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   return SCIP_OKAY;
}

SCIP_RETCODE FC_findCuts(
   SCIP* scip,                              /**< SCIP data structure */
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
   int i;

   FC_info fc;
   
   fc.myscip = scip;
   fc.num_nodes = psd->n;
   fc.adj_matrix = ccs->adj_matrix;

   (*nGen) = 0;
   (*found_efficacious_ptr) = FALSE;

   SCIP_CALL( allocDataStructures(&fc) );
   initAdjMatrix(&fc, sol, psd);
   SCIP_CALL( findSuccessors(&fc) );

   for( i = 0; i < psd->n; i++ )
   {
      SCIP_Bool is_cycle_cut = runDijkstra(&fc, scip, i);
      if( is_cycle_cut )
      {
         int* cycle = NULL;
         int length;
         SCIP_Bool added;
         SCIP_Bool efficacious = FALSE;
         extractNodesInCycle(&fc, i, &cycle, &length);
         SCIP_CALL( addClusterCutOverNodes(scip, conshdlr, sol, ccs->included_cluster, ccs->excluded_cluster, psd, cycle, length, &added, forcecuts, &efficacious, cutoff) );
         if( *cutoff )
         {
            freeDataStructures(&fc);
            return SCIP_OKAY;
         }
         if( added )
            (*nGen) = (*nGen) + 1;
         (*found_efficacious_ptr) = ((*found_efficacious_ptr) || efficacious);
      }
   }

   freeDataStructures(&fc);

   return SCIP_OKAY;
}


