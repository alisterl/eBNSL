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
 *  Contains functions related to managing ParentSetData.
 */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "parent_set_data.h"
#include "utils.h"
#include "vector.h"
#include "vectorlist.h"
#include "stack.h"
#include "scip/pub_var.h"


/** releases all variables from given hash map 
 *  NB. a similar function can be found in SCIP's
 *  heur_dualval.c
 */

static
SCIP_RETCODE releaseVariableHashmapEntries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         hashmap             /**< hashmap */
   )
{
   int nentries;
   int i;

   assert(scip != NULL);
   assert(hashmap != NULL);

   nentries = SCIPhashmapGetNEntries(hashmap);

   for( i = 0; i < nentries; ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      entry = SCIPhashmapGetEntry(hashmap, i);

      if( entry != NULL )
      {
            SCIP_VAR* var;
            var = (SCIP_VAR*) SCIPhashmapEntryGetImage(entry);
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }
   }

   return SCIP_OKAY;
}


/** Deallocates the memory associated with a ParentSetData structure.
 *
 * releasevars should only be used after solving is complete and the 
 * variables embedded in psd are no longer needed
 *
 *  @param scip The SCIP instance which the data refers to.
 *  @param psd A pointer to the data structure to deallocate.
 *  @param releasevars Whether to release the variables (family, arrow and edge)
 *  @return SCIP_OKAY if the allocation succeeded.
 */
SCIP_RETCODE PS_deallocateParentSetData(
   SCIP* scip, 
   ParentSetData** psd_pointer,
   SCIP_Bool releasevars
   )
{
   int i;
   int k;
   ParentSetData* psd;

   assert(scip != NULL);
   assert(psd_pointer != NULL);

   psd = *psd_pointer;

   assert(psd != NULL);
   
   if( releasevars )
   {
      SCIP_CALL( releaseVariableHashmapEntries(scip, psd->arrow) );
      SCIP_CALL( releaseVariableHashmapEntries(scip, psd->edge) );
      for( i = 0; i < psd->n; i++ )
         for( k = 0;  k < psd->nParentSets[i]; k++ )
            SCIP_CALL( SCIPreleaseVar(scip, &(psd->PaVars[i][k])) );
   }
      
   for( i = 0; i < psd->n; i++ )
   {
      for( k = 0;  k < psd->nParentSets[i]; k++ )
         SCIPfreeMemoryArray(scip, &(psd->ParentSets[i][k]));

      SCIPfreeMemoryArray(scip, &(psd->nParents[i]));
      SCIPfreeMemoryArray(scip, &(psd->PaVars[i]));
      SCIPfreeMemoryArray(scip, &(psd->ParentSets[i]));
      SCIPfreeMemoryArray(scip, &(psd->nodeNames[i]));
   }
   SCIPfreeMemoryArray(scip, &(psd->nParentSets));
   SCIPfreeMemoryArray(scip, &(psd->nParents));
   SCIPfreeMemoryArray(scip, &(psd->ParentSets));
   SCIPfreeMemoryArray(scip, &(psd->nodeNames));
   SCIPfreeMemoryArray(scip, &(psd->PaVars));
   
   SCIP_CALL( hashtablefreeArrow(scip, psd) );
   SCIPfreeMemory(scip, psd_pointer);

   return SCIP_OKAY;
}



/** Makes a copy of a ParentSetData structure.
 *
 *  @param scip The SCIP instance where the copy ('duplicate') will "live"
 *  @param original The original data structure.
 *  @param duplicate_pointer A pointer to the duplicated data structure.
 *  @param original_vars Whether to put the vars embedded in original into duplicate
 *  @return SCIP_OKAY if copying suceeded.
 */
static
SCIP_RETCODE PS_copyParentSetData_ExtraArgt(
   SCIP* scip, 
   ParentSetData* original, 
   ParentSetData** duplicate_pointer,
   SCIP_Bool original_vars
   )
{
   int i;
   int j;
   int k;
   int l;
   
   SCIP_VAR* arrow_i_j;
   SCIP_VAR* edge_i_j;

   ParentSetData* duplicate;
   
   assert(original != NULL);
   assert(original->nParents != NULL);
   assert(original->PaVars != NULL);
   assert(original->ParentSets != NULL);
   assert(original->arrow != NULL);
   assert(original->edge != NULL);

   /* Allocate the memory */
   SCIP_CALL( SCIPallocMemory(scip, duplicate_pointer) );

   duplicate = *duplicate_pointer;
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->nParentSets), original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->nParents),    original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->ParentSets),  original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->nodeNames),   original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->PaVars),      original->n) );
   for( i = 0; i < original->n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->nParents[i]),   original->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->PaVars[i]),     original->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->ParentSets[i]), original->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->nodeNames[i]),  SCIP_MAXSTRLEN) );
      for( k = 0;  k < original->nParentSets[i]; k++ )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(duplicate->ParentSets[i][k]), original->nParents[i][k]) );
   }

   /* Copy the data across */
   duplicate->n = original->n;
   SCIP_CALL( hashtableCreateArrow(scip, duplicate) );

   for( i = 0; i < original->n; i++ )
   {
      duplicate->nParentSets[i] = original->nParentSets[i];
      for( k = 0;  k < original->nParentSets[i]; k++ )
      {
         if( original_vars )
            duplicate->PaVars[i][k] = original->PaVars[i][k];
         duplicate->nParents[i][k] = original->nParents[i][k];
         for( l = 0;  l < original->nParents[i][k]; l++ )
            duplicate->ParentSets[i][k][l] = original->ParentSets[i][k][l];
      }
      strcpy(duplicate->nodeNames[i], original->nodeNames[i]);

      if( original_vars )
      {
         /* set up hash tables for arrow and edge variables */
         for( j = 0;  j < original->n; j++ )
         {
            arrow_i_j = get_arrow(original,i,j);
            if( arrow_i_j != NULL )
               SCIP_CALL( put_arrow(scip, duplicate, i, j, arrow_i_j) );

            /* avoid attempting to add an edge variable twice */
            if( j > i )
            {
               edge_i_j = get_edge(original,i,j);
               if( edge_i_j != NULL )
                  SCIP_CALL( put_edge(scip, duplicate, i, j, edge_i_j) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** Makes a copy of a ParentSetData structure.
 *
 *  @param scip The SCIP instance where the copy ('duplicate') will "live"
 *  @param original The original data structure.
 *  @param duplicate_pointer A pointer to the duplicated data structure.
 *  @return SCIP_OKAY if copying suceeded.
 */
SCIP_RETCODE PS_copyParentSetData(
   SCIP* scip, 
   ParentSetData* original, 
   ParentSetData** duplicate_pointer
   )
{
   SCIP_CALL( PS_copyParentSetData_ExtraArgt(scip, original, duplicate_pointer, TRUE) );
   return SCIP_OKAY;
}



/** Makes a deep copy of a ParentSetData structure.
 *
 *  @param scip The SCIP instance to which the data belongs.
 *  @param original The original data structure.
 *  @param duplicate A pointer to the duplicated data structure.
 *  @return SCIP_OKAY if copying suceeded.
 */
SCIP_RETCODE PS_deepCopyParentSetData(
   SCIP*                 scip,               /**< target SCIP data structure */
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   ParentSetData*        original,           /**< The original data structure.*/
   ParentSetData**       duplicate_pointer,  /**< A pointer to the duplicated data structure. */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid               /**< pointer to store if the copying was valid */
   )
{
   int i;
   int j;
   int k;
   
   SCIP_VAR* arrow_i_j;
   SCIP_VAR* edge_i_j;
   SCIP_VAR* arrowcopy_i_j;
   SCIP_VAR* edgecopy_i_j;
   ParentSetData* duplicate;
   
   *valid = TRUE;

   SCIP_CALL( PS_copyParentSetData_ExtraArgt(scip, original, duplicate_pointer, FALSE) );

   duplicate = *duplicate_pointer;
   
   for( i = 0; i < original->n; i++ )
   {
      for( k = 0;  k < original->nParentSets[i]; k++ )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, original->PaVars[i][k],
               &(duplicate->PaVars[i][k]), varmap, consmap, global, valid) );
         assert(*valid && (duplicate->PaVars[i][k] != NULL));    
      }

      for( j = 0;  j < original->n; j++ )
      {
         arrow_i_j = get_arrow(original,i,j);
         if( arrow_i_j != NULL )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, arrow_i_j, &arrowcopy_i_j, varmap, consmap, global, valid) );
            assert(*valid && (arrowcopy_i_j != NULL));
            SCIP_CALL( put_arrow(scip, duplicate,i,j,arrowcopy_i_j) );
         }

         if( j > i )
         {
            edge_i_j = get_edge(original,i,j);
            if( edge_i_j != NULL )
            {
               SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, edge_i_j, &edgecopy_i_j, varmap, consmap, global, valid) );
               assert(*valid && (edgecopy_i_j != NULL));
               SCIP_CALL( put_edge(scip, duplicate,i,j,edgecopy_i_j) );
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** This function transforms the variables in the original parent set data and
    puts them into the transformed parent set data structure */
SCIP_RETCODE PS_transformParentSetData(
   SCIP* scip, 
   ParentSetData* original, 
   ParentSetData** transformed)
{
   int i;
   int j;
   int k;
   
   SCIP_VAR* arrow_i_j;
   SCIP_VAR* edge_i_j;
   SCIP_VAR* arrowtrans_i_j;
   SCIP_VAR* edgetrans_i_j;

   SCIP_CALL( PS_copyParentSetData_ExtraArgt(scip, original, transformed, FALSE) );

   for( i = 0; i < original->n; i++ )
      for( k = 0;  k < original->nParentSets[i]; k++ )
         SCIP_CALL( SCIPgetTransformedVar(scip, original->PaVars[i][k], &(*transformed)->PaVars[i][k]) );

   /* set up hash tables for arrow and edge variables */
   for( j = 0;  j < original->n; j++ )
   {
      arrow_i_j = get_arrow(original,i,j);
      if( arrow_i_j != NULL )
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, arrow_i_j, &arrowtrans_i_j) );
         assert(arrowtrans_i_j != NULL);
         SCIP_CALL( put_arrow(scip, *transformed, i, j, arrowtrans_i_j) );
      }
      if( j > i )
      {
         edge_i_j = get_edge(original,i,j);
         if( edge_i_j != NULL )
         {
            SCIP_CALL( SCIPgetTransformedVar(scip, edge_i_j, &edgetrans_i_j) );
            assert(edgetrans_i_j != NULL);
            SCIP_CALL( put_edge(scip, *transformed, i, j, edgetrans_i_j) );
         }
      }
   }

   return SCIP_OKAY;
}


/** Creates a subset of parent set data mentioning only the given nodes.
 *
 *  @param scip The SCIP instance the parent set data is used in.
 *  @param original The original data from which the subset will be taken.
 *  @param nodes The nodes which the new subset should be limited to.
 *  @param num_nodes The number of nodes the subset will be limited to.
 *  @param specialisation The subset of the parent set data containing only references to the requested nodes.
 *  @return SCIP_OKAY if the procedure worked, or an appropriate error otherwise.
 */
SCIP_RETCODE PS_specialiseFor(
   SCIP* scip, 
   ParentSetData* original, 
   int* nodes, 
   int num_nodes, 
   ParentSetData** specialisation
   )
{
   int i;
   int j;
   int k;

   int* label_map;

   /* Intermediate Storage */
   int n;
   int* nParentSets;
   int** nParents;
   int*** ParentSets;
   SCIP_VAR*** PaVars;
   char** nodeNames;

   SCIP_VAR* arrow_i_j;
   SCIP_VAR* edge_i_j;

   /* Allocate the intermediate storage */
   SCIP_CALL( SCIPallocMemoryArray(scip, &label_map,   original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParentSets, original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParents,    original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ParentSets,  original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeNames,   original->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &PaVars,      original->n) );
   for( i = 0; i < original->n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nParents[i]),   original->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(PaVars[i]),     original->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i]), original->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nodeNames[i]),  SCIP_MAXSTRLEN) );
      for( j = 0;  j < original->nParentSets[i]; j++ )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i][j]), original->nParents[i][j]) );
   }

   /* label_map[i] = -1 if i is to be discarded, otherwise label_map[i] is the position for i in the answer */
   for( i = 0; i < original->n; i++ )
      label_map[i] = -1;
   for( i = 0; i < num_nodes; i++ )
      if( nodes[i] <  original->n )
         label_map[nodes[i]] = i;

   /* Find the data to preserve */
   n = 0;
   for( i = 0; i < original->n; i++ )
   {
      if( label_map[i] != -1 )
      {
         n++;
         nParentSets[i] = 0;
         for( j = 0;  j < original->nParentSets[i]; j++ )
         {
            PaVars[i][j] = original->PaVars[i][j];
            nParents[i][j] = 0;
            for( k = 0;  k < original->nParents[i][j]; k++ )
            {
               if( label_map[original->ParentSets[i][j][k]] != -1 )
                  nParents[i][j]++;
               ParentSets[i][j][k] = label_map[original->ParentSets[i][j][k]];
            }
            /*            if ( nParents[i][j] > 0 ) */
            nParentSets[i]++;
         }
         strcpy(nodeNames[i], original->nodeNames[i]);
      }
   }

   /* Allocate the specialisation's storage */
   SCIP_CALL( SCIPallocMemory(scip, specialisation) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->nParentSets), n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->nParents),    n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->ParentSets),  n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->nodeNames),   n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->PaVars),      n) );
   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->nParents[i]),   nParentSets[nodes[i]]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->PaVars[i]),     nParentSets[nodes[i]]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->ParentSets[i]), nParentSets[nodes[i]]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->nodeNames[i]),  SCIP_MAXSTRLEN) );
      for( j = 0;  j < nParentSets[nodes[i]]; j++ )
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*specialisation)->ParentSets[i][j]), nParents[nodes[i]][j]) );
   }

   /* Copy the data across */
   (*specialisation)->n = n;
   SCIP_CALL( hashtableCreateArrow(scip, *specialisation) );

   for( i = 0; i < n; i++ )
   {
      (*specialisation)->nParentSets[i] = nParentSets[nodes[i]];
      for( j = 0;  j < nParentSets[nodes[i]]; j++ )
      {
         int current_k_slot = 0;
         (*specialisation)->nParents[i][j] = nParents[nodes[i]][j];
         (*specialisation)->PaVars[i][j] = PaVars[nodes[i]][j];
         for( k = 0; k < original->nParents[nodes[i]][j]; k++ )
            if( ParentSets[nodes[i]][j][k] != -1 )
            {
               (*specialisation)->ParentSets[i][j][current_k_slot] = ParentSets[nodes[i]][j][k];
               current_k_slot += 1;
            }
      }
      for( j = 0; j < n; j++ )
      {
         arrow_i_j = get_arrow(original,nodes[i],nodes[j]);
         if( arrow_i_j != NULL )
            SCIP_CALL( put_arrow(scip, *specialisation,i,j,arrow_i_j) );

         /* avoid attempting to add an edge variable twice */
         if( j > i )
         {
            edge_i_j = get_edge(original,nodes[i],nodes[j]);
            if( edge_i_j != NULL )
               SCIP_CALL( put_edge(scip, *specialisation,i,j,edge_i_j) );
         }
      }
      strcpy((*specialisation)->nodeNames[i], nodeNames[nodes[i]]);
   }

   /* Free the intermediate data structures */
   for( i = 0; i < original->n; i++ )
   {
      for( k = 0;  k < original->nParentSets[i]; k++ )
         SCIPfreeMemoryArray(scip, &(ParentSets[i][k]));
      SCIPfreeMemoryArray(scip, &(nParents[i]));
      SCIPfreeMemoryArray(scip, &(PaVars[i]));
      SCIPfreeMemoryArray(scip, &(ParentSets[i]));
      SCIPfreeMemoryArray(scip, &(nodeNames[i]));
   }
   SCIPfreeMemoryArray(scip, &(nParentSets));
   SCIPfreeMemoryArray(scip, &(nParents));
   SCIPfreeMemoryArray(scip, &(ParentSets));
   SCIPfreeMemoryArray(scip, &(nodeNames));
   SCIPfreeMemoryArray(scip, &(PaVars));
   SCIPfreeMemoryArray(scip, &(label_map));

   return SCIP_OKAY;
}

/** Finds strongly connected components in an adjacency matrix.
 *
 *  The algorithm used is Tarjan's Algorithm.  See the Wikipedia page
 *  for details.  See PS_splitToComponents() to see how the intermediate
 *  data structures should be intialised.
 *
 *  @param n The number of nodes.
 *  @param current_node The node to consider next.
 *  @param current_index The current depth index.
 *  @param adj_matrix The adjacency matrix to calculate on.
 *  @param index_array An intermediate data structure.
 *  @param lowlink An intermediate data structure.
 *  @param s An intermediate datastructure.
 *  @param components The strongly connected components found in the matrix.
 *  @return SCIP_OKAY if the procedure worked correctly, or an error otherwsie.
 */
static 
SCIP_RETCODE findStronglyConnectedComponents(
   int n, 
   int current_node, 
   int* current_index, 
   SCIP_Bool** adj_matrix, 
   Vector* index_array, 
   Vector* lowlink, 
   Stack* s, 
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
            findStronglyConnectedComponents(n, i, current_index, adj_matrix, index_array, lowlink, s, components);
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
      VectorListAppend(components, new_component);
   }

   return SCIP_OKAY;
}
/** Splits a single set of parent set data in to its strongly connected components.
 *
 *  If any variables are fixed at zero locally, the parent set the
 *  variable refers to is treated as if it did not exist for the
 *  purposes of finding the components.
 *
 *  @param scip The SCIP instance being used.
 *  @param original The parent set data to search for components in.
 *  @param num_components The number of strongly connected components found.
 *  @param components The strongly connected components of the data.
 *  @return SCIP_OKAY if the procedure worked, or an appropriate error message otherwise.
 */
SCIP_RETCODE PS_splitToComponents(
   SCIP* scip, 
   ParentSetData* original, 
   int* num_components, 
   ParentSetData*** components
   )
{
   int i;
   int j;
   int k;
   SCIP_Bool** adj_matrix;
   Vector* index_array = VectorCreate(original->n);
   Vector* lowlink = VectorCreate(original->n);
   VectorList* scc = VectorListCreate(original->n);
   Stack* s = StackCreate(original->n);

   /* Set up data structures */
   SCIP_CALL( SCIPallocMemoryArray(scip, &adj_matrix, original->n) );
   for( i = 0; i < original->n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(adj_matrix[i]), original->n) );

   for( i = 0; i < original->n; i++ )
   {
      index_array->items[i] = -1;
      lowlink->items[i] = -1;
      for( j = 0; j < original->n; j++ )
         adj_matrix[i][j] = FALSE;
   }

   /* Set adjacency matrix entries based on the current solution */
   for( i = 0; i < original->n; i++ )
      for( j = 0; j < original->nParentSets[i]; j++ )
         if( SCIPvarGetUbLocal(original->PaVars[i][j]) > 0.5 )
            for( k = 0; k < original->nParents[i][j]; k++ )
               adj_matrix[original->ParentSets[i][j][k]][i] = 1;

   /* Find the strongly connected components */
   j = 0;
   for( i = 0; i < original->n; i++ )
      if( index_array->items[i] == -1 )
         SCIP_CALL( findStronglyConnectedComponents(original->n, i, &j, adj_matrix, index_array, lowlink, s, scc) );

   /* Create the new parent set data structures */
   (*num_components) = scc->size;
   SCIP_CALL( SCIPallocMemoryArray(scip, components, (*num_components)) );
   for( i = 0; i < (*num_components); i++ )
      SCIP_CALL( PS_specialiseFor(scip, original, scc->items[i]->items, scc->items[i]->size, &((*components)[i])) );

   /* Deallocate the memory */
   for( i = 0; i < original->n; i++ )
      SCIPfreeMemoryArray(scip, &(adj_matrix[i]));
   SCIPfreeMemoryArray(scip, &adj_matrix);
   VectorDelete(&index_array);
   VectorDelete(&lowlink);
   VectorListDelete(&scc);
   StackDelete(&s);

   return SCIP_OKAY;
}

/** Writes a ParentSetData structure to file.
 *
 *  @param scip The SCIP instance the data belongs to.
 *  @param file The file to write to.
 *  @param psd The data to write.
 *  @return SCIP_OKAY if the writing succeeded.
 */
SCIP_RETCODE PS_writeToFile(
   SCIP* scip, 
   FILE* file, 
   ParentSetData* psd
   )
{
   if( psd != NULL )
   {
      int i;
      int j;

      int** arrows_idx;
      int** edges_idx;
      SCIP_VAR*** arrows_var;
      SCIP_VAR*** edges_var;
      SCIP_VAR* var;
      int* n_arrows_idx;
      int* n_edges_idx;


      SCIP_CALL( SCIPallocMemoryArray(scip, &arrows_idx, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &edges_idx, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &arrows_var, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &edges_var, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &n_arrows_idx, psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &n_edges_idx, psd->n) );


      SCIPinfoMessage(scip, file, "%c", UT_LIST_START);
      SCIPinfoMessage(scip, file, "%d", psd->n);
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeStringArray(scip, file, psd->nodeNames, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeIntArray(scip, file, psd->nParentSets, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeIntArrayArray(scip, file, psd->nParents, psd->nParentSets, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      SCIP_CALL( UT_writeIntArrayArrayArray(scip, file, psd->ParentSets, psd->nParents, psd->nParentSets, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);
      if( psd->PaVars != NULL )
         SCIP_CALL( UT_writeVarArrayArray(scip, file, psd->PaVars, psd->nParentSets, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);

      /* construct arrays and print them out */
      for (i = 0; i < psd->n; i++)
      {
         n_arrows_idx[i] = 0;
         n_edges_idx[i] = 0;

         SCIP_CALL( SCIPallocMemoryArray(scip, &(arrows_idx[i]), psd->n) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(edges_idx[i]), psd->n) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(arrows_var[i]), psd->n) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(edges_var[i]), psd->n) );

         for (j = 0; j < psd->n; j++)
         {
            var = get_arrow(psd, i, j);
            if( var != NULL )
            {
               arrows_idx[i][n_arrows_idx[i]] = j;
               arrows_var[i][n_arrows_idx[i]++] = var;
            }
         }
         for (j = 0; j < psd->n; j++)
         {
            var = get_edge(psd, i, j);
            if( var != NULL )
            {
               edges_idx[i][n_edges_idx[i]] = j;
               edges_var[i][n_edges_idx[i]++] = var;
            }
         }
      }
      SCIP_CALL( UT_writeIntArray(scip, file, n_arrows_idx, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);  
      SCIP_CALL( UT_writeIntArrayArray(scip, file, arrows_idx, n_arrows_idx, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);  
      SCIP_CALL( UT_writeVarArrayArray(scip, file, arrows_var, n_arrows_idx, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);  
      SCIP_CALL( UT_writeIntArray(scip, file, n_edges_idx, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);  
      SCIP_CALL( UT_writeIntArrayArray(scip, file, edges_idx, n_edges_idx, psd->n) );
      SCIPinfoMessage(scip, file, "%c ", UT_LIST_SEP);  
      SCIP_CALL( UT_writeVarArrayArray(scip, file, edges_var, n_edges_idx, psd->n) ); 

      SCIPinfoMessage(scip, file, "%c", UT_LIST_END);

      for( i = 0; i < psd->n; i++ )
      {
         SCIPfreeMemoryArray(scip, &(arrows_idx[i]));
         SCIPfreeMemoryArray(scip, &(arrows_var[i]));
         SCIPfreeMemoryArray(scip, &(edges_idx[i]));
         SCIPfreeMemoryArray(scip, &(edges_var[i]));
      }
      
      SCIPfreeMemoryArray(scip, &n_arrows_idx);
      SCIPfreeMemoryArray(scip, &arrows_idx);
      SCIPfreeMemoryArray(scip, &arrows_var);
      SCIPfreeMemoryArray(scip, &n_edges_idx);
      SCIPfreeMemoryArray(scip, &edges_idx);
      SCIPfreeMemoryArray(scip, &edges_var);

   }
   return SCIP_OKAY;
}
/** Parses a ParentSetData structure from a sting.
 *
 *  @param scip The SCIP instance the data will belong to.
 *  @param str The string to parse.
 *  @param psd A pointer to the data structure resulting from parsing.
 *  @return SCIP_OKAY if parsing succeeded.
 */
SCIP_RETCODE PS_parse(
   SCIP* scip, 
   char* str, 
   ParentSetData** psd
   )
{
   int i;
   int j;
   char** all_data = NULL;

   int n;
   char** nodeNames = NULL;
   int* nParentSets = NULL;
   int** nParents = NULL;
   int*** ParentSets = NULL;
   SCIP_VAR*** PaVars = NULL;

   int** arrows_idx;
   int** edges_idx;
   SCIP_VAR*** arrows_var;
   SCIP_VAR*** edges_var;
   int* n_arrows_idx;
   int* n_edges_idx;


   /* Get the data as string arrays */
   SCIP_CALL( SCIPallocMemoryArray(scip, &all_data, 12) );
   for( i = 0; i < 12; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(all_data[i]), strlen(str)) );
   SCIP_CALL( UT_parseArray((char*)str, &all_data) );

   /* Find the number of nodes */
   assert(all_data != NULL);
   assert(all_data[0] != NULL);
   sscanf(all_data[0], "%d", &n);

   /* Get the node names */
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodeNames, n) );
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nodeNames[i]), strlen(str)) );
   SCIP_CALL( UT_parseStringArray(all_data[1], &nodeNames, n) );

   /* Get the number of parent sets for each node */
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParentSets, n) );
   SCIP_CALL( UT_parseIntArray(all_data[2], &nParentSets, n) );

   /* Get the number of parents in each parent set */
   SCIP_CALL( SCIPallocMemoryArray(scip, &nParents, n) );
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(nParents[i]), nParentSets[i]) );
   SCIP_CALL( UT_parseIntArrayArray(all_data[3], &nParents, nParentSets, n) );

   /* Get the parents in each parent set */
   SCIP_CALL( SCIPallocMemoryArray(scip, &ParentSets, n) );
   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i]), nParentSets[i]) );
      for( j = 0; j < nParentSets[i]; j++ )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i][j]), nParents[i][j]) );
   }
   SCIP_CALL( UT_parseIntArrayArrayArray(all_data[4], &ParentSets, nParents, nParentSets, n) );

   /* Get the parent variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &PaVars, n) );
   for( i = 0; i < n; i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &(PaVars[i]), nParentSets[i]) );
   SCIP_CALL( UT_parseVarArrayArray(scip, all_data[5], &PaVars, nParentSets, n) );

   /* Get the arrow variables */

   SCIP_CALL( SCIPallocMemoryArray(scip, &arrows_idx, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &arrows_var, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_arrows_idx, n) );

   SCIP_CALL( UT_parseIntArray(all_data[6], &n_arrows_idx, n) );

   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(arrows_idx[i]), n_arrows_idx[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(arrows_var[i]), n_arrows_idx[i]) );
   }
   SCIP_CALL( UT_parseIntArrayArray(all_data[7], &arrows_idx, n_arrows_idx, n) );
   SCIP_CALL( UT_parseVarArrayArray(scip, all_data[8], &arrows_var, n_arrows_idx, n) );

   /* Get the edge variables */

   SCIP_CALL( SCIPallocMemoryArray(scip, &edges_idx, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &edges_var, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_edges_idx, n) );

   SCIP_CALL( UT_parseIntArray(all_data[9], &n_edges_idx, n) );

   for( i = 0; i < n; i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(edges_idx[i]), n_edges_idx[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(edges_var[i]), n_edges_idx[i]) );
   }
   SCIP_CALL( UT_parseIntArrayArray(all_data[10], &edges_idx, n_edges_idx, n) );
   SCIP_CALL( UT_parseVarArrayArray(scip, all_data[11], &edges_var, n_edges_idx, n) );

   /* Set the constraint data */
   SCIP_CALL( SCIPallocMemory(scip, psd) );
   (*psd)->n = n;
   (*psd)->nParentSets = nParentSets;
   (*psd)->nParents = nParents;
   (*psd)->ParentSets = ParentSets;
   (*psd)->nodeNames = nodeNames;
   (*psd)->PaVars = PaVars;

   /* set up hash table for arrow and edge variables */
   SCIP_CALL( hashtableCreateArrow(scip, *psd) );

   /* insert arrow variables */
   for( i = 0; i < n; i++ )
      for( j = 0; j < n_arrows_idx[i]; j++ )
         SCIP_CALL( put_arrow(scip, *psd, i, arrows_idx[i][j], arrows_var[i][j]) );
   /* insert edge variables */
   for( i = 0; i < n; i++ )
      for( j = 0; j < n_edges_idx[i]; j++ )
         SCIP_CALL( put_edge(scip, *psd, i, edges_idx[i][j], edges_var[i][j]) );

   /* Clean up */
   for( i = 0; i < 12; i++ )
      SCIPfreeMemoryArray(scip, &(all_data[i]));
   
   for( i = 0; i < n; i++ )
   {
      SCIPfreeMemoryArray(scip, &(arrows_idx[i]));
      SCIPfreeMemoryArray(scip, &(arrows_var[i]));
      SCIPfreeMemoryArray(scip, &(edges_idx[i]));
      SCIPfreeMemoryArray(scip, &(edges_var[i]));
   }

   SCIPfreeMemoryArray(scip, &n_arrows_idx);
   SCIPfreeMemoryArray(scip, &arrows_idx);
   SCIPfreeMemoryArray(scip, &arrows_var);
   SCIPfreeMemoryArray(scip, &n_edges_idx);
   SCIPfreeMemoryArray(scip, &edges_idx);
   SCIPfreeMemoryArray(scip, &edges_var);


   SCIPfreeMemoryArray(scip, &all_data);

   return SCIP_OKAY;
}
