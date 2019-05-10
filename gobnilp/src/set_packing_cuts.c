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
 *  Functions to add "set packing" constraints for acyclicity constraints.
 */

#include <string.h>
#include "set_packing_cuts.h"
#include "utils.h"
#include "scip/scipdefplugins.h"


static 
SCIP_Bool add_spc_constraint(
   SCIP* scip, 
   ParentSetData* psd, 
   SCIP_Bool*** store, 
   int* cluster, 
   int cluster_size, 
   int highest_index, 
   SCIP_Bool propagate
   )
{
   int i;
   int k;
   int kk;
   int ci;
   int ci2;
   int* k_indices;
   int n_k_indices;
   int k_indices_size = 0;

   int parent;

   char cluster_name[SCIP_MAXSTRLEN];
   char tmp_str[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;

   SCIP_Bool first;
   SCIP_Bool ok = TRUE;

   int* new_cluster;

   /* find which element of the cluster has the most parent sets */
   for( ci = 0 ; ci < cluster_size ; ++ci )
      if( psd->nParentSets[cluster[ci]] > k_indices_size )
         k_indices_size = psd->nParentSets[cluster[ci]];
   SCIP_CALL( SCIPallocMemoryArray(scip, &k_indices, k_indices_size) );

   /* construct name for the constraint */
   (void) SCIPsnprintf(cluster_name, SCIP_MAXSTRLEN, "spc(");
   for( ci = 0 ; ci < cluster_size ; ++ci )
   {
      (void) SCIPsnprintf(tmp_str, SCIP_MAXSTRLEN, "%d,", cluster[ci]);
      (void) strcat(cluster_name, tmp_str);
   }
   (void) strcat(cluster_name, ")");

   /* lazily create linear constraint, SCIP will upgrade it to set packing */
   SCIP_CALL(SCIPcreateConsLinear(scip, &cons, cluster_name, 0, NULL, NULL,
                                  -SCIPinfinity(scip),
                                  1,
                                  TRUE,
                                  TRUE,  /* separation here */
                                  FALSE,  /* no need to enforce */
                                  FALSE,  /* no need to check */
                                  propagate,
                                  FALSE, FALSE, FALSE, FALSE, FALSE));


   for( ci = 0; ci < cluster_size; ++ci )
   {
      i = cluster[ci];

      /* find all (indices of) parent sets for i containing
         all of the other elements in cluster
      */
      first = TRUE;
      n_k_indices = 0;
      for( ci2 = 0; ci2 < cluster_size; ++ci2 )
      {
         if( ci == ci2 )
            continue;

         parent = cluster[ci2];
         if( first )
         {
            /* for first 'parent' just copy ( indices of ) parent
               sets containing the 'parent' into k_indices
            */
            for( k = 0; k < psd->nParentSets[i]; ++k )
               if( store[i][k][parent] )
                  k_indices[n_k_indices++] = k;
            first = FALSE;
         }
         else
         {
            /* for non-first parents, remove from k_indices
               those ( indices of ) parents sets not containing parent
            */
            kk = 0;
            while( kk < n_k_indices )
            {
               if( !store[i][k_indices[kk]][parent] )
               {
                  /* parent missing, remove entry from k_indices */
                  n_k_indices--;
                  if( kk < n_k_indices )
                     k_indices[kk] = k_indices[n_k_indices];
               }
               else
                  kk++;
            }
         }
      }
      if( n_k_indices == 0 )
      {
         ok = FALSE;
         break;
      }
      /* now have correct indices in k_indices*/
      for( kk = 0; kk < n_k_indices; ++kk )
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, psd->PaVars[i][k_indices[kk]], 1) );
   }

   if( ok )
   {
      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
   }
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIPfreeMemoryArray(scip, &k_indices);

   /* if ok looking for 'bigger' spc constraints */
   if( ok )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &new_cluster, cluster_size + 1) );
      for( ci = 0; ci < cluster_size; ++ci )
         new_cluster[ci] = cluster[ci];
      for( i = highest_index + 1; i < psd->n; ++i )
      {
         new_cluster[cluster_size] = i;
         add_spc_constraint(scip, psd, store, new_cluster, cluster_size + 1, i, propagate);
      }
      SCIPfreeMemoryArray(scip, &new_cluster);
   }

   return ok;
}


SCIP_RETCODE SP_add_spc_constraints(
   SCIP* scip, 
   ParentSetData* psd, 
   SCIP_Bool*** store, 
   SCIP_Bool propagate
   )

{   
   int i;
   int j;
   int k;
   int cluster[3];
   
   /* start with clusters of size 3, clusters of size 2 already dealt
      with by edge variables
   */
   for( i = 0; i < psd->n; i++ )
   {
      cluster[0] = i;
      for( j = i + 1; j < psd->n; j++ )
      {
         /* if any arrows missing then we will not get a set-packing constraint */
         if( get_arrow(psd,i,j) == NULL || get_arrow(psd,j,i) == NULL )
            continue;
         
         cluster[1] = j;
         for( k = j + 1; k < psd->n; k++ )
         {
            /* if any arrows missing then we will not get a set-packing constraint */
            if( get_arrow(psd,i,k) == NULL || get_arrow(psd,k,i) == NULL 
               || get_arrow(psd,j,k) == NULL || get_arrow(psd,k,j) == NULL )
               continue;

            cluster[2] = k;
            add_spc_constraint(scip, psd, store, cluster, 3, k, propagate);
         }
      }
   }
   return SCIP_OKAY;
}
