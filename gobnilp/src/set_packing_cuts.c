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
#include "probdata_bn.h"
#include "scip/scipdefplugins.h"





/** Adds a set packing constraint and stores it in the problem data structure
 */
static 
SCIP_Bool add_spc_constraint(
   SCIP* scip,
   ParentSetData* psd, 
   int* cluster, 
   int cluster_size, 
   int highest_index, 
   SCIP_Bool propagate,
   int* size,
   SCIP_Bool releasecons
   )
{
   int i;
   int k;
   int l;
   int ci;
   int ci2;

   int parent;

   char cluster_name[SCIP_MAXSTRLEN];
   char tmp_str[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;

   int* new_cluster;

   int nvars = 0;
   int maxnvars = 0;
   SCIP_VAR** vars;
   SCIP_Bool found_one;
   SCIP_Bool ok;

   SCIP_PROBDATA* probdata;
   
   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   
   /* allocate (typically more than enough)  space for all variables in the constraint */
   for( ci = 0 ; ci < cluster_size ; ++ci )
      maxnvars += psd->nParentSets[cluster[ci]];
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvars) );

   /* construct name for the constraint */
   (void) SCIPsnprintf(cluster_name, SCIP_MAXSTRLEN, "spc(");
   for( ci = 0 ; ci < cluster_size ; ++ci )
   {
      (void) SCIPsnprintf(tmp_str, SCIP_MAXSTRLEN, "%s,", psd->nodeNames[cluster[ci]]);
      (void) strcat(cluster_name, tmp_str);
   }
   (void) strcat(cluster_name, ")");

   for( ci = 0; ci < cluster_size; ++ci )
   {
      i = cluster[ci];

      /* find all (parent sets for i containing
         all of the other elements in cluster
      */
      found_one = FALSE;
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         SCIP_Bool all_there = TRUE;
         for( ci2 = 0; ci2 < cluster_size; ++ci2 )
         {
            if( ci == ci2 )
               continue;

            parent = cluster[ci2];
            ok = FALSE;
            for( l = 0; l < psd->nParents[i][k]; l++)
            {
               if( psd->ParentSets[i][k][l] == parent )
               {
                  ok = TRUE;
                  break;
               }
            }
            
            if( !ok )
            {
               all_there = FALSE;
               break;
            }
         }
         
         if( all_there )
         {
            found_one = TRUE;
            vars[nvars++] = psd->PaVars[i][k];
         }
      }

      if( !found_one )
      {
         /* no parent set for i in the constraint
            so give up on this cluster */
         break;
      }
   }

   if( found_one )
   {
      SCIP_CALL( SCIPcreateConsSetpack(
            scip,
            &cons,
            cluster_name,
            nvars,
            vars,   
            TRUE,        /*initial,*/
            TRUE,        /* separate, */
            TRUE,        /* enforce */
            TRUE,        /* check */
            propagate,   /* propagate */
            FALSE,       /* local */
            FALSE,       /* modifiable */
            FALSE,       /* dynamic */
            FALSE,       /* removable */
            FALSE        /* stickingatnode */
            ));

      SCIP_CALL( SCIPaddCons(scip, cons) );
      /*SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  );*/
      
      /* ensure enough space */
      if( probdata->nspc_conss + 1 > *size)
      {
          int newsize;

          newsize = SCIPcalcMemGrowSize(scip, probdata->nspc_conss + 1);
          SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->spc_conss, *size, newsize) );
          *size = newsize;
      }
      
      /* store the constraint, and associated info */
      SCIP_CALL( storeclustercut(scip, &(probdata->spc_conss[probdata->nspc_conss++]),
            cons, NULL, cluster_size, cluster, cluster_size - 1) );

      if( releasecons )
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &new_cluster, cluster_size + 1) );
      for( ci = 0; ci < cluster_size; ++ci )
         new_cluster[ci] = cluster[ci];
      for( i = highest_index + 1; i < psd->n; ++i )
      {
         new_cluster[cluster_size] = i;
         add_spc_constraint(scip, psd, new_cluster, cluster_size + 1, i, propagate, size, releasecons);
      }
      SCIPfreeBlockMemoryArray(scip, &new_cluster, cluster_size + 1);
   }

   /* clean up */
   SCIPfreeBufferArray(scip, &vars);

   return found_one;
}



SCIP_RETCODE SP_add_spc_constraints(
   SCIP* scip, 
   ParentSetData* psd
   )

{   
   int i;
   int j;
   int k;
   int cluster[3];
   int size = 0;
   SCIP_Bool propagate;
   SCIP_Bool addspc;
   SCIP_Bool pricing;
   
   SCIP_PROBDATA* probdata;

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/propagatespc", &propagate) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/addspc", &addspc) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/pricing", &pricing) );

   if( !addspc )
      return SCIP_OKAY;

   
   /* start with clusters of size 3, clusters of size 2 already dealt
      with by arrow constraints
   */
   for( i = 0; i < psd->n; i++ )
   {
      cluster[0] = i;
      for( j = i + 1; j < psd->n; j++ )
      {
         cluster[1] = j;
         for( k = j + 1; k < psd->n; k++ )
         {
            cluster[2] = k;
            /* if we are not pricing then the cons should be released here */
            add_spc_constraint(scip, psd, cluster, 3, k, propagate, &size, !pricing);
         }
      }
   }

   /* shrink arrays appropriately */

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->spc_conss, size, probdata->nspc_conss) );
   
   return SCIP_OKAY;
}
