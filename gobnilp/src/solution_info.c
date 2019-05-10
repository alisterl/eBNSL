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

#include "solution_info.h"

/** @file
 *  Contains cached information about a solution to make it more faster
 *  and more convenient to reaccess this information.
 */


/** information on which non-empty parent sets are positive in a sol
    (typically, not always, an LP solution)
   is used by a number of separators so here we find and store it
   Sometimes augment with additional information about 'the graph'
*/
SCIP_RETCODE SI_setsolinfo(
   SCIP* scip,                 /**< SCIP instance */
   SolutionInfo* solinfo,      /**< solution information */
   ParentSetData* psd,         /**< parent set data */
   SCIP_SOL* sol,              /**< solution from which to extract information */
   SCIP_Bool augmented,        /**< whether to compute and store 'extra' info */
   SCIP_Bool dummysol          /**< whether to ignore 'sol' and pretend all family variables set to 1 */
   )
{
   int i;
   int j;
   int k;
   int l;
   int* parent_set;

   assert( !dummysol || sol == NULL);
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->pa), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->npa), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->ch), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->nch), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->ispa), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->posvars), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->nposvars), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->lpsolvals), psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->pa[i]), psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->ch[i]), psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->ispa[i]), psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->posvars[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(solinfo->lpsolvals[i]), psd->nParentSets[i]) );
   }

   for( i = 0 ; i < psd->n ; ++i )
   {
      solinfo->nposvars[i] = 0;
      if( augmented )
      {
         for( j = 0 ; j < psd->n ; ++j )
            solinfo->ispa[i][j] = FALSE;
         solinfo->npa[i] = 0;
         solinfo->nch[i] = 0;
      }
   }

   for( i = 0 ; i < psd->n ; ++i )
   {
      if( dummysol )
      {
         for( k = 0 ; k < psd->nParentSets[i]; ++k )
            solinfo->lpsolvals[i][k] = 1.0;
      }
      else
      {
         SCIP_CALL( SCIPgetSolVals(scip, sol, psd->nParentSets[i], psd->PaVars[i], solinfo->lpsolvals[i]) );
      }
      for( k = 0 ; k < psd->nParentSets[i]; ++k )
      {
         if( psd->nParents[i][k] == 0 )
            continue;
         if( SCIPisPositive(scip, solinfo->lpsolvals[i][k]) )
         {
            solinfo->posvars[i][(solinfo->nposvars[i])++] = k;
            if( augmented )
            {
               parent_set = psd->ParentSets[i][k];
               for( l = 0; l < psd->nParents[i][k]; ++l )
               {
                  j = parent_set[l];
                  if( !solinfo->ispa[i][j] )
                  {
                     solinfo->ispa[i][j] = TRUE;
                     assert(solinfo->npa[j] > -1 && solinfo->npa[j] < psd->n);
                     solinfo->pa[i][(solinfo->npa[i])++] = j;
                     assert(solinfo->nch[j] > -1 && solinfo->nch[j] < psd->n);
                     solinfo->ch[j][(solinfo->nch[j])++] = i;
                  }
               }
            }
         }
      }
   }
   return SCIP_OKAY;
}

void SI_freesolinfo(
   SolutionInfo* solinfo, 
   int n
   )
{
   int i;
   for( i = 0; i < n; i++ )
   {
      SCIPfreeMemoryArray(scip, &(solinfo->pa[i]));
      SCIPfreeMemoryArray(scip, &(solinfo->ch[i]));
      SCIPfreeMemoryArray(scip, &(solinfo->ispa[i]));
      SCIPfreeMemoryArray(scip, &(solinfo->posvars[i]));
      SCIPfreeMemoryArray(scip, &(solinfo->lpsolvals[i]));
   }
   SCIPfreeMemoryArray(scip, &(solinfo->ispa));
   SCIPfreeMemoryArray(scip, &(solinfo->pa));
   SCIPfreeMemoryArray(scip, &(solinfo->npa));
   SCIPfreeMemoryArray(scip, &(solinfo->ch));
   SCIPfreeMemoryArray(scip, &(solinfo->nch));
   SCIPfreeMemoryArray(scip, &(solinfo->posvars));
   SCIPfreeMemoryArray(scip, &(solinfo->nposvars));
   SCIPfreeMemoryArray(scip, &(solinfo->lpsolvals));
}
