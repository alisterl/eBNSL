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
 *  Code using a sub-IP to find "cluster constraint" cutting planes
 */

/*#define SCIP_DEBUG*/
#include "subip_cuts.h"
#include "scip/scipdefplugins.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

#if 0
/** Suppose a cluster cut has been found for cluster \f$\{a,b,c,d,e\}\f$. This corresponds to a uniform matroid
 * of rank 1 with circuits \f$\{ab,ac,ad,ae,bc,bd,be,cd,ce,de\}\f$. Given any pair from this cluster a connected matroid
 * of rank 2 can be constructed by joining together circuits containing either of the pair. For example, if the pair were a and b the 
 * relevant connected matroid (as represented by its circuits) would be \f$\{abc,abd,abe,cd,ce,de\}\f$. This function constructs 
 * such matroid inequalities in the hope they are useful cuts
 */
static SCIP_RETCODE AddMatroidPairCuts(
   SCIP*           scip,                  /**< (Main) SCIP instance */
   const ParentSetData*  psd,             /**< family variable information */
   SCIP_CONSHDLR*  conshdlr,              /**< the constraint handler responsible for adding these cuts (will be 'dagcluster') */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   SCIP_Bool*      incluster,             /**< the cluster itself: incluster[i] = TRUE iff i is in the cluster */
   int             cluster_size,          /**< the size of the found cluster */
   SCIP_Bool*      found_efficacious_ptr, /**< for recording whether the cutting plane is efficacious */
   SCIP_Bool*      cutoff                 /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
)
{
   SCIP_ROW* cut;
   int rhs;
   int i;
   int k;
   int l;
   int n_included;
   int n_excluded;
   SCIP_VAR** included;
   SCIP_VAR** excluded;
   SCIP_Bool include_in_cut;
   SCIP_Bool founda;
   SCIP_Bool foundb;
   SCIP_Bool foundx;


   int* parent_set;
   int parent;
   
   int a;   /* first element in the pair */
   int b;   /* second element in the pair */

   int ca;
   int cb;
   int ci;
   
   int* cluster;
   int c;

   int maxnpas;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &cluster, cluster_size) );
   c = 0;
   maxnpas = 0;
   for( i = 0 ; i < psd->n ; ++i )
   {
      if( incluster[i] )
         cluster[c++] = i;
      maxnpas = max(maxnpas,psd->nParentSets[i]);
   }
   assert(c ==  cluster_size);
   SCIP_CALL( SCIPallocBufferArray(scip, &included, maxnpas) );
   SCIP_CALL( SCIPallocBufferArray(scip, &excluded, maxnpas) );

   for( ca = 0 ; ca < cluster_size ; ++ca )
   {
      a = cluster[ca];
         
      for( cb = ca+1 ; cb < cluster_size ; ++cb )
         {

            b = cluster[cb];
            /* since matroid is of rank 2 get a tighter RHS for inequality */
            rhs = cluster_size-2;
            
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "clustercut_ab",
                  -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
            
               for( ci = 0 ; ci < cluster_size ; ++ci )
               {
                  i = cluster[ci];

                  n_included = 0;
                  n_excluded = 0;

                  for( k = 0;  k < psd->nParentSets[i]; ++k )
                  {
                     include_in_cut = FALSE;
                     founda = FALSE;
                     foundb = FALSE;
                     foundx = FALSE;
                     
                     parent_set = psd->ParentSets[i][k];
                     for( l = 0; l < psd->nParents[i][k]; ++l )
                     {
                        parent = parent_set[l];
                        if( parent == a )
                        {
                           founda = TRUE;
                        }
                        else if ( parent == b )
                        {
                           foundb = TRUE;
                        }
                        else if( incluster[parent] )
                        {
                           foundx = TRUE;
                        }
                        
                        if( foundx )
                        {
                           /* always need at least one non-pair cluster member
                              in the parent set to be included
                           */
                           if( i == a )
                           {
                              if( foundb )
                              {
                                 /* a always needs b in the parent set */
                                 include_in_cut = TRUE;
                                 break;
                              }
                           }
                           else if( i == b )
                           {
                              if( founda )
                              {
                                 /* b always needs a in the parent set */
                                 include_in_cut = TRUE;
                                 break;
                              }
                           }
                           else
                           {
                              /* cluster elements not in the pair just need one other
                                 non-pair cluster element in the parent set to be included
                              */
                              include_in_cut = TRUE;
                              break;
                           }
                        }
                     }
                     if( include_in_cut )
                        included[n_included++] = psd->PaVars[i][k];
                     else
                        excluded[n_excluded++] = psd->PaVars[i][k];
                  }

                  /* Use convexity constraint to reduce number of variables in the cut */
                  if( n_excluded < psd->nParentSets[i]  / 2 )
                  {
                     SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_excluded, excluded, -1.0) );
                     rhs--;
                  }
                  else
                     SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_included, included, 1.0) );

               }

               SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
               assert(SCIPisIntegral(scip, rhs));
               SCIPdebugMessage(" -> Matroid-pair-cut <clustercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                  SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                  SCIPgetCutEfficacy(scip, NULL, cut),
                  SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                  SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
               SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));

               if( SCIPisCutEfficacious(scip, sol, cut) )
               {
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
                  if( *cutoff )
                  {
                     SCIPdebugMessage("Matroid-pair cut led to cutoff\n");
                     SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                     return SCIP_OKAY;
                  }
                  *found_efficacious_ptr = TRUE;
               }
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         }
   }
   SCIPfreeBufferArray(scip, &excluded);
   SCIPfreeBufferArray(scip, &included);
   SCIPfreeBufferArray(scip, &cluster);
   return SCIP_OKAY;
}
#endif

/* finds family variables involved in a given cluster and create a specialised parent set data structure
 * 
 * also construct a non-Boolean representation of the cluster.
 * this function is very similar to PS_specialiseFor. They should be merged at some point
 */
static
SCIP_RETCODE FindFamilyVarsinClusterCut(
   SCIP*                scip,         /**< (Main) SCIP instance */
   const ParentSetData* psd,          /**< family variable information */
   const SCIP_Bool*     incluster,    /**< the cluster itself: incluster[i] = TRUE iff i is in the cluster */
   int                  cluster_size, /**< the size of the found cluster */
   int**                cluster,      /**< (*cluster)[i] will be the ith cluster member */
   ParentSetData** outpsd          /**< family variables in the cut */
   )
{

   int i;
   int k;
   int l;
   int c;
   int kc;
   int* tmp;
   int tmp_size;
   int n_tmp;
   int* invcluster;
   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, cluster, cluster_size) );
   SCIP_CALL( SCIPallocBlockMemory(scip, outpsd) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->nParentSets), cluster_size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->nParents), cluster_size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->ParentSets), cluster_size) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->PaVars), cluster_size) );
   (*outpsd)->nodeNames = NULL;
   (*outpsd)->arrow = NULL;
   (*outpsd)->edge = NULL;
   (*outpsd)->n = cluster_size;

   
   /* 10 is big enough for most problem instances, but we reallocate it not */
   tmp_size = 10;
   SCIP_CALL( SCIPallocBufferArray(scip, &tmp, tmp_size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &invcluster, psd->n) );

   c = 0;
   for( i = 0 ; i < psd->n ; ++i )
   {
      if( incluster[i] )
      {
         (*cluster)[c] = i;
         invcluster[i] = c++;
      }
      else
         invcluster[i] = -1;
   }

   c = 0;
   for( i = 0 ; i < psd->n ; ++i )
   {
      if( !incluster[i] )
         continue;
      kc = 0;
      
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->nParents[c]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->ParentSets[c]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->PaVars[c]), psd->nParentSets[i]) );
      
      for( k = 0;  k < psd->nParentSets[i]; ++k )
      {
         n_tmp = 0;
         if( psd->nParents[i][k] > tmp_size )
         {
            tmp_size = psd->nParents[i][k];
            SCIP_CALL( SCIPreallocBufferArray(scip, &tmp, tmp_size) );
         }
         for( l = 0; l < psd->nParents[i][k]; ++l )
            if ( incluster[psd->ParentSets[i][k][l]] ) 
               tmp[n_tmp++]  = invcluster[psd->ParentSets[i][k][l]];

         (*outpsd)->nParents[c][kc] = n_tmp;
         
         /* include all parents sets with at least one parent in the cluster */
         if( n_tmp > 0 )
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*outpsd)->ParentSets[c][kc]), n_tmp) );
            for( l = 0; l < n_tmp; ++l )
               (*outpsd)->ParentSets[c][kc][l] = tmp[l];
            
            (*outpsd)->PaVars[c][kc] = psd->PaVars[i][k];
            kc++;
         }
      }
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &((*outpsd)->nParents[c]), psd->nParentSets[i], kc) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &((*outpsd)->ParentSets[c]), psd->nParentSets[i], kc) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &((*outpsd)->PaVars[c]), psd->nParentSets[i], kc) );
      (*outpsd)->nParentSets[c] = kc;
      c++;
   }

   SCIPfreeBufferArray(scip, &invcluster);
   SCIPfreeBufferArray(scip, &tmp);
   
   return SCIP_OKAY;
}

static
SCIP_Bool stays(
   const ParentSetData* psd,
   int a,
   int b,
   int i,
   int k
   )
{
   int l;
   int pa;
   
   if( psd->nParents[i][k] == 1 )
   {
      if( i == a|| i == b
         || psd->ParentSets[i][k][0] == a
         || psd->ParentSets[i][k][0] == b )
         return FALSE;
      else
         return TRUE;
   }
   else
   {
      if( i != a && i != b )
         return TRUE;
      for( l = 0; l < psd->nParents[i][k]; ++l )
      {
         pa = psd->ParentSets[i][k][l];
         if( pa == a || pa == b )
            return TRUE;
      }
      return FALSE;
   }
}
      

static void inc_loss(
   SCIP_Real** loss,
   int a,
   int b,
   SCIP_Real inc
   )
{
   if(a < b)
      loss[a][b-a-1] += inc;
   else if (b < a)
      loss[b][a-b-1] += inc;
}


/** Find which variables in a cluster can be paired to give an efficacious cut */
static SCIP_RETCODE FindPairs(
   SCIP*           scip,                  /**< (Main) SCIP instance */
   const ParentSetData* old_psd,          /**< original family variable information */
   const ParentSetData* psd,              /**< specialised family variable information */
   const int* cluster,                    /**< the cluster */
   SCIP_CONSHDLR*  conshdlr,              /**< the constraint handler responsible for adding these cuts (will be 'dagcluster') */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   int             maxcuts,               /**< the maximum number of cuts to add (in this function) */
   SCIP_Bool       addtopool,             /**< whether to add the cut to the global cut pool */
   SCIP_Bool       forcecuts,             /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr, /**< for recording whether the cutting plane is efficacious */
   SCIP_Real       excess,
   SCIP_Bool*      cutoff                 /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
   )
{
   int i;
   int i0;
   int j;
   int k;
   int l;
   SCIP_Real val;
   SCIP_Real** loss;
   int pa;

   SCIP_ROW* cut;
   int old_i;
   SCIP_VAR** included;
   int n_included;
   SCIP_Real rhs;

   int ij;
   int ij0;
   int n_ij;
   int* i_i;
   int* j_i;
   int* k_i;
   SCIP_Real* flat_loss;
   const int npairs = (psd->n)*((psd->n)-1)/2;
   
   assert( psd != NULL );
   assert( psd->n > 0 );
   
   SCIP_CALL( SCIPallocBufferArray(scip, &loss, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(loss[i]), psd->n - i - 1) );
      for( j = 0; j < psd->n - i - 1; ++j )
         loss[i][j] = 0.0;
   }

   for( i = 0; i < psd->n; ++i )
   {
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {         
         val = SCIPgetSolVal(scip, sol, psd->PaVars[i][k]);

         if( psd->nParents[i][k] == 1 )
         {
            /* all pair matroid cuts where i is in the pair
               will not have variable psd->PaVars[i][k]) */
            
            /* all pair matroid cuts where the single parent is in the pair
               will not have variable psd->PaVars[i][k]) */

            pa = psd->ParentSets[i][k][0];
            for( j = 0; j < psd->n; ++j )
            {
               inc_loss(loss,i,j,val);
               inc_loss(loss,pa,j,val);
            }
            /* correct for double counting */
            inc_loss(loss,i,pa,-val);
         }
         else
         {
            /* loss for all pairs where i is in the pair and the other
               member is not in this parent set */
            for( j = 0; j < psd->n; ++j )
               inc_loss(loss,i,j,val);
            for( l = 0; l < psd->nParents[i][k]; ++l )
               inc_loss(loss,i,psd->ParentSets[i][k][l],-val);
         }
      }
   }

   n_ij = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &i_i, npairs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &j_i, npairs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &k_i, npairs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &flat_loss, npairs) );
   for( i = 0; i < psd->n; ++i )
      for( j = i+1; j < psd->n; ++j )
      {
         if( SCIPisLT(scip,loss[i][j-i-1],1.0+excess) )
         {
            /* worth adding a cut for pair {i,j} */
            i_i[n_ij] = i;
            j_i[n_ij] = j;
            k_i[n_ij] = n_ij;
            flat_loss[n_ij++] = loss[i][j-i-1];
         }
      }

   /* for( ij0 = 0; ij0 < n_ij; ++ij0 ) */
   /*    printf("loss is now %g\n", flat_loss[ij0]); */
   
   SCIPsortRealInt(flat_loss,k_i,n_ij);

   /* for( ij0 = 0; ij0 < n_ij; ++ij0 ) */
   /*    printf("loss is later %g\n", flat_loss[ij0]); */

   
   /* add up to min(maxcuts,n_ij) cuts, best first */
   for( ij0 = 0; ij0 < min(maxcuts,n_ij); ++ij0 )
   {
      ij = k_i[ij0];
      i = i_i[ij];
      j = j_i[ij];

      /* printf("loss is %g\n", flat_loss[ij0]); */
      assert( j > i );
      /* printf(" foo %g %g\n", flat_loss[ij0],loss[i][j-i-1]); */
      assert( SCIPisEQ(scip, flat_loss[ij0],loss[i][j-i-1]) );
      
      rhs = (psd->n)-2;
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "clustercut_ab",
            -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
            
      for( i0 = 0; i0 < psd->n; ++i0 )
      {
         n_included = 0;
         SCIP_CALL( SCIPallocBufferArray(scip, &included, psd->nParentSets[i0]) );
         for( k = 0; k < psd->nParentSets[i0]; ++k )
            if( stays(psd,i,j,i0,k) )
               included[n_included++] = psd->PaVars[i0][k];

         old_i = cluster[i0];
         
         if( n_included < old_psd->nParentSets[old_i] / 2 )
            SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_included, included, 1.0) );
         else
         {
            /* if too many variables to include, use all other variables in the original
               set of parent set variables for this child
               
               at present doing this poorly by scanning for variables
            */
            SCIP_VAR* var;
            SCIP_Bool foundvar;
            
            rhs--;
            
            SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );
            for( k = 0; k < old_psd->nParentSets[old_i]; ++k )
            {
               var = old_psd->PaVars[old_i][k];
               foundvar = FALSE;
               for( l = 0; l < n_included; ++l)
               {
                  if( included[l] == var )
                  {
                     foundvar = TRUE;
                     break;
                  }
               }
               if( !foundvar )
                  SCIP_CALL( SCIPaddVarToRow(scip, cut, var, -1.0) );
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, cut) );
         }
         SCIPfreeBufferArray(scip, &included);
      }
      
      SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
      assert(SCIPisIntegral(scip, rhs));
      SCIPdebugMessage(" -> Matroid-cut <clustercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
         SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
         SCIPgetCutEfficacy(scip, NULL, cut),
         SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
         SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
      SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));
      
      SCIP_CALL( SCIPaddRow(scip, cut, forcecuts, cutoff) );
      if( *cutoff )
      {
         SCIPdebugMessage("Matroid pair cut led to cutoff\n");
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         return SCIP_OKAY;
      }
      
      if( addtopool )
         SCIP_CALL( SCIPaddPoolCut(scip, cut) );
      
      if( SCIPisCutEfficacious(scip, sol, cut) )
         *found_efficacious_ptr = TRUE;
      
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }
   
   SCIPfreeBufferArray(scip, &flat_loss);
   SCIPfreeBufferArray(scip, &k_i);
   SCIPfreeBufferArray(scip, &j_i);
   SCIPfreeBufferArray(scip, &i_i);

   for( i = (psd->n)-1; i > -1; --i )
   {
      SCIPfreeBufferArray(scip, &(loss[i]));
   }
   SCIPfreeBufferArray(scip, &loss);

   
   return SCIP_OKAY;
}


/** Given a cluster, finds family variables for that cluster and adds the cluster cut */
static SCIP_RETCODE AddClusterCut(
   SCIP*           scip,                  /**< (Main) SCIP instance */
   const ParentSetData*  psd,             /**< family variable information */
   SCIP_CONSHDLR*  conshdlr,              /**< the constraint handler responsible for adding these cuts (will be 'dagcluster') */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   SCIP_Bool*      incluster,             /**< the cluster itself: incluster[i] = TRUE iff i is in the cluster */
   int             cluster_size,          /**< the size of the found cluster */
   int             kval,                  /**< kval = 1 for normal cluster cuts, kval = k for a 'k-cluster cut' */
   SCIP_Bool       addtopool,             /**< whether to add the cut to the global cut pool */
   SCIP_Bool       forcecuts,             /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr, /**< for recording whether the cutting plane is efficacious */
   SCIP_Bool       ci_cut,                /**< whether we are adding a CI cut */
   SCIP_Bool       matroidpaircuts,
   int             matroidpaircutslimit,
   int             matroidpairmaxcuts,
   SCIP_Real       excess,
   SCIP_Bool*      cutoff                 /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
)
{
   SCIP_ROW* cut;
   int rhs = cluster_size - kval;
   int i;
   int k;
   int l;
   int n_included;
   int n_excluded;
   SCIP_VAR** included;
   SCIP_VAR** excluded;
   int overlap;
   SCIP_Bool include_in_cut;
   int* parent_set;

   int maxnpas;

   ParentSetData*  outpsd;
   int* cluster;
   
   assert(psd != NULL);
   assert(psd->PaVars != NULL);

   /* CI cuts are tighter than cluster cuts */
   if( ci_cut )
      rhs--;

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "clustercut", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );

   maxnpas = 0;
   for( i = 0 ; i < psd->n ; ++i )
      if( incluster[i] )
         maxnpas = max(maxnpas,psd->nParentSets[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &included, maxnpas) );
   SCIP_CALL( SCIPallocBufferArray(scip, &excluded, maxnpas) );

   for( i = 0 ; i < psd->n ; ++i )
   {
      if( !incluster[i] )
         continue;

      n_included = 0;
      n_excluded = 0;


      /* include all parents sets with at least kval parents in cluster */
      for( k = 0;  k < psd->nParentSets[i]; ++k )
      {
         include_in_cut = FALSE;
         overlap = 0;
         parent_set = psd->ParentSets[i][k];
         for( l = 0; l < psd->nParents[i][k]; ++l )
         {
            /* if ( incluster[psd->ParentSets[i][k][l]] ) */
            if( incluster[parent_set[l]] )
            {
               overlap++;
               if( SCIPisGE(scip, overlap, kval) )
               {
                  include_in_cut = TRUE;
                  break;
               }
            }
         }

         if( include_in_cut )
            included[n_included++] = psd->PaVars[i][k];
                  else
            excluded[n_excluded++] = psd->PaVars[i][k];
      }
      
      /* Use convexity constraint to reduce number of variables in the cut */
      if( n_excluded < psd->nParentSets[i]  / 2 )
      {
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_excluded, excluded, -1.0) );
         rhs--;
      }
      else
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_included, included, 1.0) );

   }

   SCIPfreeBufferArray(scip, &excluded);
   SCIPfreeBufferArray(scip, &included);

   SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
   assert(SCIPisIntegral(scip, rhs));
   SCIPdebugMessage(" -> Cluster-cut <clustercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                    SCIPgetCutEfficacy(scip, NULL, cut),
                    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
   SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));

   SCIP_CALL( SCIPaddRow(scip, cut, forcecuts, cutoff) );
   if( *cutoff )
   {
      SCIPdebugMessage("Cluster cut led to cutoff\n");
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      return SCIP_OKAY;
   }

   if( addtopool )
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   if( SCIPisCutEfficacious(scip, sol, cut) )
      *found_efficacious_ptr = TRUE;

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   if( matroidpaircuts && cluster_size < matroidpaircutslimit )
   {
      SCIP_CALL( FindFamilyVarsinClusterCut(scip, psd, incluster, cluster_size, &cluster, &outpsd) );
      SCIP_CALL( FindPairs(scip, psd, outpsd, cluster, conshdlr, sol, matroidpairmaxcuts,
            addtopool, forcecuts, found_efficacious_ptr, excess, cutoff) );

      /* deallocate here, could move to separate function */
      SCIPfreeBlockMemoryArray(scip, &cluster, cluster_size);
      for( i = 0; i < cluster_size; ++i )
      {
         for( k = 0; k < outpsd->nParentSets[i]; ++k )
            SCIPfreeBlockMemoryArray(scip, &(outpsd->ParentSets[i][k]), outpsd->nParents[i][k]);
         
         SCIPfreeBlockMemoryArray(scip, &(outpsd->nParents[i]), outpsd->nParentSets[i]);
         SCIPfreeBlockMemoryArray(scip, &(outpsd->ParentSets[i]), outpsd->nParentSets[i]);
         SCIPfreeBlockMemoryArray(scip, &(outpsd->PaVars[i]), outpsd->nParentSets[i]);
      }
      SCIPfreeBlockMemoryArray(scip, &(outpsd->nParentSets), cluster_size);
      SCIPfreeBlockMemoryArray(scip, &(outpsd->nParents), cluster_size);
      SCIPfreeBlockMemoryArray(scip, &(outpsd->ParentSets), cluster_size);
      SCIPfreeBlockMemoryArray(scip, &(outpsd->PaVars), cluster_size);
      SCIPfreeBlockMemory(scip, &outpsd);
      
      /* SCIP_CALL( AddMatroidPairCuts(scip, psd, conshdlr, sol, incluster, cluster_size,  */
      /*       found_efficacious_ptr, cutoff) ); */
   }
   
   return SCIP_OKAY;
}

#if 0
/** Given the family variables for a given cluster, just adds the cluster cut */
static SCIP_RETCODE JustAddClusterCut(
   SCIP*           scip,                  /**< (Main) SCIP instance */
   SCIP_CONSHDLR*  conshdlr,              /**< the constraint handler responsible for adding these cuts (will be 'dagcluster') */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   int             cluster_size,          /**< the size of the found cluster */
   SCIP_Bool       addtopool,             /**< whether to add the cut to the global cut pool */
   SCIP_Bool       forcecuts,             /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr, /**< for recording whether the cutting plane is efficacious */
   SCIP_Bool*      cutoff,                 /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
   const int* n_included,     /**< n_included[c] is the number of family variables with cth cluster member included in the cut */
   const int* n_excluded,     /**< n_excluded[c] is the number of family variables with cth cluster member excluded from the cut */
   SCIP_VAR*** includedvars,  /**< includedvars[c][j] is the index of the jth family variable for cth cluster member included in the cut */
   SCIP_VAR*** excludedvars   /**< excludedvars[c][j] is the index of the jth family variable for cth cluster member excluded from the cut */
)
{

   int c;
   SCIP_ROW* cut;
   int rhs = cluster_size - 1;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( n_included != NULL );
   assert( n_excluded != NULL );
   assert( includedvars != NULL );
   assert( excludedvars != NULL );
   
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "clustercut", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );

   for( c = 0 ; c < cluster_size ; ++c )
   {
      if( n_excluded[c] < n_included[c] )
      {
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_excluded[c], excludedvars[c], -1.0) );
         rhs--;
      }
      else
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_included[c], includedvars[c], 1.0) );
   }
   SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
   assert(SCIPisIntegral(scip, rhs));
   SCIPdebugMessage(" -> Cluster-cut <clustercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                    SCIPgetCutEfficacy(scip, NULL, cut),
                    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
   SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));

   SCIP_CALL( SCIPaddRow(scip, cut, forcecuts, cutoff) );
   if( *cutoff )
   {
      SCIPdebugMessage("Cluster cut led to cutoff\n");
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      return SCIP_OKAY;
   }

   if( addtopool )
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   if( SCIPisCutEfficacious(scip, sol, cut) )
      *found_efficacious_ptr = TRUE;

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** Adds a found cluster cut ( 1-cluster, not a CI cut )*/
static SCIP_RETCODE AddBasicClusterCut(
   SCIP*           scip,                  /**< (Main) SCIP instance */
   const ParentSetData*  psd,             /**< family variable information */
   SCIP_CONSHDLR*  conshdlr,              /**< the constraint handler responsible for adding these cuts (will be 'dagcluster') */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   SCIP_Bool*      incluster,             /**< the cluster itself: incluster[i] = TRUE iff i is in the cluster */
   int             cluster_size,          /**< the size of the found cluster */
   SCIP_Bool       addtopool,             /**< whether to add the cut to the global cut pool */
   SCIP_Bool       forcecuts,             /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr, /**< for recording whether the cutting plane is efficacious */
   SCIP_Bool*      cutoff                 /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
)
{
   SCIP_ROW* cut;
   int rhs = cluster_size - 1;
   int i;
   int k;
   int l;
   int n_included;
   int n_excluded;
   SCIP_VAR** included;
   SCIP_VAR** excluded;
   SCIP_Bool include_in_cut;
   int maxnpas;
   
   assert(psd != NULL);
   assert(psd->PaVars != NULL);

   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "clustercut", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );

   maxnpas = 0;
   for( i = 0 ; i < psd->n ; ++i )
      if( incluster[i] )
         maxnpas = max(maxnpas,psd->nParentSets[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &included, maxnpas) );
   SCIP_CALL( SCIPallocBufferArray(scip, &excluded, maxnpas) );

   for( i = 0 ; i < psd->n ; ++i )
   {
      if( !incluster[i] )
         continue;

      n_included = 0;
      n_excluded = 0;

      /* include all parents sets with at least one parent in the cluster */
      for( k = 0;  k < psd->nParentSets[i]; ++k )
      {
         include_in_cut = FALSE;
         for( l = 0; l < psd->nParents[i][k]; ++l )
         {
            if ( incluster[psd->ParentSets[i][k][l]] ) 
            {
               include_in_cut = TRUE;
               break;
            }
         }

         if( include_in_cut )
            included[n_included++] = psd->PaVars[i][k];
         else
            excluded[n_excluded++] = psd->PaVars[i][k];
      }
      
      /* Use convexity constraint to reduce number of variables in the cut */
      if( n_excluded < psd->nParentSets[i]  / 2 )
      {
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_excluded, excluded, -1.0) );
         rhs--;
      }
      else
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_included, included, 1.0) );
   }

   SCIPfreeBufferArray(scip, &excluded);
   SCIPfreeBufferArray(scip, &included);

   SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
   assert(SCIPisIntegral(scip, rhs));
   SCIPdebugMessage(" -> Cluster-cut <clustercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                    SCIPgetCutEfficacy(scip, NULL, cut),
                    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
   SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));

   SCIP_CALL( SCIPaddRow(scip, cut, forcecuts, cutoff) );
   if( *cutoff )
   {
      SCIPdebugMessage("Cluster cut led to cutoff\n");
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
      return SCIP_OKAY;
   }

   if( addtopool )
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   if( SCIPisCutEfficacious(scip, sol, cut) )
      *found_efficacious_ptr = TRUE;

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}
#endif

/** main routine for looking for cutting planes

The number of found cutting planes is recorded in *nGen. A positive value indicates that the current solution is infeasible.
 */
extern SCIP_RETCODE IP_findCuts(
   SCIP*           scip,                    /**< SCIP data structure */
   ParentSetData*  psd,                     /**< family variable information */
   SolutionInfo*   solinfo,                 /**< information about the solution to be separated */
   SCIP_SOL*       sol,                     /**< solution to be separated */
   int*            nGen,                    /**< *nGen is number of cutting planes added ( even non-efficacious ones are added ) */
   int             k_lb,                    /**< lowerbound on 'k' values for k-cluster searching, always positive */
   int             k_ub,                    /**< upperbound on 'k' values for k-cluster searching */
   SCIP_CONSHDLR*  conshdlr,                /**< constraint handler */
   SCIP_Bool       addtopool,               /**< whether to add any found cut to the global cut pool */
   SCIP_Bool       forcecuts,               /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,   /**< to return whether an efficacious cutting plane was found */
   SCIP_Real       limits_time,             /**< limit on how long to spend sub-IP solving */
   SCIP_Real       limits_gap,              /**< limit on size of gap in sub-IP */
   SCIP_Real       limits_absgap,           /**< limit on size of the absolute gap in sub-IP */
   SCIP_Bool       incumbent_cons,          /**< whether to consider only cutting planes on which the incumbent lies */
   int*            must_be_included,        /**< set of nodes which must be included in any found cluster */
   int             n_must_be_included,      /**< size of the set of nodes which must be included in any found cluster */
   int*            must_be_excluded,        /**< set of nodes which must be excluded from any found cluster */
   int             n_must_be_excluded,      /**< size of the set of nodes which must be excluded from any found cluster */
   SCIP_Bool       ci_cut,                  /**< TRUE if we are looking for conditional independence cuts */
   SCIP_Bool       matroidpaircuts,
   int             matroidpaircutslimit,
   int             matroidpairmaxcuts,
   SCIP_Bool*      cutoff                   /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
)
{

   int i;
   int k;
   int l;
   SCIP_STATUS status;
   int s;
   int nsols;
   SCIP_SOL** subscip_sols;
   SCIP_SOL* subscip_sol;
   SCIP_SOL* best_sol;

   char consname[SCIP_MAXSTRLEN];
   char varname[SCIP_MAXSTRLEN];

   SCIP_Real kval;
   SCIP_Real val;

   SCIP_VAR** clausevars;
   int nvars;

   int cluster_size;
   SCIP_Bool* incluster;

   int ki;

   int* parent_set;

   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   SCIP*           auxipdata_subscip;          /**< sub MIP for finding good clusters for cutting planes */
   SCIP_VAR***     auxipdata_family;           /**< subIP variables: family[i][k] = 1 
                                                * if kth parent set of ith variable is one of those in cluster cut */
   SCIP_VAR**      auxipdata_incluster;        /**< subIP variables: incluster[i] if variable i in cluster */
   SCIP_VAR*       auxipdata_kvar;             /**< lower bound on number of parents to be in cluster for family variable to be set */
   SCIP_CONS***    auxipdata_clausecons;       /**< clausecons[i][k] is the constraint  "if family[i][k]=1 then incluster[i]=1" */
   SCIP_CONS***    auxipdata_overlapcons;      /**< overlapcons[i][k] is the constraint  "if family[i][k]=1 then \sum_{u \in W} >= kvar" */
   SCIP_CONS*      auxipdata_ck_cons;          /**< constraint for lb for cluster size (depends on kvar) */
   SCIP_CONS*      auxipdata_incumbentcons;    /**< (optional) constraint that incumbent must lie on cluster cut */

   SCIP_Real excess;
   
   /* check called with sensible 'k' values */

   assert(k_lb > 0);
   assert(k_ub >= k_lb);

   /* check called with sensible must_be_included, must_be_excluded values */

   assert(must_be_included != NULL || n_must_be_included == 0);
   assert(must_be_excluded != NULL || n_must_be_excluded == 0);

   /* check that if a CI cut then k_lb = k_ub = 1 */

   assert(!ci_cut || (k_lb == 1 && k_ub == 1));

   (*nGen) = 0;
   (*found_efficacious_ptr) = FALSE;


   /* create and initialise auxiliary IP data structure */

   /* allocate temporary memory for subscip elements */

   SCIP_CALL( SCIPallocBufferArray(scip, &auxipdata_family, psd->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxipdata_incluster, psd->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxipdata_clausecons, psd->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &auxipdata_overlapcons, psd->n) );

   for( i = 0 ; i < psd->n ; ++i )
   {
      auxipdata_incluster[i] = FALSE;
      SCIP_CALL( SCIPallocBufferArray(scip, &(auxipdata_family[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(auxipdata_clausecons[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(auxipdata_overlapcons[i]), psd->nParentSets[i]) );
      for( k = 0 ; k < psd->nParentSets[i] ; ++k )
      {
         auxipdata_family[i][k] = NULL;
         auxipdata_clausecons[i][k] = NULL;
         auxipdata_overlapcons[i][k] = NULL;
      }
   }

   /* allocate temporary memory for building clausal constraints */

   SCIP_CALL( SCIPallocBufferArray(scip, &clausevars, (psd->n) + 1) );

   /* create and initialise subscip */

   SCIP_CALL( SCIPcreate(&(auxipdata_subscip)) );


   SCIP_CALL( SCIPincludeDefaultPlugins(auxipdata_subscip) );

   SCIP_CALL( SCIPcreateProb(auxipdata_subscip, "DAG cluster separating MIP", NULL, NULL , NULL , NULL , NULL , NULL , NULL) );

   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetCharParam(auxipdata_subscip, "nodeselection/childsel", 'd') );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "limits/maxsol", 100000) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "limits/maxorigsol", 2000) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "nodeselection/dfs/stdpriority", 536870911) );
   SCIP_CALL( SCIPsetHeuristics(auxipdata_subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "lp/solvefreq", 1) );

   SCIP_CALL( SCIPsetRealParam(auxipdata_subscip, "limits/time", limits_time) );
   SCIP_CALL( SCIPsetRealParam(auxipdata_subscip, "limits/gap", limits_gap) );
   SCIP_CALL( SCIPsetRealParam(auxipdata_subscip, "limits/absgap", limits_absgap) );

   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/closecuts/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/cgmip/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/cmir/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/flowcover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/impliedbounds/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/intobj/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/mcf/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/oddcycle/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/rapidlearning/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/strongcg/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/zerohalf/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "separating/clique/freq", -1) );

   /* experimental */
   /*   SCIP_CALL(  SCIPsetIntParam(auxipdata_subscip, "separating/gomory/freq", -1)  );*/

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "heuristics/rins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata_subscip, "heuristics/crossover/freq", -1) );

   SCIP_CALL( SCIPsetBoolParam(auxipdata_subscip, "constraints/logicor/negatedclique", FALSE) );


   /* create subscip family variables for each main-problem family variable that is positive in the current solution */
   /* This solution will typically be the solution to the linear relaxation */
   /* use the same name for both */

   for( i = 0 ; i < psd->n ; ++i )
   {

      assert(solinfo->nposvars[i] > -1);
      for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
      {
         k = solinfo->posvars[i][ki];
         assert(k > -1 && k < psd->nParentSets[i]);
         /* essential not to consider the empty parent set */
         if( psd->nParents[i][k] == 0 )
            continue;

         /* val = SCIPgetSolVal(scip, sol, psd->PaVars[i][k]); */
         val = solinfo->lpsolvals[i][k];
         /*if ( !SCIPisPositive(scip,val) )
           printf("foo %f, %d, %d\n",val,i,ki);*/
         assert(SCIPisPositive(scip, val));
         SCIP_CALL( SCIPcreateVarBasic(auxipdata_subscip, &(auxipdata_family[i][k]), SCIPvarGetName(psd->PaVars[i][k]),
               0.0, 1.0, val, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(auxipdata_subscip, auxipdata_family[i][k]) );
      }
      
      /* create variables to identify clusters */
      /* would make it a dummy variable if no non-empty parent sets have positive value, but this oddly causes a slow down (why?) */

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "incluster#%d", i);
      SCIP_CALL( SCIPcreateVarBasic(auxipdata_subscip, &(auxipdata_incluster[i]), varname,
            0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(auxipdata_subscip, auxipdata_incluster[i]) );
      SCIP_CALL( SCIPchgVarBranchPriority(auxipdata_subscip, auxipdata_incluster[i], 10) );
   }

   /* create variable for lower bound */
   /* convenient to create it, even if k_lb=k_ub=1 */
   SCIP_CALL( SCIPcreateVarBasic(auxipdata_subscip, &(auxipdata_kvar), "kvar",
         k_lb, k_ub, 1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(auxipdata_subscip, auxipdata_kvar) );

   /* constraint that incumbent lies on the cutting plane */
   if( incumbent_cons )
   {
      best_sol = SCIPgetBestSol(scip);
      if( best_sol != NULL )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "incumbent");
         SCIP_CALL(SCIPcreateConsLinear(auxipdata_subscip, &(auxipdata_incumbentcons), consname, 0, NULL, NULL,
                                        -1, -1,
                                        TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
         for( i = 0 ; i < psd->n ; ++i )
         {
            SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_incumbentcons, auxipdata_incluster[i], -1) );
            for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
            {
               k = solinfo->posvars[i][ki];
               if( solinfo->lpsolvals[i][k] > 0.5 )
               {
                  if( psd->nParents[i][k] > 0 )
                     SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_incumbentcons, auxipdata_family[i][k], 1) );
                  break;
               }
            }
         }
         SCIP_CALL( SCIPaddCons(auxipdata_subscip, auxipdata_incumbentcons) );
         SCIP_CALL( SCIPreleaseCons(auxipdata_subscip, &(auxipdata_incumbentcons)) );
      }
   }

   /* if family[i][k]=1 then incluster[i]=1 */
   /* (ie in consistent notation with constraints below: If I(W->v)=1 then incluster[v] */
   /*  ~family[i][k]=1 + incluster[i] >= 1 */
   for( i = 0 ; i < psd->n ; ++i )
      for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
      {
         k = solinfo->posvars[i][ki];
         /* essential not to consider the empty parent set */
         if( psd->nParents[i][k] == 0 )
            continue;

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "clause#%d#%d", i, k);
         SCIP_CALL( SCIPgetNegatedVar(auxipdata_subscip, auxipdata_family[i][k], &(clausevars[0])) );
         clausevars[1] = auxipdata_incluster[i];
         SCIP_CALL(SCIPcreateConsLogicor(auxipdata_subscip, &(auxipdata_clausecons[i][k]), consname, 2, clausevars,
                                         TRUE, TRUE, TRUE, TRUE,
                                         TRUE, /*propagate*/
                                         FALSE, FALSE, FALSE, FALSE, FALSE));
         SCIP_CALL( SCIPaddCons(auxipdata_subscip, auxipdata_clausecons[i][k]) );
         SCIP_CALL( SCIPreleaseCons(auxipdata_subscip, &(auxipdata_clausecons[i][k])) );
      }



   /* k_ub*I(W->v) <= \sum_{u \in W} - k + k_ub */
   /* if I(W->v)=1 this becomes \sum_{u \in W} >= k */
   /* if I(W->v)=0 this becomes vacuous */
   /* just post as a normal linear constraint:
      -inf <= k_ub*I(W->v) - \sum_{u \in W} + k <= k_ub

      note: in the code below 'k' has a different meaning. It indexes
      parent sets
      ' u \in W' is represented by the binary variable incluster[psd->ParentSets[i][k][l]]

      if k_ub == 1 use an equivalent logicor representation

   */

   if( k_ub == 1 )
   {
      for( i = 0 ; i < psd->n ; ++i )
         for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
         {
            k = solinfo->posvars[i][ki];
            /* essential not to consider the empty parent set */
            if( psd->nParents[i][k] == 0 )
               continue;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "overlap#%d#%d", i, k);

            SCIP_CALL( SCIPgetNegatedVar(auxipdata_subscip, auxipdata_family[i][k], &(clausevars[0])) );
            nvars = 1;
            parent_set = psd->ParentSets[i][k];
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               /* clausevars[nvars++] = auxipdata_incluster[psd->ParentSets[i][k][l]]; */
               clausevars[nvars++] = auxipdata_incluster[parent_set[l]];
            }
            SCIP_CALL(SCIPcreateConsLogicor(auxipdata_subscip, &(auxipdata_clausecons[i][k]), consname, nvars, clausevars,
                                            TRUE, TRUE, TRUE, TRUE,
                                            TRUE, /* propagate*/
                                            FALSE, FALSE, FALSE, FALSE, FALSE));
            SCIP_CALL( SCIPaddCons(auxipdata_subscip, auxipdata_clausecons[i][k]) );
            SCIP_CALL( SCIPreleaseCons(auxipdata_subscip, &(auxipdata_clausecons[i][k])) );
         }
   }
   else
   {
      for( i = 0 ; i < psd->n ; ++i )
         for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
         {
            k = solinfo->posvars[i][ki];

            /* essential not to consider the empty parent set */
            if( psd->nParents[i][k] == 0 )
               continue;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "overlap#%d#%d", i, k);
            SCIP_CALL(SCIPcreateConsLinear(auxipdata_subscip, &(auxipdata_overlapcons[i][k]), consname, 0, NULL, NULL,
                                           -SCIPinfinity(scip), k_ub,
                                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
            SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_overlapcons[i][k], auxipdata_family[i][k], k_ub) );
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_overlapcons[i][k],
                     auxipdata_incluster[psd->ParentSets[i][k][l]], -1) );
            }
            SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_overlapcons[i][k], auxipdata_kvar, 1) );
            SCIP_CALL( SCIPaddCons(auxipdata_subscip, auxipdata_overlapcons[i][k]) );
            SCIP_CALL( SCIPreleaseCons(auxipdata_subscip, &(auxipdata_overlapcons[i][k])) );
         }

   }


   /* 2 <= |C|-k  <= inf : for the added cut to make sense */
   /* for k=1, this becomes 3 <= |C| <= inf */
   /* can ignore clusters of size two, since arrow variables already deal with them
      and no use for conditional independence cuts */
   if( k_ub == 1 )
      SCIP_CALL(SCIPcreateConsLinear(auxipdata_subscip, &(auxipdata_ck_cons), "ck_constraint", 0, NULL, NULL,
                                     2, SCIPinfinity(scip),
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   else
      SCIP_CALL(SCIPcreateConsLinear(auxipdata_subscip, &(auxipdata_ck_cons), "ck_constraint", 0, NULL, NULL,
                                     1, SCIPinfinity(scip),
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

   for( i = 0 ; i < psd->n ; ++i )
   {
      SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_ck_cons, auxipdata_incluster[i], 1) );
   }
   if( k_ub != 1 )
      SCIP_CALL( SCIPaddCoefLinear(auxipdata_subscip, auxipdata_ck_cons, auxipdata_kvar, -1) );
   SCIP_CALL( SCIPaddCons(auxipdata_subscip, auxipdata_ck_cons) );
   SCIP_CALL( SCIPreleaseCons(auxipdata_subscip, &(auxipdata_ck_cons)) );

   /* all constraints posted - free temporary memory */

   SCIPfreeBufferArray(scip, &clausevars);

   /* let I(u) denote u is in the cluster, then
      objective function is \sum_{v,W}I(W->v) - \sum_{u}I(u) + k
      If a feasible solution has a positive objective value,
      then a cutting plane has been found,
      so maximise and rule out non-positive solutions
   */

   SCIP_CALL_ABORT(SCIPsetObjsense(auxipdata_subscip, SCIP_OBJSENSE_MAXIMIZE));

   /* for cluster cuts rule out non-positive solutions using SCIPsetObjlimit
      for CI cuts -1 is enough */
   if( ci_cut )
      SCIP_CALL( SCIPsetObjlimit(auxipdata_subscip, -1) );
   else
      SCIP_CALL( SCIPsetObjlimit(auxipdata_subscip, 0) );

   /*SCIP_CALL( SCIPwriteOrigProblem(auxipdata_subscip,NULL,NULL,FALSE)  );
     SCIP_CALL( SCIPwriteParams(auxipdata_subscip,NULL,FALSE,TRUE)  );
    */

   for( i = 0; i < n_must_be_included; ++i )
   {
      SCIP_CALL( SCIPfixVar(auxipdata_subscip, auxipdata_incluster[must_be_included[i]], 1.0, &infeasible, &fixed) );
      assert(!infeasible && fixed);
   }
   for( i = 0; i < n_must_be_excluded; ++i )
   {
      SCIP_CALL( SCIPfixVar(auxipdata_subscip, auxipdata_incluster[must_be_excluded[i]], 0.0, &infeasible, &fixed) );
      assert(!infeasible && fixed);
   }

   if( ci_cut )
      SCIPdebugMessage("Looking for a conditional independence cluster cut using a subIP...\n");
   else
      SCIPdebugMessage("Looking for a cluster cut using a subIP...\n");

   SCIP_CALL( SCIPsolve(auxipdata_subscip) );

   status = SCIPgetStatus(auxipdata_subscip);
   /*SCIP_CALL( SCIPprintStatus(auxipdata_subscip,NULL) );*/

   if( status == SCIP_STATUS_INFEASIBLE )
   {
      SCIPdebugMessage("There is no cluster cut for this solution.\n");
      goto TERMINATE;
   }

      
   if( status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFORUNBD )
   {
      /* /\* print out LP relaxation that could not be separated *\/ */
      /*SCIP_CALL(  SCIPwriteMIP(scip,"foo",TRUE,TRUE)  );
      SCIP_CALL(  SCIPprintSol(scip,NULL,NULL,FALSE)  );
      exit(1);*/

      /*printf("infeasible.\n");
      SCIP_CALL(  SCIPprintSol(scip,NULL,NULL,FALSE)  );*/

      SCIPdebugMessage("could not find a cluster cut.\n");
      goto TERMINATE;
   }

   /* if there are feasible solutions but the best has objective value not better that
      0, then we have not found a cutting plane.
      This code snippet from Timo Berthold
   */
   nsols = SCIPgetNSols(auxipdata_subscip);
   if( nsols > 0 && SCIPisFeasLE(auxipdata_subscip, SCIPgetSolOrigObj(auxipdata_subscip, SCIPgetBestSol(auxipdata_subscip)), 0.0) )
   {
      /*printf("obj value too low.\n");
        SCIP_CALL(  SCIPprintSol(scip,NULL,NULL,FALSE)  );
      */
      SCIPdebugMessage("could not find a cluster cut: best cluster objective too low\n");
      goto TERMINATE;
   }

   if( status != SCIP_STATUS_SOLLIMIT && status != SCIP_STATUS_GAPLIMIT &&
      status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_NODELIMIT  && status != SCIP_STATUS_TIMELIMIT )
   {
      SCIPerrorMessage("Solution of subscip for DAG cluster separation returned with invalid status %d.\n", status);
      goto TERMINATE;
   }

   /* To get here a cutting plane must have been found */

   subscip_sols = SCIPgetSols(auxipdata_subscip);

   SCIP_CALL( SCIPallocBufferArray(scip, &incluster, psd->n) );

   for( s = 0; s <  nsols; ++s )
   {
      subscip_sol = subscip_sols[s];

      /*SCIP_CALL(  SCIPprintSol(auxipdata_subscip,subscip_sol,NULL,FALSE)  );*/

      excess = SCIPgetSolOrigObj(auxipdata_subscip, subscip_sol);

      /* only solutions which correspond to cuts should be considered */
      if( !SCIPisPositive(auxipdata_subscip,excess) )
         break;
      
#ifdef SCIP_DEBUG
      if( ci_cut )
         SCIPdebugMessage("found conditional independence cut for this cluster: ");
      else
         SCIPdebugMessage("found cut for this cluster: ");
      for( i = 0 ; i < psd->n ; ++i )
         if( SCIPisPositive(scip, SCIPgetSolVal(auxipdata_subscip, subscip_sol, auxipdata_incluster[i])) )
            SCIPdebugPrintf("%d,", i);
      SCIPdebugPrintf(" %f", excess);
      SCIPdebugPrintf("\n");
#endif

      kval = SCIPgetSolVal(auxipdata_subscip, subscip_sol, auxipdata_kvar);
      cluster_size = 0;
      for( i = 0 ; i < psd->n ; ++i )
         if( SCIPgetSolVal(auxipdata_subscip, subscip_sol, auxipdata_incluster[i]) > 0.5 )
         {
            incluster[i] = TRUE;
            cluster_size++;
         }
         else
            incluster[i] = FALSE;

      SCIP_CALL( AddClusterCut(scip, psd, conshdlr, sol, incluster, cluster_size, kval,
            addtopool, forcecuts, found_efficacious_ptr, ci_cut, matroidpaircuts, matroidpaircutslimit, matroidpairmaxcuts, excess, cutoff) );
      
      if( *cutoff )
      {
         SCIPdebugMessage("Adding cluster cut led to cutoff\n"); 
         SCIPfreeBufferArray(scip, &incluster);
         goto TERMINATE;
      }
      (*nGen)++;
   }

   SCIPfreeBufferArray(scip, &incluster);

   SCIPdebugMessage("added %d cluster cuts.\n", *nGen);

TERMINATE:
   /* deallocate in reverse order */
   for( i = (psd->n)-1 ; i > -1 ; --i )
   {
      for(k = 0; k < psd->nParentSets[i]; k++ )
      {
         if( auxipdata_family[i][k] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(auxipdata_subscip, &(auxipdata_family[i][k])) );
         }
      }
      SCIP_CALL( SCIPreleaseVar(auxipdata_subscip, &(auxipdata_incluster[i])) );
      SCIPfreeBufferArray(scip, &(auxipdata_overlapcons[i]));
      SCIPfreeBufferArray(scip, &(auxipdata_clausecons[i]));
      SCIPfreeBufferArray(scip, &(auxipdata_family[i]));
            
   }
   SCIPfreeBufferArray(scip, &auxipdata_overlapcons);
   SCIPfreeBufferArray(scip, &auxipdata_clausecons);
   SCIPfreeBufferArray(scip, &auxipdata_incluster);
   SCIPfreeBufferArray(scip, &auxipdata_family);

   SCIP_CALL( SCIPreleaseVar(auxipdata_subscip, &(auxipdata_kvar)) ); 
   
   if( auxipdata_subscip != NULL )
      SCIP_CALL( SCIPfree(&(auxipdata_subscip)) );

   return SCIP_OKAY;
}


