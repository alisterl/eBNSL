/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2018 James Cussens, Mark Bartlett        *
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

/**@file   pricer_family.c
 * @brief  family variable pricer
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>

#include "pricer_family.h"
#include "cons_dagcluster.h"
#include "probdata_bn.h"
#include "scip/pub_lp.h"
#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"

#define INITDUALINFOSIZE 1000

#define PRICER_NAME            "family"
#define PRICER_DESC            "pricer for new family variables"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */




/*
 * Data structures
 */

/* TODO: fill in the necessary variable pricer data */

/** variable pricer data */
struct SCIP_PricerData
{
   int n;                              /**< number of BN variables */
   SCIP_CONS** one_parent_set_conss;   /**< constraints stating that there is at most one parent set for any child */
   SCIP_CONS** arrow_conss;            /**< arrow_conss[n*i + j] is the cons where i<-j is an upper bound on sum of
                                          relevant family variables */
   SCIP_CONS*  dagcluster_cons;        /**< the single dagcluster constraint */
   int nspc_conss;                     /**< number of set packing constraints */
   CLUSTER_CUT** spc_conss;            /**< 'set packing' constraints */
   DUALINFO* dualinfo;                 /**< dualinfo (updated after solving each LP) */
};




/*
 * Local methods
 */

/* static */
/* SCIP_RETCODE */
/* storeDualInfo( */
/*    SCIP* scip, */
/*    int n; */
/*    DUALINFO** dualinfo to populate */
/*    ) */
/* { */
   
/* } */



/*
 * Callback methods of variable pricer
 */

/* TODO: Implement all necessary variable pricer methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for pricer plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRICERCOPY(pricerCopyFamily)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of family variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerCopyFamily NULL
#endif

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeFamily)
{  /*lint --e{715}*/

   SCIP_PRICERDATA* pricerdata;
   int n;
   
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   if( pricerdata != NULL && pricerdata->one_parent_set_conss != NULL)
   {
      n = pricerdata->n;

      /* free memory ( all should be non-NULL ) */

      SCIPfreeBlockMemoryArray(scip, &pricerdata->one_parent_set_conss, n);
      SCIPfreeBlockMemoryArray(scip, &pricerdata->arrow_conss, n*n);
      SCIPfreeBlockMemoryArray(scip, &pricerdata->spc_conss, pricerdata->nspc_conss);

      SCIPfreeBlockMemoryArray(scip, &pricerdata->dualinfo->childpenalties, n);
      SCIPfreeBlockMemoryArray(scip, &pricerdata->dualinfo->arrowpenalties, n*n);
      SCIPfreeBlockMemoryArray(scip, &pricerdata->dualinfo->nclusters, n);
      SCIPfreeBlockMemoryArray(scip, &pricerdata->dualinfo->clusters, n);
      /* no need to free pricerdata->dualinfo->clusters[i], this is freed after every pricing round */
      SCIPfreeBlockMemory(scip, &pricerdata->dualinfo);
      
      SCIPfreeBlockMemory(scip, &pricerdata);
   }
   return SCIP_OKAY;
}

/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
/* static */
/* SCIP_DECL_PRICERINITSOL(pricerInitsolFamily) */
/* {  /\*lint --e{715}*\/ */

/*    SCIP_PRICERDATA* pricerdata; */
/*    SCIP_CONS* cons; */
/*    int c; */

/*    assert(scip != NULL); */
/*    assert(pricer != NULL); */

/*    pricerdata = SCIPpricerGetData(pricer); */
/*    assert(pricerdata != NULL); */


/*    /\* get transformed set packing constraints  */
/*       and replace original constraints with the transformed ones */
/*    *\/ */
/*    for( c = 0; c < pricerdata->nspc_conss; ++c ) */
/*    { */
/*       cons = pricerdata->spc_conss[c]; */

/*       assert(cons != NULL); */
      
/*       /\* release original constraint *\/ */
/*       SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->spc_conss[c]) ); */

/*       /\* get transformed constraint *\/ */
/*       SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->spc_conss[c]) ); */

/*       /\* capture transformed constraint *\/ */
/*       SCIP_CALL( SCIPcaptureCons(scip, pricerdata->spc_conss[c]) ); */
/*    } */


/*    return SCIP_OKAY; */
/* } */



/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitFamily)
{  /*lint --e{715}*/

   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS* cons;
   int c;
   int i;
   int j;
   int n;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   n = pricerdata->n;
   
   /* get transformed one parent set constraints 
      and replace original constraints with the transformed ones
   */
   for( c = 0; c < n; ++c )
   {
      cons = pricerdata->one_parent_set_conss[c];

      assert(cons != NULL);

      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->one_parent_set_conss[c]) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->one_parent_set_conss[c]) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, pricerdata->one_parent_set_conss[c]) );
   }
   
   /* get transformed arrow  constraints 
      and replace original constraints with the transformed ones
   */
   for( i = 0; i < n; ++i )
      for( j = 0; j < n; ++j )
      {
         if(i == j)
            continue;

         c = n*i+j;
         cons = pricerdata->arrow_conss[c];
         
         assert(cons != NULL);
         
         /* release original constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->arrow_conss[c]) );
         
         /* get transformed constraint */
         SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->arrow_conss[c]) );
         
         /* capture transformed constraint */
         SCIP_CALL( SCIPcaptureCons(scip, pricerdata->arrow_conss[c]) );
      }

   /* get transformed set packing constraints 
      and replace original constraints with the transformed ones
   */
   for( c = 0; c < pricerdata->nspc_conss; ++c )
   {
      cons = (pricerdata->spc_conss[c])->cons;

      assert(cons != NULL);
      
      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &((pricerdata->spc_conss[c])->cons)) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &((pricerdata->spc_conss[c])->cons)) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, (pricerdata->spc_conss[c]->cons)) );
   }

   /* get transformed dagcluster constraint 
      and replace original constraint with the transformed one
   */
   
   cons = pricerdata->dagcluster_cons;

   assert(cons != NULL);
   
   /* release original constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->dagcluster_cons) );
   
   /* get transformed constraint */
   SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->dagcluster_cons) );

   /* capture transformed constraint */
   SCIP_CALL( SCIPcaptureCons(scip, pricerdata->dagcluster_cons) );
   
   return SCIP_OKAY;
}

/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolFamily)
{  /*lint --e{715}*/

   SCIP_PRICERDATA* pricerdata;
   int i;
   int j;
   int c;
   int n;
   
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   n = pricerdata->n;
   
   /* get release constraints */
   for( c = 0; c < n; ++c )
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->one_parent_set_conss[c])) );
   }

   /* get release constraints */
   for( i = 0; i < n; ++i )
      for( j = 0; j < n; ++j )
      {
         if( i == j )
            continue;

         c = n*i+j;

         /* release constraint */
         SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->arrow_conss[c])) );
      }
   
   /* get release constraints */
   for( c = 0; c < pricerdata->nspc_conss; ++c )
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &((pricerdata->spc_conss[c])->cons)) );
   }

   SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->dagcluster_cons)) );
   
   return SCIP_OKAY;
}



/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostFamily)
{

   int nclustercuts;                  /**< number of cluster cuts generated */
   CLUSTER_CUT** clustercuts;         /**< cluster cuts generated */
   /* int nspcuts;                  /\**< number of sp cuts generated *\/ */
   /* CLUSTER_CONS** spcuts;         /\**< sp cuts generated *\/ */
   int i;
   int j;


#ifdef SCIP_DEBUG
   SCIP_ROW** rows;
   SCIP_ROW* row;
   int nrows;
#endif
   
   SCIP_PRICERDATA* pricerdata;
   int n;
   DUALINFO* dualinfo;
   int c;

   int n_allclusters;
   int elt;
   
   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   /* grab old dualinfo from pricer data -
      it will get updated by this function
   */
   dualinfo = pricerdata->dualinfo;
   n = pricerdata->n;

   /* get child-specific dual penalties */ 
   for( c = 0; c < n; c++ )
   {
      
      dualinfo->childpenalties[c] = -SCIPgetDualsolSetppc(scip,pricerdata->one_parent_set_conss[c]);
      assert(!SCIPisPositive(scip,-dualinfo->childpenalties[c]));
      SCIPdebugMessage("Dual value for child %d is %f\n", c, -dualinfo->childpenalties[c]);
   }

   /* get arrow-specific dual penalties */ 
   for( i = 0; i < n; i++ )
   {
      for( j = 0; j < n; j++ )
      {
         if( i == j )
            continue;

         c = n*i + j;
         
         dualinfo->arrowpenalties[c] = -SCIPgetDualsolSetppc(scip,pricerdata->arrow_conss[c]);
         assert(!SCIPisPositive(scip,-dualinfo->arrowpenalties[c]));

         SCIPdebugMessage("Dual value for arrow %d<-%d is %f\n", i, j, -dualinfo->arrowpenalties[c]);
      }
   }


   /* get cluster cuts (but process after SPC constraints) */
   SCIP_CALL( DC_getclustercuts(scip, &nclustercuts, &clustercuts) );
   assert( nclustercuts == 0 || clustercuts != NULL );
   n_allclusters = nclustercuts + pricerdata->nspc_conss;
   
   /* initialise for cluster information */
   for(i = 0; i < n; i++)
   {
      /* allocated enough space even if variable i involved in all possible clusters */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(dualinfo->clusters[i]), n_allclusters) );
      dualinfo->nclusters[i] = 0;
   }
   
   /* get dual penalties for SPC constraints */
   for( c = 0; c < pricerdata->nspc_conss; c++ )
   {
      CLUSTER_CUT* cluster_cut;

      cluster_cut = pricerdata->spc_conss[c];
      
      assert(cluster_cut != NULL);
      assert(cluster_cut->nelts > 2);
      assert(cluster_cut->elts != NULL);

      /* a cons not a cut */
      assert(cluster_cut->cons != NULL);
      assert(cluster_cut->row == NULL);

      /* update dualsol */
      cluster_cut->dualsol = -SCIPgetDualsolSetppc(scip,cluster_cut->cons);
      assert(!SCIPisPositive(scip,-cluster_cut->dualsol));
      
#ifdef SCIP_DEBUG
      printf("Dual value for spc constraint (");
      for(i = 0; i < cluster_cut->nelts; i++)
         printf("%d ", cluster_cut->elts[i]);
      printf(") is %f\n", cluster_cut->dualsol);
#endif

      /* only consider active conss */
      if( !SCIPisZero(scip,cluster_cut->dualsol) )
      {
         for(i = 0; i < cluster_cut->nelts; i++)
         {
            elt = cluster_cut->elts[i];
            dualinfo->clusters[elt][dualinfo->nclusters[elt]++] = cluster_cut;
         }
      }
   }

   /* now process cuts */
   for(i = 0; i < nclustercuts; i++)
   {
      CLUSTER_CUT* cluster_cut;

      cluster_cut = clustercuts[i];
      
      assert(cluster_cut != NULL);
      assert(cluster_cut->elts != NULL);
      assert(cluster_cut->nelts > 2);

      /* a cut not a cons */
      assert(cluster_cut->cons == NULL);
      assert(cluster_cut->row != NULL);

#ifdef SCIP_DEBUG
      printf("cut %d k=%d row %p ", i, cluster_cut->k, (void*) cluster_cut->row);
      for(j = 0; j < cluster_cut->nelts; j++)
         printf("%d,",cluster_cut->elts[j]);
#endif

      /* put in correct dualsol */
      
      if( cluster_cut->row == NULL )
      {
         SCIPdebugMessage("NULL row\n");
         cluster_cut->dualsol = 0.0;

      }
      else
      {
         cluster_cut->dualsol = -SCIProwGetDualsol(cluster_cut->row);
         SCIPdebugMessage(" : dual val = %f\n", cluster_cut->dualsol);
      }

      /* only consider active cuts */
      if( !SCIPisZero(scip,cluster_cut->dualsol) )
      {
         for(j = 0; j < cluster_cut->nelts; j++)
         {
            elt = cluster_cut->elts[j];
            dualinfo->clusters[elt][dualinfo->nclusters[elt]++] = cluster_cut;
         }
      }
   }
   
#ifdef SCIP_DEBUG
   rows = SCIPgetLPRows(scip);
   nrows = SCIPgetNLPRows(scip);
   printf("START dual values for %d rows\n", nrows);
   for(i = 0; i < nrows; i++)
   {
      SCIP_Real val;

      row = rows[i];
      val = SCIProwGetDualsol(row);
      printf("Row %d, dualval=%f:\n",i, val);
      SCIPprintRow(scip,row,NULL);
   }
   printf("END dual values for %d rows\n", nrows);
#endif

   /* shrink to fit */
   for( i = 0; i < n; i++)
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(dualinfo->clusters[i]), n_allclusters, dualinfo->nclusters[i]) );

#ifdef SCIP_DEBUG
   printf("Dualinfo:\n");
   for( i = 0; i < n; i++)
      printf("Child penalty %d: %f\n", i, dualinfo->childpenalties[i]);

   for( i = 0; i < n; i++)
      for( j = 0; j < n; j++)
      {
         if(i == j)
            continue;

         c = n*i+j;

         printf("Arrow penalty %d<-%d: %f\n", i, j,  dualinfo->arrowpenalties[c]);
      }

   for( i = 0; i < n; i++)
   {
      printf("Cluster penalties for %d:\n", i);
      for( j = 0; j < dualinfo->nclusters[i]; j++)
      {
         CLUSTER_CUT* cluster_cut = dualinfo->clusters[i][j];
         printf("Dual value for k = %d cluster ", cluster_cut->k);
         for(c = 0; c < cluster_cut->nelts; c++)
            printf("%d ", cluster_cut->elts[c]);
         printf(" is %f\n", cluster_cut->dualsol);
      }
   }
#endif
   
   /* todo, send dualinfo to scoring to find new variables */
   
   
   (*result) = SCIP_SUCCESS;


   for( i = 0; i < n; i++)
      SCIPfreeBlockMemoryArray(scip, &(dualinfo->clusters[i]), dualinfo->nclusters[i]);

   
   return SCIP_OKAY;
}


/** Farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasFamily)
{  /*lint --e{715}*/
   
   SCIPwarningMessage(scip, "Current LP is infeasible, but Farkas pricing was not implemented\n");
   SCIPABORT();

   return SCIP_OKAY; /*lint !e527*/
}






/*
 * variable pricer specific interface methods
 */

/** creates the family variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerFamily(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create family variable pricer data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );

   pricerdata->n = 0;
   pricerdata->one_parent_set_conss = NULL;
   pricerdata->arrow_conss = NULL;
   pricerdata->dagcluster_cons = NULL;
   pricerdata->nspc_conss = 0;
   pricerdata->spc_conss = NULL;
   pricerdata->dualinfo = NULL;

   /* include variable pricer */
   /* use SCIPincludePricerBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostFamily, pricerFarkasFamily, pricerdata) );
   assert(pricer != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeFamily) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitFamily) );
   /* SCIP_CALL( SCIPsetPricerInitsol(scip, pricer, pricerInitsolFamily) ); */
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolFamily) );


   /* add family variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerFamilyActivate(
   SCIP* scip,                           /**< SCIP data structure */
   int n,                                /**< number of BN variables */
   SCIP_CONS** one_parent_set_conss,     /**< constraints stating that there is at most one parent set for any child */
   SCIP_CONS** arrow_conss,              /**< arrow_conss[n*i + j] is the cons where i<-j is an upper bound on sum of
                                            relevant family variables */
   SCIP_CONS*  dagcluster_cons,          /**< the single dagcluster constraint */
   int nspc_conss,                       /**< number of set packing constraints */
   CLUSTER_CUT** spc_conss               /**< 'set packing' constraints */ 
   )
{

   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   
   assert(scip != NULL);
   assert(n > 0);
   assert(one_parent_set_conss != NULL);
   assert(arrow_conss != NULL);
   assert(dagcluster_cons != NULL);
   assert(nspc_conss == 0 || spc_conss != NULL );
   
   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* copy arrays like in Binpacking */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->one_parent_set_conss, one_parent_set_conss, n) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->arrow_conss, arrow_conss, n*n) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->spc_conss, spc_conss, nspc_conss) );

   /* create space for dual info */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata->dualinfo) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->dualinfo->childpenalties, n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->dualinfo->arrowpenalties, n*n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->dualinfo->nclusters, n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->dualinfo->clusters, n) );

   pricerdata->n = n;
   pricerdata->nspc_conss = nspc_conss;
   pricerdata->dagcluster_cons = dagcluster_cons;

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );
   
   return SCIP_OKAY;
}
