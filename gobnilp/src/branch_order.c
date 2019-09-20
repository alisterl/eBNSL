
/**@file   branch_order.c
 * @brief  order branching rule
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>

#include "branch_order.h"
#include "probdata_bn.h"


#define BRANCHRULE_NAME            "order"
#define BRANCHRULE_DESC            "order using position variables branching rule"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
   SCIP_VAR** posvars;
   int n;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BRANCHCOPY(branchCopyOrder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of order branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyOrder NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeOrder)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;
   
   assert(scip != NULL);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);

   SCIPfreeBlockMemory(scip, &branchruledata);
   
   SCIPbranchruleSetData(branchrule, NULL);
   
   return SCIP_OKAY;
}



/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitOrder)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   /* just copy the pointer */
   branchruledata->posvars = probdata->posvars;
   branchruledata->n = probdata->psd->n;
   
   return SCIP_OKAY;
}



/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BRANCHEXIT(branchExitOrder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of order branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitOrder NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolOrder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of order branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolOrder NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolOrder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of order branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolOrder NULL
#endif


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpOrder)
{  /*lint --e{715}*/


   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** posvars;
   int lowestlb;
   int i;
   int n;
   int besti;
   SCIP_Real bestlpsolval;
   int lb;
   SCIP_Longint nodenum;
   
   *result = SCIP_DIDNOTRUN;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->posvars != NULL);
   
   posvars = branchruledata->posvars;
   n = branchruledata->n;

   lowestlb = n;
   besti = -1;
   bestlpsolval = 0.0;
   for( i = 0; i < n; i++ )
   {
      SCIP_VAR* var;
      SCIP_Real lpsolval;
      
      var = posvars[i];
      lb = SCIPconvertRealToInt(scip, SCIPfeasCeil(scip,SCIPvarGetLbLocal(var)));
      lpsolval = SCIPvarGetLPSol(var);

      assert(lb > -1);
      assert(lb < n );
      
      if( lb <= lowestlb && lpsolval > bestlpsolval )
      {
         /* cannot branch on fixed variables */
         if( lb < SCIPconvertRealToInt(scip, SCIPfeasFloor(scip,SCIPvarGetUbLocal(var))) )
         {
            lowestlb = lb;
            bestlpsolval = lpsolval;
            besti = i;
         }
      }
   }


   /* get current node number */
   nodenum = SCIPgetNNodes(scip);

   if( besti != -1 )
   {
      SCIPdebugMsg(scip, "Branching on variable <%s> at node <%d> using value <%d> (fractional LP)\n",
         SCIPvarGetName(posvars[besti]), nodenum, lowestlb);
      *result = SCIP_BRANCHED;
      SCIP_CALL( SCIPbranchVarVal(scip, posvars[besti], lowestlb, NULL, NULL, NULL) );
   }
   else
   {
      SCIPdebugMsg(scip, "Could not branch on a position variable at node <%d> (fractional LP)\n", nodenum);
   }

   return SCIP_OKAY;
}



/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextOrder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of order branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextOrder NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsOrder)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** posvars;
   int lowestlb;
   int i;
   int n;
   int besti;
   SCIP_Real bestlpsolval;
   int lb;
   SCIP_Longint nodenum;
   
   *result = SCIP_DIDNOTRUN;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->posvars != NULL);
   
   posvars = branchruledata->posvars;
   n = branchruledata->n;

   lowestlb = n;
   besti = -1;
   bestlpsolval = 0.0;
   for( i = 0; i < n; i++ )
   {
      SCIP_VAR* var;
      SCIP_Real lpsolval;

      var = posvars[i];
      lb = SCIPconvertRealToInt(scip, SCIPfeasCeil(scip,SCIPvarGetLbLocal(var)));
      lpsolval = SCIPvarGetLPSol(var);

      assert(lb > -1);
      assert(lb < n );
      
      if( lb < lowestlb && lpsolval > bestlpsolval )
      {
         /* cannot branch on fixed variables */
         if( lb < SCIPconvertRealToInt(scip, SCIPfeasFloor(scip,SCIPvarGetUbLocal(var))) )
         {
            lowestlb = lb;
            bestlpsolval = lpsolval;
            besti = i;
         }
      }
   }

   /* get current node number */
   nodenum = SCIPgetNNodes(scip);

   if( besti != -1 )
   {
      SCIPdebugMsg(scip, "Branching on variable <%s> at node <%d> using value <%d> (pseudo solution) \n",
         SCIPvarGetName(posvars[besti]), nodenum, lowestlb);
      *result = SCIP_BRANCHED;
      SCIP_CALL( SCIPbranchVarVal(scip, posvars[besti], lowestlb, NULL, NULL, NULL) );
   }
   else
   {
      SCIPdebugMsg(scip, "Could not branch on a position variable at node <%d> (pseudo solution)\n",nodenum);
#ifdef SCIP_DEBUG
      for( i = 0; i < n; i++ )
      {
         SCIP_VAR* var;

         var = posvars[i];
         SCIPdebugMsg(scip, "<%s>, lb=<%g>, ub=<%g>\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      }
#endif
   }


   return SCIP_OKAY;
}



/*
 * branching rule specific interface methods
 */

/** creates the order branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleOrder(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create order branching rule data */
   branchruledata = NULL;
   /* TODO: (optional) create branching rule specific data here */

   SCIP_CALL( SCIPallocBlockMemory(scip, &branchruledata) );

   branchrule = NULL;

   /* include branching rule */
   /* use SCIPincludeBranchruleBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   /* SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyOrder) ); */
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeOrder) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitOrder) ); 
   /* SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitOrder) ); */
   /* SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolOrder) ); */
   /* SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolOrder) ); */
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpOrder) ); 
   /* SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextOrder) ); */
   SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsOrder) ); 


   /* add order branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
