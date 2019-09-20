
/**@file   branch_dualarrow.c
 * @brief  dualarrow branching rule
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>

#include "branch_dualarrow.h"
#include "probdata_bn.h"
#include "scip/cons_setppc.h"


#define BRANCHRULE_NAME            "dualarrow"
#define BRANCHRULE_DESC            "branching rule template"
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
   int n;                        /** number of BN variables in problem */
   SCIP_CONS** arrow_conss;      /**< arrow_conss[n*i + j] is the cons where i<-j is an upper bound on sum of
                                    relevant family variables */
   SCIP_VAR** arrow_vars;        /**< arrow_vars[i][j] is the indicator variable for i<-j, such variables also accessible via
                                    the hashmap in `psd` (or NULL if this array not used) */
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
SCIP_DECL_BRANCHCOPY(branchCopyDualarrow)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of dualarrow branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyDualarrow NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeDualarrow)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;
   int n;
   
   assert(scip != NULL);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);

   if( branchruledata != NULL )
   {
      n = branchruledata->n;
      if( branchruledata->arrow_conss != NULL)
         SCIPfreeBlockMemoryArray(scip, &branchruledata->arrow_conss, n*n);
      if( branchruledata->arrow_vars != NULL)
      SCIPfreeBlockMemoryArray(scip, &branchruledata->arrow_vars, n*n);

      SCIPfreeBlockMemory(scip, &branchruledata);
      
      SCIPbranchruleSetData(branchrule, NULL);
   }

   return SCIP_OKAY;
}



/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitDualarrow)
{  /*lint --e{715}*/

   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_PROBDATA* probdata;
   
   SCIP_CONS* cons;
   SCIP_VAR* var;
   int c;
   int i;
   int j;
   int n;

   assert(scip != NULL);
   assert(branchrule != NULL);

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);


   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   /* if these missing cannot use this rule */
   if( probdata->arrow_conss == NULL ||  probdata->arrow_vars == NULL )
   {
      branchruledata->n = 0;
      branchruledata->arrow_conss = NULL;
      branchruledata->arrow_vars = NULL;
      return SCIP_OKAY;
   }
   
   n = probdata->psd->n;
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &branchruledata->arrow_conss, probdata->arrow_conss, n*n) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &branchruledata->arrow_vars, probdata->arrow_vars, n*n) );
   branchruledata->n = n;


   /* get transformed arrow  constraints 
      and replace original constraints with the transformed ones
   */

   for( i = 0; i < n; ++i )
      for( j = 0; j < n; ++j )
      {
         if(i == j)
            continue;

         c = n*i+j;
         cons = branchruledata->arrow_conss[c];
         var = branchruledata->arrow_vars[c];
         assert(cons != NULL);

         
         /* get transformed constraint */
         SCIP_CALL( SCIPgetTransformedCons(scip, cons, &branchruledata->arrow_conss[c]) );

         /* get transformed variable */
         SCIP_CALL( SCIPgetTransformedVar(scip, var, &branchruledata->arrow_vars[c]) );  
      }

   return SCIP_OKAY;
}



/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BRANCHEXIT(branchExitDualarrow)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of dualarrow branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitDualarrow NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolDualarrow)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of dualarrow branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolDualarrow NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolDualarrow)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of dualarrow branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolDualarrow NULL
#endif


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpDualarrow)
{  /*lint --e{715}*/


   int i;
   int j;
   int n;
   int c;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_VAR* bestvar = NULL;

   SCIP_Real bestdualval = 0.0;
   SCIP_Real dualval;

   SCIP_BRANCHRULEDATA* branchruledata;

   *result = SCIP_DIDNOTRUN;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   n = branchruledata->n;

   /* just skip if not set up */
   if( n == 0 )
      return SCIP_OKAY;
   
   for( i = 0; i < n; ++i )
      for( j = 0; j < n; ++j )
      {
         if(i == j)
            continue;

         c = n*i+j;
         var = branchruledata->arrow_vars[c];
         assert(var != NULL);


         if( SCIPisFeasIntegral(scip, SCIPvarGetLPSol(var)) )
            continue;
         
         /* if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var) ) ) */
         /*    continue; */
         
         cons = branchruledata->arrow_conss[c];
         assert(cons != NULL);

         dualval = SCIPgetDualsolSetppc(scip,cons);
         if( dualval < bestdualval )
         {
            bestdualval = dualval;
            bestvar = branchruledata->arrow_vars[c];
         }
      }

   if( bestvar == NULL)
   {
      SCIPdebugMsg(scip, "No variable with non-zero dual value found\n");
      return SCIP_OKAY;
   }
   else
      SCIPdebugMsg(scip, "Branching on variable <%s> with associated dual value of <%g>\n", SCIPvarGetName(bestvar), bestdualval);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, bestvar, NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}



/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextDualarrow)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of dualarrow branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextDualarrow NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 0
static
SCIP_DECL_BRANCHEXECPS(branchExecpsDualarrow)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of dualarrow branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsDualarrow NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the dualarrow branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleDualarrow(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create dualarrow branching rule data */
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
   /* SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyDualarrow) ); */
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeDualarrow) );
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitDualarrow) );
   /* SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitDualarrow) ); */
   /* SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolDualarrow) ); */
   /* SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolDualarrow) ); */
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpDualarrow) );
   /* SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextDualarrow) ); */
   /* SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsDualarrow) ); */


   /* add dualarrow branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
