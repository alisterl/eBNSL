
/**@file   branch_ancestral.c
 * @brief  ancestral branching rule
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "branch_ancestral.h"
#include "probdata_bn.h"


#define BRANCHRULE_NAME            "ancestral"
#define BRANCHRULE_DESC            "ancestral branching"
#define BRANCHRULE_PRIORITY        0
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

#define ABS(x)        ((x) >= 0 ? (x) : -(x))

/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
      ParentSetData* psd;
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */
#if 0
static
int getancestral(
   SCIP_Bool* ancestral,
   int* selected,
   ParentSetData* psd
   )
{

   SCIP_Bool made_progress;
   SCIP_Bool all_ancestral;
   
   int i;
   int k;
   int l;

   int ancestralsetsize = 0;
   
   assert( ancestral != NULL);
   assert( psd != NULL);
   
   /* find any selected parent sets and mark vertices with selected empty parent sets as ancestral */
   made_progress = FALSE;
   for( i = 0; i < psd->n; ++i )
   {
      ancestral[i] = FALSE;
      selected[i] = -1;
      for( k = 0 ; k < psd->nParentSets[i]; ++k )
      {
         if( SCIPvarGetLbLocal(psd->PaVars[i][k]) > 0.5 )
         {
            if( psd->nParents[i][k] == 0 )
            {
               ancestral[i] = TRUE;
               ancestralsetsize++;
               made_progress = TRUE;
            }
            selected[i] = k;
            break;
         }
      }
   }
      
   while( made_progress )
   {
      made_progress = FALSE;
      for( i = 0; i < psd->n; ++i )
      {
         
         k = selected[i];
         if( k == -1 || ancestral[i] )
            continue;

         all_ancestral = TRUE;
         for(l = 0; l < psd->nParents[i][k]; ++l )
         {
            if( !ancestral[psd->ParentSets[i][k][l]] )
            {
               all_ancestral = FALSE;
               break;
            }
         }

         if( all_ancestral )
         {
            ancestral[i] = TRUE;
            ancestralsetsize++;
            made_progress = TRUE;
         }
      }
   }
   return ancestralsetsize;
}
#endif   

/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BRANCHCOPY(branchCopyAncestral)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ancestral branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyAncestral NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeAncestral)
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
SCIP_DECL_BRANCHINIT(branchInitAncestral)
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
   branchruledata->psd = probdata->psd;
   
   return SCIP_OKAY;
}



/** deinitialization method of branching rule (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_BRANCHEXIT(branchExitAncestral)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ancestral branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitAncestral NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_BRANCHINITSOL(branchInitsolAncestral)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ancestral branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolAncestral NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolAncestral)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ancestral branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolAncestral NULL
#endif


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpAncestral)
{  /*lint --e{715}*/

   ParentSetData* psd;
   SCIP_Bool* ancestral;                /* will represent the (locally) maximal ancestral set */
   int* selected;                       /* if parent set for child i already chosen to be k then selected[i]=k, else -1 */
   SCIP_Bool all_ancestral;
   
   int i;
   int k;
   int l;

   SCIP_Real bestscore;
   SCIP_VAR* bestvar = NULL;
   
   SCIP_BRANCHRULEDATA* branchruledata;

   *result = SCIP_DIDNOTRUN;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);
   assert(branchruledata->psd != NULL);
   
   psd = branchruledata->psd;

   SCIP_CALL( SCIPallocBufferArray(scip, &ancestral, psd->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selected, psd->n) );


   for( i = 0; i < psd->n; ++i )
   {
      if( selected[i] != -1 )
         continue;
         
      for( k = 0 ; k < psd->nParentSets[i]; ++k )
      {
         SCIP_Real lpsolval;
         SCIP_Real score;
         SCIP_VAR* var = psd->PaVars[i][k];

         if( SCIPvarGetStatus(SCIPvarGetTransVar(var)) == SCIP_VARSTATUS_MULTAGGR )
            continue;
         
         /* ignore variables fixed to 0 */
         if( SCIPvarGetUbLocal(var) < 0.5 )
            continue;

         /* cannot be fixed to 1 */
         assert( SCIPvarGetLbLocal(var) < 0.5 );
         
         all_ancestral = TRUE;
         for(l = 0; l < psd->nParents[i][k]; ++l )
         {
            if( !ancestral[psd->ParentSets[i][k][l]] )
            {
               all_ancestral = FALSE;
               break;
            }
         }
         if( !all_ancestral )
            continue;

         /* go for most infeasible */
         lpsolval = SCIPvarGetLPSol(var);
         if( lpsolval > 0.5 )
            score = 1.0 - lpsolval;
         else
            score = lpsolval;

         if( score > bestscore )
         {
            bestscore = score;
            bestvar = var;
         }
      }
   }

   if( bestvar == NULL)
   {
      SCIPdebugMsg(scip, "No 'ancestral' fractional variable found\n");
      return SCIP_OKAY;
   }
   else
      SCIPdebugMsg(scip, "Branching on variable <%s> with fractionality score of <%g>\n", SCIPvarGetName(bestvar), bestscore);

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, bestvar, NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}



/** branching execution method for external candidates */
#if 0
static
SCIP_DECL_BRANCHEXECEXT(branchExecextAncestral)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ancestral branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextAncestral NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#if 0
static
SCIP_DECL_BRANCHEXECPS(branchExecpsAncestral)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of ancestral branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsAncestral NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the ancestral branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleAncestral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create ancestral branching rule data */
   branchruledata = NULL;

   /* create branching rule specific data here */

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
   /* SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyAncestral) ); */
   SCIP_CALL( SCIPsetBranchruleFree(scip, branchrule, branchFreeAncestral) ); 
   SCIP_CALL( SCIPsetBranchruleInit(scip, branchrule, branchInitAncestral) );
   /* SCIP_CALL( SCIPsetBranchruleExit(scip, branchrule, branchExitAncestral) ); */
   /* SCIP_CALL( SCIPsetBranchruleInitsol(scip, branchrule, branchInitsolAncestral) ); */
   /* SCIP_CALL( SCIPsetBranchruleExitsol(scip, branchrule, branchExitsolAncestral) ); */
   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpAncestral) ); 
   /* SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextAncestral) ); */
   /* SCIP_CALL( SCIPsetBranchruleExecPs(scip, branchrule, branchExecpsAncestral) ); */


   /* add ancestral branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
