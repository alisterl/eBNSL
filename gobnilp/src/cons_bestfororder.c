/**@file   cons_bestfororder.c
 * @brief  constraint handler for bestfororder constraints
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>

#include "cons_bestfororder.h"
#include <string.h>


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "bestfororder"
#define CONSHDLR_DESC          "best DAG for a given order constraint handler"
#define CONSHDLR_ENFOPRIORITY         -10 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */




/* TODO: (optional) enable linear or nonlinear constraint upgrading */
#if 0
#include "scip/cons_linear.h"
#include "scip/cons_nonlinear.h"
#define LINCONSUPGD_PRIORITY          0 /**< priority of the constraint handler for upgrading of linear constraints */
#define NONLINCONSUPGD_PRIORITY       0 /**< priority of the constraint handler for upgrading of nonlinear constraints */
#endif


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for bestfororder constraints */
struct SCIP_ConsData
{
   ParentSetData* psd;              /** parent set data */
   SCIP_VAR** posvars;              /** position variables */
};

/* /\** constraint handler data *\/ */
/* struct SCIP_ConshdlrData */
/* { */
/* }; */


/*
 * Local methods
 */

/** creates constraint data for bestfororder constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   ParentSetData*        psd,                /** parent set data */
   SCIP_VAR**            posvars             /** position variables */
   )
{

   assert(consdata != NULL);
   assert(psd != NULL);
   assert(posvars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* just store the pointers */
   (*consdata)->psd = psd;
   (*consdata)->posvars = posvars;

   return SCIP_OKAY;
}

static
SCIP_Bool contains(
   SCIP_Bool* seen,
   int* paset,
   int n
   )
{
   int i;
   
   for( i = 0; i < n; i++ )
   {
      if( !seen[paset[i]] )
         return FALSE;
   }
   return TRUE;
}

/** checks bestfororder constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int k;
   int n;

   SCIP_Real* posvals;
   int* tmp;
   SCIP_Bool* seen;
   int pos;
   
   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->psd != NULL);
   assert(consdata->psd->nParentSets != NULL);
   assert(consdata->psd->nParents != NULL);
   assert(consdata->psd->ParentSets != NULL);
   assert(consdata->psd->PaVars != NULL);

   
   *violated = FALSE;
   n = consdata->psd->n;

   SCIP_CALL( SCIPallocBufferArray(scip, &posvals, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmp, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &seen, n) );
   for( i = 0; i < n; i++)
   {
      tmp[i] = i;
      seen[i] = FALSE;
   }
   
   /* get values in this solution */
   SCIP_CALL( SCIPgetSolVals(scip, sol, n, consdata->posvars, posvals) );

   /* order variables according to solution value */
   SCIPsortRealInt(posvals, tmp, n);

   for(pos = 0; pos < n && !(*violated); pos++)
   {
      /* variable i is in position pos */
      i = tmp[pos];
      assert( i > -1 );
      assert( i < n );

      assert(consdata->psd->nParents[i] != NULL);
      assert(consdata->psd->ParentSets[i] != NULL);
      assert(consdata->psd->PaVars[i] != NULL);
      
      /* find best=first parent set for i such that all parents
         already seen
      */

      for(k = 0; k < consdata->psd->nParentSets[i] && !(*violated); k++)
      {
         assert(consdata->psd->ParentSets[i][k] != NULL);
         
         /* find out whether all parents in kth parent set come earlier in order */
         if( contains(seen, consdata->psd->ParentSets[i][k], consdata->psd->nParents[i][k]) )
         {
            /* just check that best parent set for order is set to 1,
               Don't bother to check that others set to 0, that is done elsewhere */
            if( !SCIPisEQ(scip, SCIPgetSolVal(scip, sol, consdata->psd->PaVars[i][k]), 1.0) )
            {
               *violated = TRUE;
            }
            seen[i] = TRUE;
            break;
         }
      }
   }
            
   SCIPfreeBufferArray(scip, &seen);
   SCIPfreeBufferArray(scip, &tmp);
   SCIPfreeBufferArray(scip, &posvals);

   SCIPdebugMsg(scip, "checked bestfororder constraint <%s> for feasibility of solution %p, violated = %d\n",
      SCIPconsGetName(cons), (void*)sol, *violated);

   return SCIP_OKAY;
   
}

               
      
/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a bestfororder constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdBestfororder)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to bestfororder constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMsg(scip, "upgrading constraint <%s> to bestfororder constraint\n", SCIPconsGetName(cons));

      /* create the bin Bestfororder constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsBestfororder(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyBestfororder NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeBestfororder NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitBestfororder NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitBestfororder NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreBestfororder NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreBestfororder NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolBestfororder NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolBestfororder NULL
#endif


/** frees specific constraint data */
#if 0
static
SCIP_DECL_CONSDELETE(consDeleteBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeleteBestfororder NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransBestfororder NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpBestfororder NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpBestfororder NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolBestfororder NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpBestfororder)
{  /*lint --e{715}*/

   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && *result == SCIP_FEASIBLE; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
      }
   }


   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for relaxation solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxBestfororder)
{  /*lint --e{715}*/

   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && *result == SCIP_FEASIBLE; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
      }
   }


   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsBestfororder)
{  /*lint --e{715}*/

   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && *result == SCIP_FEASIBLE; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
      }
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckBestfororder)
{  /*lint --e{715}*/

   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;

   /* method is called only for integral solutions, because the enforcing priority is negative */
   for( i = 0; i < nconss && *result == SCIP_FEASIBLE; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropBestfororder)
{  /*lint --e{715}*/

   int c;
   
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      int i;
      int n;
      int* lbs;
      int* ubs;
      int* lbsort;
      int* ubsort;
      SCIP_Bool proplbs;
      SCIP_Bool propubs;
      SCIP_Bool* earlier;

      SCIP_Bool found;
      int val;
      int j;
      int k;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "propagating bestforder constraint <%s>.\n", SCIPconsGetName(cons));

      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->posvars != NULL );
      assert( consdata->psd != NULL );

      n = consdata->psd->n;

      SCIP_CALL( SCIPallocBufferArray(scip, &lbs, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ubs, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lbsort, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ubsort, n) );
      SCIP_CALL( SCIPallocBufferArray(scip, &earlier, n) );
      
      proplbs = FALSE;
      propubs = FALSE;
      for( i = 0; i < n; i++ )
      {
         lbs[i] = SCIPconvertRealToInt(scip, SCIPfeasCeil(scip,SCIPvarGetLbLocal(consdata->posvars[i])));
         ubs[i] = SCIPconvertRealToInt(scip, SCIPfeasFloor(scip,SCIPvarGetUbLocal(consdata->posvars[i])));
         lbsort[i] = i;
         ubsort[i] = i;

         if( lbs[i] == n-1 )
            proplbs = TRUE;
         if( ubs[i] == 0 )
            propubs = TRUE;
      }

      if( propubs )
      {
         SCIPsortIntInt(ubs,ubsort,n);
         for(i = 0; i < n; i++ )
            earlier[i] = FALSE;

         /* ubs is array of upper bounds in non-decreasing order */
         /* ubsort[j] is the BN variable whose upper bound is ubs[j] */
         for(j = 0; j < n; j++ )
         {
            if( ubs[j] <= j )
            {
               
               i = ubsort[j];
               earlier[i] = TRUE;

               /* can select the parent set for variable i,
                it is the best=first one whose parents are all
                earlier in the ordering, fix all other parent sets to 0 */

               found = FALSE;
               for(k = 0; k < consdata->psd->nParentSets[i]; k++)
               {
                  /* find out whether all parents in kth parent set come earlier in order */
                  if( !found && contains(earlier, consdata->psd->ParentSets[i][k], consdata->psd->nParents[i][k]) )
                  {
                     SCIP_CALL( SCIPtightenVarLb(scip,consdata->psd->PaVars[i][k],1,TRUE,&infeasible,&tightened) );
                     found = TRUE;
                     val = 1;
                  }
                  else
                  {
                     SCIP_CALL( SCIPtightenVarUb(scip,consdata->psd->PaVars[i][k],0,TRUE,&infeasible,&tightened) );
                     val = 0;
                  }
                  
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible (tried for fix <%s> to <%d>, using upper bounds).\n",
                        SCIPvarGetName(consdata->psd->PaVars[i][k]),val);
                     *result = SCIP_CUTOFF;
                     goto TIDYUP;
                  }
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> Fixed variable <%s> to %d \n",
                        SCIPvarGetName(consdata->psd->PaVars[i][k]), val);
                     *result = SCIP_REDUCEDDOM;
                  }
               }
            }
            else
            {
               /* jth smallest ub was not small enough, no more propagations */
               break;
            }
         }
      }

      if( proplbs )
      {
         SCIPsortDownIntInt(lbs,lbsort,n);
         for(i = 0; i < n; i++ )
            earlier[i] = TRUE;

         /* lbs is array of lower bounds in non-increasing order */
         /* lbsort[j] is the BN variable whose lower bound is lbs[j] */
         for(j = 0; j < n; j++ )
         {
            if( lbs[j] >= (n-1)-j )
            {
               i = lbsort[j];
               earlier[i] = FALSE;
               /* can select the parent set for variable i,
                it is the best=first one whose parents are all
                earlier in the ordering, fix all other parent sets to 0 */

               found = FALSE;
               for(k = 0; k < consdata->psd->nParentSets[i]; k++)
               {
                  /* find out whether all parents in kth parent set come earlier in order */
                  if( !found && contains(earlier, consdata->psd->ParentSets[i][k], consdata->psd->nParents[i][k]) )
                  {
                     SCIP_CALL( SCIPtightenVarLb(scip,consdata->psd->PaVars[i][k],1,TRUE,&infeasible,&tightened) );
                     found = TRUE;
                     val = 1;
                  }
                  else
                  {
                     SCIP_CALL( SCIPtightenVarUb(scip,consdata->psd->PaVars[i][k],0,TRUE,&infeasible,&tightened) );
                     val = 0;
                  }
                  
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible (tried for fix <%s> to %d, using lower bounds).\n",
                        SCIPvarGetName(consdata->psd->PaVars[i][k]),val);
                     *result = SCIP_CUTOFF;
                     goto TIDYUP;
                  }
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> Fixed variable <%s> to %d \n",
                        SCIPvarGetName(consdata->psd->PaVars[i][k]), val);
                     *result = SCIP_REDUCEDDOM;
                  }
               }
            }
            else
            {
               /* jth biggest lb was not big enough, no more propagations */
               break;
            }
         }
      }

   TIDYUP:
      SCIPfreeBufferArray(scip, &earlier);
      SCIPfreeBufferArray(scip, &ubsort);
      SCIPfreeBufferArray(scip, &lbsort);
      SCIPfreeBufferArray(scip, &ubs);
      SCIPfreeBufferArray(scip, &lbs);
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }   
         
   return SCIP_OKAY;
}



/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolBestfororder NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropBestfororder NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockBestfororder)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;
   int i;
   int k;
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);


   /* lock everything in both directions */
   for( i = 0; i < consdata->psd->n; i++)
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->posvars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      for( k = 0; k < consdata->psd->nParentSets[i]; k++ )
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->psd->PaVars[i][k], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveBestfororder NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveBestfororder NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableBestfororder NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableBestfororder NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsBestfororder NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintBestfororder NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyBestfororder NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseBestfororder NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsBestfororder NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsBestfororder NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsBestfororder)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of bestfororder constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsBestfororder NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for bestfororder constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrBestfororder(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create bestfororder constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */

   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpBestfororder, consEnfopsBestfororder, consCheckBestfororder, consLockBestfororder,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   /* SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyBestfororder, consCopyBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolBestfororder, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) ); */
   /* SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintBestfororder) ); */
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropBestfororder, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, 
         CONSHDLR_PROP_TIMING) ); 
   /* SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropBestfororder) ); */
   /* SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpBestfororder, consSepasolBestfororder, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) ); */
   /* SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransBestfororder) ); */
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxBestfororder) ); 

#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdBestfororder, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add bestfororder constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a bestfororder constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBestfororder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,                /**< parent set data structure */
   SCIP_VAR**            posvars,            /**< array with position variables */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsBestfororder() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the bestfororder constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("bestfororder constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   SCIP_CALL( consdataCreate(scip, &consdata, psd, posvars) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a bestfororder constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicBestfororder(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,                /**< parent set data structure */
   SCIP_VAR**            posvars            /**< array with position variables */
   )
{
   SCIP_CALL( SCIPcreateConsBestfororder(scip, cons, name, psd, posvars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
