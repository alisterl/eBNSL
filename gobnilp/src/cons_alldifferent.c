
/**@file   cons_alldifferent.c
 * @brief  constraint handler for alldifferent constraints
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>

#include "cons_alldifferent.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "alldifferent"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY         -10 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         1 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
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

#define min(A,B) ((A) > (B) ? (B) : (A))
#define max(A,B) ((A) > (B) ? (A) : (B))


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for alldifferent constraints */
struct SCIP_ConsData
{
   int nvars;                  /**< number of variables in the constraint */
   SCIP_VAR** vars;            /**< variables in the constraint */
};

/** constraint handler data */
/* struct SCIP_ConshdlrData */
/* { */
/* }; */


/*
 * Local methods
 */


static
SCIP_Bool naiveprop(
   SCIP* scip,
   int n,
   int* lbs,
   int* ubs
   )
{
   int i;
   int j;
   int jj;
   int iu;
   int jl;

   SCIP_Bool result;

   int* outside;
   int* ius;
   int* jls;
   int* tmplbs;
   int* tmpubs;

   int previousub;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &outside, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ius, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &jls, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmplbs, n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpubs, n) );
   for( i = 0; i < n; i++ )
   {
      ius[i] = i;
      jls[i] = i;
      tmplbs[i] = lbs[i];
      tmpubs[i] = ubs[i];
   }
   
   SCIPsortDownIntInt(tmplbs,jls,n);
   SCIPsortDownIntInt(tmpubs,ius,n);

   SCIPfreeBufferArray(scip, &tmpubs);
   SCIPfreeBufferArray(scip, &tmplbs);

   result = TRUE;

   previousub = -1;
   for( iu = 0; iu < n && result; iu++ )
   {
      int lb;
      int ub;
      int ninside;
      int noutside;
      
      i = ius[iu];
      /*NB ubs change in this function */
      ub = ubs[i];

      /* if this ub unchanged from last iteration
         don't bother looking for a Hall set */
      if( ub == previousub )
         continue;
      else
         previousub = ub;

      
      ninside = 0;
      noutside = 0;
      lb = n;     /* dummy value */
      for( jl = 0; jl < n && result; jl++ )
      {
         int capacity;

         j = jls[jl];

         /* since lbs change, need to ensure that 'lb'
            is a lower bound for all variables seen so far */
         lb = min(lb,lbs[j]);

         if( ub < lb )
            continue;

         SCIPdebugMsg(scip, "lb=%d, ub=%d\n", lb, ub);

         /* is variable j in the interval? */
         /* note that variable j *and all earlier variables in this loop*
            have lb <= lbs[j] */
         if( ubs[j] <= ub )
            ninside++;
         else
            outside[noutside++] = j;
         
         capacity = ub - lb + 1;
         
         if( ninside > capacity )
            result = FALSE;
         else if (ninside == capacity )
         {
            /* remove, if possible, [lb,ub] from domains of those outside */
            /* these will be (1) {j: j = jls[jj], jj > jl} - 'later variables' */
            /* and (2) all those 'earlier variables' whose upper bound > ub */

            /* later variables */
            for(jj = jl+1; jj < n; jj++)
            {
               int j2;
               
               j2 = jls[jj];
               if( lb <= lbs[j2] && ub + 1 > lbs[j2])
                  lbs[j2] = ub + 1;
               if( ubs[j2] <= ub && lb - 1 < ubs[j2])
                  ubs[j2] = lb - 1;
            }

            /* earlier variables */
            for(jj = 0; jj < noutside; jj++)
            {
               int j2;
               
               j2 = outside[jj];
               if( lb <= lbs[j2] && ub + 1 > lbs[j2])
                  lbs[j2] = ub + 1;
               if( ubs[j2] <= ub && lb - 1 < ubs[j2])
                  ubs[j2] = lb - 1;
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &jls);
   SCIPfreeBufferArray(scip, &ius);
   SCIPfreeBufferArray(scip, &outside);

   return result;
}

/** separate by adding facets */
static
SCIP_RETCODE ADseparate(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_SOL*             sol,                /**< solution to be separated */
   int*                  nGen,               /**< output: pointer to store number of added rows */
   SCIP_Bool*            cutoff              /**< output: pointer to store whether we detected a cutoff */
   )
{
   char s[SCIP_MAXSTRLEN];
   int i;
   int j;
   SCIP_Real* vals;
   SCIP_VAR** tmpvars;
   SCIP_Real valsum;
   int lb;
   const int sum = nvars*(nvars-1)/2;
   
   assert( scip != NULL );
   assert( vars != NULL );
   assert( nGen != NULL );
   assert( cutoff != NULL );

   *cutoff = FALSE;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, nvars) );
   for(i = 0; i < nvars; i++)
      tmpvars[i] = vars[i];

   /* get values in this solution */
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, vals) );

   /* order variables according to solution value */
   SCIPsortRealPtr(vals, (void**)tmpvars, nvars);

   valsum = 0.0;
   lb = 0;
   for(i = 0; i < nvars; i++)
   {
      lb += i;
      valsum += vals[i];
      if( SCIPisFeasLT(scip, valsum, lb) )
      {
         SCIP_ROW *row;

         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "lb#%d", i);

         if( i < nvars / 2 )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, lb, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
            for( j = 0; j <= i; j++ )
               SCIP_CALL( SCIPaddVarToRow(scip, row, tmpvars[j], 1.0) );
         }
         else
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), sum - lb, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
            for( j = i+1; j < nvars; j++ )
               SCIP_CALL( SCIPaddVarToRow(scip, row, tmpvars[j], 1.0) );
         }
         SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
         SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
         SCIP_CALL( SCIPreleaseRow(scip, &row));
         ++(*nGen);
      }
   }
   
   SCIPfreeBufferArray(scip, &tmpvars);
   SCIPfreeBufferArray(scip, &vals);
   
   return SCIP_OKAY;
}


/** creates constraint data for alldifferent constraint */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   int                   nvars,              /**< number of variables in the alldifferent operation */
   SCIP_VAR**            vars                /**< variables in alldifferent operation */
   )
{

   assert(consdata != NULL);
   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->vars, vars, nvars) );
   (*consdata)->nvars = nvars;

   /* get transformed variables, if we are in the transformed problem */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVars(scip, (*consdata)->nvars, (*consdata)->vars, (*consdata)->vars) );
   }

   return SCIP_OKAY;
}

/** frees constraint data for alldifferent constraint */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->nvars);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** checks alldifferent constraint for feasibility of given solution: returns TRUE iff constraint is feasible */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int i;
   int j;
   
   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violated = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, consdata->nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, consdata->nvars, consdata->vars, vals) );
   
   for( i = 0; i < consdata->nvars && !(*violated); i++ )
   {
      for( j = i+1; j < consdata->nvars; j++ )
      {
         if( SCIPisEQ(scip, vals[i], vals[j]) )
         {
            SCIPdebugMsg(scip, "%d=%g, %d=%g\n",i,vals[i],j,vals[j]);
            *violated = TRUE;
            break;
         }
      }
   }

   SCIPfreeBufferArray(scip, &vals);

   SCIPdebugMsg(scip, "checked alldifferent constraint <%s> for feasibility of solution %p, violated = %d\n",
      SCIPconsGetName(cons), (void*)sol, *violated);

   
   return SCIP_OKAY;
}

/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a alldifferent constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdAlldifferent)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to alldifferent constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMsg(scip, "upgrading constraint <%s> to alldifferent constraint\n", SCIPconsGetName(cons));

      /* create the bin Alldifferent constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsAlldifferent(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyAlldifferent NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeAlldifferent NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitAlldifferent NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitAlldifferent NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreAlldifferent NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreAlldifferent NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolAlldifferent NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolAlldifferent NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteAlldifferent)
{  /*lint --e{715}*/

   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}



/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransAlldifferent NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpAlldifferent)
{  /*lint --e{715}*/


   char s[SCIP_MAXSTRLEN];
   int c;
   int nGen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR** vars;
      SCIP_ROW* row;
      int nvars;
      int sum;
      
      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMsg(scip, "adding initial rows for all different constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      nvars = consdata->nvars;
      vars = consdata->vars;
      sum = nvars*(nvars-1)/2;

      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sum#%d#%d", nvars, sum);
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, sum, sum, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, row, nvars, vars, 1.0) );
#ifdef SCIP_DEBUG
      SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
      SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row));
      ++nGen;
      
      /* cannot handle infeasible case here - just exit */
      if ( *infeasible )
         return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "added %d equations.\n", nGen);

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpAlldifferent)
{  /*lint --e{715}*/

   int nGen = 0;
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
      SCIP_Bool cutoff;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "separating LP solution for all different constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( ADseparate(scip, conshdlr, consdata->nvars, consdata->vars, NULL, &nGen, &cutoff) );
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   if (nGen > 0)
      *result = SCIP_SEPARATED;
   SCIPdebugMsg(scip, "separated %d cuts.\n", nGen);


   return SCIP_OKAY;
}



/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolAlldifferent)
{  /*lint --e{715}*/
   int nGen = 0;
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
      SCIP_Bool cutoff;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "separating solution for alldifferent constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( ADseparate(scip, conshdlr, consdata->nvars, consdata->vars, sol, &nGen, &cutoff) );
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   if (nGen > 0)
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpAlldifferent)
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
SCIP_DECL_CONSENFORELAX(consEnforelaxAlldifferent)
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
SCIP_DECL_CONSENFOPS(consEnfopsAlldifferent)
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
SCIP_DECL_CONSCHECK(consCheckAlldifferent)
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

         if( printreason )
         {
            int v;
            int w;
            SCIP_Real vx;
            SCIP_Real wx;
            SCIP_CONSDATA* consdata;
            SCIP_Bool ok;
            
            consdata = SCIPconsGetData(conss[i]);
            assert( consdata != NULL );

            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );

            ok = TRUE;
            for( v = 0; v < consdata->nvars && ok; ++v )
            {
               vx = SCIPgetSolVal(scip, sol, consdata->vars[v]);
               for( w = v+1; w < consdata->nvars; ++w )
               {
                  wx = SCIPgetSolVal(scip, sol, consdata->vars[w]);
                  if( SCIPisEQ(scip, vx, wx) )
                  {
                     ok = FALSE;
                     break;
                  }
               }
            }

            if( !ok )
            {
               SCIPinfoMessage(scip, NULL, ";\nviolation: <%s>=<%g> and <%s>=<%g>\n",
                  consdata->vars[v],vx,consdata->vars[w],wx);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropAlldifferent)
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
      SCIP_VAR** vars;
      int i;
      /* int j; */
      /* int k; */
      int nvars;
      /* int* count; */
      /* int** prefixvars; */
      int* ubs;
      /* int* tmp; */
      int* lbs;
      SCIP_Bool feasible;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "propagating alldifferent constraint <%s>.\n", SCIPconsGetName(cons));

      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );

      vars = consdata->vars;
      nvars = consdata->nvars;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &lbs, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ubs, nvars) ); 

      for (i = 0; i < nvars; ++i)
      {
         lbs[i] = SCIPconvertRealToInt(scip, SCIPfeasCeil(scip,SCIPvarGetLbLocal(vars[i]))); 
         ubs[i] = SCIPconvertRealToInt(scip, SCIPfeasFloor(scip,SCIPvarGetUbLocal(vars[i]))); 
      }

#ifdef SCIP_DEBUG
      for (i = 0; i < nvars; ++i)
         printf("before %d: [%d-%d]\n", i, lbs[i], ubs[i]);
#endif
      
      feasible = naiveprop(scip, nvars, lbs, ubs);

#ifdef SCIP_DEBUG
      if( feasible )
      {
         for (i = 0; i < nvars; ++i)
            printf("after %d: [%d-%d]\n", i, lbs[i], ubs[i]);
      }
      else
      {
         printf("infeasible\n");
      }
#endif

      
      if( !feasible )
      {
         SCIPdebugMessage(" -> node infeasible.\n");
         *result = SCIP_CUTOFF;
      }
      else
      {
         for (i = 0; i < nvars; ++i)
         {
            SCIP_Bool infeasible; 
            SCIP_Bool tightened; 

            SCIP_CALL( SCIPtightenVarLb(scip,vars[i],lbs[i],TRUE,&infeasible,&tightened) );
            if( infeasible )
            {
               SCIPdebugMessage(" -> node infeasible.\n");
               *result = SCIP_CUTOFF;
               break;
            }
            if( tightened ) 
            { 
               SCIPdebugMessage(" -> lower bound of <%s> increased to %d \n", SCIPvarGetName(vars[i]),lbs[i]); 
               *result = SCIP_REDUCEDDOM;
            }

            SCIP_CALL( SCIPtightenVarUb(scip,vars[i],ubs[i],TRUE,&infeasible,&tightened) );
            if( infeasible )
            {
               SCIPdebugMessage(" -> node infeasible.\n");
               *result = SCIP_CUTOFF;
               break;
            }
            if( tightened ) 
            { 
               SCIPdebugMessage(" -> upper bound of <%s> decreased to %d \n", SCIPvarGetName(vars[i]),ubs[i]); 
               *result = SCIP_REDUCEDDOM;
            }
         }
      }
      SCIPfreeBlockMemoryArray(scip, &ubs, nvars);
      SCIPfreeBlockMemoryArray(scip, &lbs, nvars);
   }
   return SCIP_OKAY;
}


/*       SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) ); */
/*       SCIP_CALL( SCIPallocBufferArray(scip, &tmp, nvars) ); */
      
/*       for (i = 0; i < nvars; ++i) */
/*       { */
/*          ubs[i] = SCIPconvertRealToInt(scip, SCIPfeasFloor(scip,SCIPvarGetUbLocal(vars[i]))); */
/*          tmp[i] = i; */
/*       } */
      
/*       SCIPsortIntInt(ubs,tmp,nvars); */

/*       for( i = 0; i < nvars; ++i ) */
/*       { */
/*          if( ubs[i] > i ) */
/*             break; */
/*       } */

/*       /\* variables from i onwards can have their lower bounds increased to i *\/ */
/*       for( j = i; j < nvars; ++j ) */
/*       { */
/*          SCIP_Bool infeasible; */
/*          SCIP_Bool tightened; */

/*          SCIP_CALL( SCIPtightenVarLb(scip,vars[tmp[j]],i,TRUE,&infeasible,&tightened) ); */

/*          if( infeasible ) */
/*          { */
/*             SCIPdebugMessage(" -> node infeasible.\n"); */
/*             *result = SCIP_CUTOFF; */
/*             break; */
/*          } */
/*          if( tightened ) */
/*          { */
/*             SCIPdebugMessage(" -> lower bound of <%s> increased to %d \n", SCIPvarGetName(vars[tmp[j]]),i); */
/*             *result = SCIP_REDUCEDDOM; */
/*          } */
/*       } */
/*       SCIPfreeBufferArray(scip, &ubs); */
/*       SCIPfreeBufferArray(scip, &tmp); */
/*       if( *result == SCIP_CUTOFF ) */
/*          return SCIP_OKAY; */
/*    } */
/*       return SCIP_OKAY; */
/* }       */
         
      
/*       /\* count[j] will be the number of vars with a local upper bound ub */
/*          such that ub <= j */

/*          if vars 2,5,6 have ubs <= 2 */
/*          then prefixvars[2] = [2,5,6]  */
/*          and count[2] = 3 */
/*       *\/ */
/*       SCIP_CALL( SCIPallocBufferArray(scip, &count, nvars) ); */
/*       SCIP_CALL( SCIPallocBufferArray(scip, &prefixvars, nvars) ); */
/*       for (j = 0; j < nvars; ++j) */
/*       { */
/*          SCIP_CALL( SCIPallocBufferArray(scip, &(prefixvars[j]), j+1) ); */
/*          count[j] = 0; */
/*       } */

/*       for (i = 0; i < nvars; ++i) */
/*       { */
/*          int ub; */

/*          ub = SCIPconvertRealToInt(scip, SCIPfeasFloor(scip,SCIPvarGetUbLocal(vars[i]))); */
/*          assert( ub >= 0 ); */
/*          assert( ub < nvars ); */
         
/*          for (j = ub; j < nvars; ++j) */
/*          { */
/*             assert( count[j] <= j+1 ); */
/*             if( count[j] == j+1 ) */
/*             { */
/*                /\* already full, that i needs to be added indicates */
/*                   infeasibility *\/ */
/* #ifdef SCIP_DEBUG */
/*                int jj; */
/*                for (jj = 0; jj < nvars; ++jj) */
/*                { */
/*                   printf("%d: ",jj); */
/*                   for (k = 0; k < count[jj]; ++k) */
/*                      printf("%d ",prefixvars[jj][k]); */
/*                   printf("\n"); */
/*                } */
/*                printf("But want to add %d to row %d\n", i, j); */
/* #endif */
/*                SCIPdebugMessage(" -> node infeasible.\n"); */
/*                *result = SCIP_CUTOFF; */
/*                goto TIDYUP; */
/*             } */
/*             prefixvars[j][count[j]++] = i; */
/*          } */
/*       }          */

/* #ifdef SCIP_DEBUG */
/*       for (j = 0; j < nvars; ++j) */
/*       { */
/*          printf("%d: ",j); */
/*          for (k = 0; k < count[j]; ++k) */
/*             printf("%d ",prefixvars[j][k]); */
/*          printf("\n"); */
/*       } */
/* #endif */
      
/*       for (j = 0; j < nvars; ++j) */
/*       { */
/*          if ( count[j] == j+1 ) */
/*          { */
/*             /\* bump up lower bound of all variables not in prefix *\/ */
/*             k = 0; */
/*             for (i = 0; i < nvars; ++i) */
/*             { */
/*                if( prefixvars[j][k] == i ) */
/*                { */
/*                   k++; */
/*                } */
/*                else */
/*                { */
/*                   SCIP_Bool infeasible; */
/*                   SCIP_Bool tightened; */
/*                   SCIP_CALL( SCIPtightenVarLb(scip,vars[i],j+1,TRUE,&infeasible,&tightened) ); */
/*                   if( infeasible ) */
/*                   { */
/*                      SCIPdebugMessage(" -> node infeasible.\n"); */
/*                      *result = SCIP_CUTOFF; */
/*                      goto TIDYUP; */
/*                   } */
/*                   if( tightened ) */
/*                   { */
/*                      SCIPdebugMessage(" -> lower bound of <%s> increased to %d \n", SCIPvarGetName(vars[i]),j+1); */
/*                      *result = SCIP_REDUCEDDOM; */
/*                   } */
/*                } */
/*             } */
/*          } */
/*       } */

/*    TIDYUP: */
/*       for (j = nvars-1; j >= 0; --j) */
/*          SCIPfreeBufferArray(scip, &(prefixvars[j])); */
/*       SCIPfreeBufferArray(scip, &prefixvars); */
/*       SCIPfreeBufferArray(scip, &count); */

/*       if( *result == SCIP_CUTOFF ) */
/*       { */
/* #ifdef SCIP_DEBUG */
/*          printf("Infeasible\n"); */
/*          for (i = 0; i < nvars; ++i) */
/*             printf("%d, %g -- %g\n",i, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])); */
/* #endif          */
/*          return SCIP_OKAY; */
/*       } */
/*    } */
   
/*    return SCIP_OKAY; */
/* } */



/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolAlldifferent NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropAlldifferent NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockAlldifferent)
{  /*lint --e{715}*/


   SCIP_CONSDATA* consdata;
   int i;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   for( i = 0; i < consdata->nvars; i++)
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[i], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveAlldifferent NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveAlldifferent NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableAlldifferent NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableAlldifferent NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsAlldifferent NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintAlldifferent NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyAlldifferent NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseAlldifferent NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsAlldifferent NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsAlldifferent NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsAlldifferent)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of alldifferent constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsAlldifferent NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for alldifferent constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrAlldifferent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create alldifferent constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpAlldifferent, consEnfopsAlldifferent, consCheckAlldifferent, consLockAlldifferent,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   /* SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyAlldifferent, consCopyAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveAlldifferent) ); */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteAlldifferent) );
   /* SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolAlldifferent) ); */
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpAlldifferent) ); 
   /* SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseAlldifferent) ); */
   /* SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolAlldifferent, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) ); */
   /* SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintAlldifferent) ); */
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropAlldifferent, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, 
         CONSHDLR_PROP_TIMING) ); 
   /* SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropAlldifferent) ); */
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpAlldifferent, consSepasolAlldifferent, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   /* SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransAlldifferent) ); */
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxAlldifferent) ); 



#ifdef LINCONSUPGD_PRIORITY
   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdAlldifferent, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }
#endif

   /* add alldifferent constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a alldifferent constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsAlldifferent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars,               /**< array with variables of constraint entries */
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsAlldifferent() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the alldifferent constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("alldifferent constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   SCIP_CALL( consdataCreate(scip, &consdata, nvars, vars) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a alldifferent constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicAlldifferent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in the constraint */
   SCIP_VAR**            vars                /**< array with variables of constraint entries */
   )
{
   SCIP_CALL( SCIPcreateConsAlldifferent(scip, cons, name, nvars, vars, 
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
