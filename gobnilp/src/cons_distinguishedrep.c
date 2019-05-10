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

/**@file   cons_distinguishedrep.c
 * @brief  constraint handler for Markov equivalence class distinguished representative constraints
 * @author James Cussens
 *
 *
 *  Doxygen documentation for the locally defined structs @c SCIP_ConsData and @c SCIP_ConshdlrData is unfortunately not available
 *  since structs with the same name are defined in other constraint handlers and Doxygen cannot handle the name clash. You will have to consult
 *  the source code to view this documentation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>
#include "cons_distinguishedrep.h"
#include "utils.h"


/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "distinguishedrep"
#define CONSHDLR_DESC          "constraint handler for distinguished Markov equivalence class representative constraints"
#define CONSHDLR_ENFOPRIORITY         -9010 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY         -90000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            1  /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */


#define MAXINITCUTS 0

/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for distinguishedrep constraints */
struct SCIP_ConsData
{
   ParentSetData* psd;
};


/*
 * Local methods
 */



/** Creates the data for a constraint.
 *  @param scip The SCIP instance to which the constraint belongs.
 *  @param consdata The location to store the new constraint data.
 *  @param psd The parent set data on which the constraint is based.
 *  @return SCIP_OKAY if successful, or an appropriate error otherwise.
 */
static
SCIP_RETCODE createConsData(
   SCIP* scip,
   SCIP_CONSDATA** consdata,
   ParentSetData* psd
   )
{

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( PS_copyParentSetData(scip, psd, &((*consdata)->psd)) );

   return SCIP_OKAY;
}

/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a distinguishedrep constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdDistinguishedrep)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to distinguishedrep constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to distinguishedrep constraint\n", SCIPconsGetName(cons));

      /* create the bin Distinguishedrep constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsDistinguishedrep(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyDistinguishedrep NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_CONSFREE(consFreeDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consFreeDistinguishedrep NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitDistinguishedrep NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitDistinguishedrep NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreDistinguishedrep NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreDistinguishedrep NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolDistinguishedrep NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolDistinguishedrep NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteDistinguishedrep)
{

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPdebugMessage("deleting distinguished rep constraint <%s>.\n", SCIPconsGetName(cons));

   SCIP_CALL( PS_deallocateParentSetData(scip, &((*consdata)->psd), FALSE) );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;

 
}



/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransDistinguishedrep NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
#if 0
static
SCIP_DECL_CONSTRANS(consInitlpDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpDistinguishedrep NULL
#endif

/** separation method of constraint handler for LP solutions */
#if 0
static
SCIP_DECL_CONSSEPALP(consSepalpDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpDistinguishedrep NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolDistinguishedrep NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDistinguishedrep)
{    

   int c;
   int nGen = 0;

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
      ParentSetData* psd;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("enforcing lp solution for  distinguished representative constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );
      psd = consdata->psd;

      if( !is_dr_feasible(scip, psd, NULL, FALSE, NULL) )
      {
         SCIP_ROW *row;
         SCIP_Bool infeasible;
         int i;
         int k;

         /* uncommenting following two lines causes a big slowdown on insurance_10000_1_3.scores */
         /* *result = SCIP_INFEASIBLE; */
	 /* return SCIP_OKAY; */

         
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "todo", -SCIPinfinity(scip), (psd->n)-1, FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
         for( i = 0; i < psd->n; i++ ) 
            for( k = 0; k < psd->nParentSets[i]; k++ )
               if( SCIPisGT(scip, SCIPgetSolVal(scip, NULL, psd->PaVars[i][k]), 0.5) )
               {
                  SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[i][k], 1.0) );
                  break;
               }
         SCIP_CALL( SCIPflushRowExtensions(scip, row) );
         SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
         SCIPdebug(SCIP_CALL( SCIPprintRow(scip, row, NULL) ));
         SCIP_CALL( SCIPreleaseRow(scip, &row));
         ++nGen;

         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

      }
      if (nGen > 0)
      {
	 *result = SCIP_SEPARATED;
	 return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all distinguished representative constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDistinguishedrep)
{  
 int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      ParentSetData* psd;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("enforcing pseudo solution for distinguished representative constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );
      psd = consdata->psd;

      if( !is_dr_feasible(scip, psd, NULL, FALSE, NULL) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all distinguished representative constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDistinguishedrep)
{  
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( conss != NULL );
   assert( result != NULL );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_CONS* cons;
      ParentSetData* psd;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("checking distinguished representative constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );
      psd = consdata->psd;

      if( !is_dr_feasible(scip, psd, sol, FALSE, NULL) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all distinguished representative constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropDistinguishedrep)
{  /*lint --e{715}*/

   int c;
   int nGen = 0;

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
      ParentSetData* psd;
      int i;
      int j;
      int k;
      int kk;
      int l_lo;
      int l_hi;
      int n;
      int child;
      SCIP_Bool children_ok;
      SCIP_Bool ok;
      SCIP_Bool hi_found;
      SCIP_VAR* arrow_var;

      int lo;
      int hi;
      SCIP_Bool ilo;
      int klo = -1;
      int khi = -1;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("propagating distinguished representative constraint <%s>.\n", SCIPconsGetName(cons));

      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );

      psd = consdata->psd;
      n = psd->n;
      
      for( i = 0; i < n; ++i)
      {
         /* find which parent set for i fixed to one, if any */
         for( k = 0; k < psd->nParentSets[i]; k++ )
            if ( SCIPvarGetLbLocal(psd->PaVars[i][k]) > 0.5 )
               break;

         if( k == psd->nParentSets[i] )
            /* no parent set fixed for vertex i, so no propagations possible */
            continue;

         /* look for a lo and hi (lo < hi) where lo = i or hi = i, such that an edge lo<-hi can be reversed */
         for( j = 0; j < n; ++j )
         {
            if( i == j )
               continue;

            if( i < j )
            {
               lo = i;
               hi = j;
               ilo = TRUE;
               klo = k;
            }
            else
            {
               lo = j;
               hi = i;
               ilo = FALSE;
               khi = k;
            }
            
            arrow_var = get_arrow(psd,lo,hi);
            if( arrow_var == NULL || SCIPvarGetUbLocal(arrow_var) < 0.5 )
               /* not possible for there to be a lo<-hi arrow to reverse, so continue */
               continue;

            children_ok = TRUE;
            for( child = 0; child < n; ++child)
            {
               if( child == lo || child == hi )
                  continue;
               
               arrow_var = get_arrow(psd,child,hi);
               if( arrow_var != NULL && SCIPvarGetUbLocal(arrow_var) > 0.5 )
               {
                  /* hi may be a parent of child, lo has to be also */
                  arrow_var = get_arrow(psd,child,lo);
                  if( arrow_var == NULL || SCIPvarGetLbLocal(arrow_var) < 0.5 )
                  {
                     /* child may not be a a child of lo, so continue */
                     children_ok = FALSE;
                     break;
                  }
               }
            }

            if( !children_ok )
               /* any covered arcs will not be reversible since hi may have children that lo does not
                  so continue search
               */
               continue;

            /* want lo to have hi's parents plus hi */
            for( kk = 0; kk < psd->nParentSets[j]; kk++ )
            {
               
               if( ilo )
                  khi = kk;
               else
                  klo = kk;
               
               if( psd->nParents[lo][klo] != psd->nParents[hi][khi] + 1 )
                     continue;

               l_hi = 0;
               ok = TRUE;
               hi_found = FALSE;

               assert( klo < psd->nParentSets[lo] );
               assert( khi < psd->nParentSets[hi] );
               for( l_lo = 0; l_lo < psd->nParents[lo][klo]; l_lo++ )
               {
                  if( psd->ParentSets[lo][klo][l_lo] == hi )
                  {
                     hi_found = TRUE;
                  }
                  else
                  {
                     assert( l_hi <= psd->nParents[hi][khi] );

                     if( l_hi == psd->nParents[hi][khi] || psd->ParentSets[lo][klo][l_lo] != psd->ParentSets[hi][khi][l_hi++] )
                     {
                        ok = FALSE;
                        break;
                     }
                  }
               }

               if( ok && hi_found )
               {
                  SCIP_Bool infeasible;
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip,psd->PaVars[j][kk],0,TRUE,&infeasible,&tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible.\n");
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  if( tightened )
                  {
                     SCIPdebugMessage("Ruling out: %s due to covered arc\n", SCIPvarGetName(psd->PaVars[j][kk]));
                     ++nGen;
                  }
                  break;
               }
            }
         }
      }
   }
   if( nGen > 0 )
      *result = SCIP_REDUCEDDOM;
   SCIPdebugMessage("propagated %d domains.\n", nGen);

   return SCIP_OKAY;
}



/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolDistinguishedrep NULL
#endif


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropDistinguishedrep NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDistinguishedrep)
{ 
   SCIP_CONSDATA* consdata;
   ParentSetData* psd;
   
   int i;
   int k;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMessage("Locking distinguished representative constraint <%s>.\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->psd != NULL );
   psd = consdata->psd;

   /* just lock everything in both directions */
   for( i = 0; i < psd->n; i++ )
      for( k = 0; k < psd->nParentSets[i]; k++ )
         SCIP_CALL( SCIPaddVarLocks(scip, psd->PaVars[i][k], nlockspos + nlocksneg, nlockspos + nlocksneg) );

   
   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveDistinguishedrep NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveDistinguishedrep NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableDistinguishedrep NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableDistinguishedrep NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsDistinguishedrep NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintDistinguishedrep NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyDistinguishedrep NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseDistinguishedrep NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsDistinguishedrep NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsDistinguishedrep NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsDistinguishedrep)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of distinguishedrep constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsDistinguishedrep NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for distinguishedrep constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDistinguishedrep(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create distinguishedrep constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */

   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpDistinguishedrep, consEnfopsDistinguishedrep, consCheckDistinguishedrep, consLockDistinguishedrep,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   /* SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyDistinguishedrep, consCopyDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveDistinguishedrep) ); */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteDistinguishedrep) );
   /* SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpDistinguishedrep) );  */
   /* SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolDistinguishedrep, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) ); */
   /* SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintDistinguishedrep) ); */
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropDistinguishedrep, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
      CONSHDLR_PROP_TIMING) );
   /* SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropDistinguishedrep) ); */
   /* SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpDistinguishedrep, consSepasolDistinguishedrep, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) ); */
   /* SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransDistinguishedrep) ); */

   return SCIP_OKAY;
}

/** creates and captures a distinguishedrep constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */

SCIP_RETCODE SCIPcreateConsDistinguishedrep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsDistinguishedrep() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the distinguishedrep constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("distinguishedrep constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* Initialise constraint data */
   SCIP_CALL( createConsData(scip, &consdata, psd) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a distinguishedrep constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicDistinguishedrep(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd
   )
{
   SCIP_CALL( SCIPcreateConsDistinguishedrep(scip, cons, name, psd,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
