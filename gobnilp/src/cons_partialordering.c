/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* uncomment for debug output: */
/* #define SCIP_DEBUG */

/**@file   cons_partialordering.c
 * @brief  example constraint handler for partial ordering constraints
 * @author Marc Pfetsch
 * @author James Cussens
 *
 * We handle the following system of partial constraints:
 * - \f$ x_{ij} + x_{ji} \leq 1 \f$         (symmetry inequalities - added initially)
 * \f$ x_{ij} + x_{jk} - x_{ik} \leq 1 \f$  (triangle inequalities)
 *
 * The partial order must be consistent with the arrow variables which are also in the constraint
 * If minimal=TRUE, the partial order must be the minimal one consistent with the arrows (ie the ancestor relation)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cons_partialordering.h>
#include "utils.h"
#include "parent_set_data.h"

#include <assert.h>
#include <string.h>


/* constraint handler properties */
#define CONSHDLR_NAME          "partialordering"
#define CONSHDLR_DESC          "partial ordering constraint handler"
#define CONSHDLR_SEPAPRIORITY       100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY      -200 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY     -100 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            10 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP


/** constraint data for partial ordering constraints */
struct SCIP_ConsData
{
   int                   n;                  /**< number of elements */
   SCIP_VAR***           vars;               /**< partial order variables */
   ParentSetData*        psd;                /**< contains arrow and parent set variables */
   SCIP_Bool             minimal;            /**< whether the partial order must minimal while being consistent with arrows */
};


/** separate symmetry equations and triangle inequalities */
static
SCIP_RETCODE PartialOrderingSeparate(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of variables */
   SCIP_SOL*             sol,                /**< solution to be separated */
   int*                  nGen,               /**< output: pointer to store number of added rows */
   SCIP_Bool*            cutoff              /**< output: pointer to store whether we detected a cutoff */
   )
{
   char s[SCIP_MAXSTRLEN];
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( vars != NULL );
   assert( nGen != NULL );
   assert( cutoff != NULL );

   *cutoff = FALSE;
   for (i = 0; i < n && ! (*cutoff); ++i)
   {
      for (j = 0; j < n && ! (*cutoff); ++j)
      {
	 SCIP_Real valIJ;
	 if (j == i)
	    continue;

	 valIJ = SCIPgetSolVal(scip, sol, vars[i][j]);

	 /* if symmetry inequalities are violated - should not be the case, if they are added in the beginning */
	 if ( ! SCIPisFeasLE(scip, valIJ + SCIPgetSolVal(scip, sol, vars[j][i]), 1.0) )
	 {
	    SCIP_ROW *row;

	    (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sym#%d#%d", i, j);

	    SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
	    SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	    SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][i], 1.0) );
	    SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	    SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	    SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
	    SCIP_CALL( SCIPreleaseRow(scip, &row));
	    ++(*nGen);

            if ( *cutoff )
               break;
	 }

	 /* check triangle inequalities */
	 for (k = 0; k < n; ++k)
	 {
	    SCIP_Real sum = 0.0;
	    if (k == i || k == j)
	       continue;

	    sum = valIJ + SCIPgetSolVal(scip, sol, vars[j][k]) - SCIPgetSolVal(scip, sol, vars[i][k]);

	    /* if sum - 1.0 > 0, i.e., the cut is violated */
	    if ( SCIPisEfficacious(scip, sum - 1.0) )
	    {
	       SCIP_ROW *row;

	       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "triangle#%d#%d#%d", i, j, k);

	       SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
	       SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][k], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][k], -1.0) );
	       SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	       SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	       SCIP_CALL( SCIPaddRow(scip, row, FALSE, cutoff) );
	       SCIP_CALL( SCIPreleaseRow(scip, &row));
	       ++(*nGen);

               if ( *cutoff )
                  break;
	    }
	 }
      }
   }

   return SCIP_OKAY;
}


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyPartialOrdering)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrPartialOrdering(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeletePartialOrdering)
{  /*lint --e{715}*/
   int i;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( consdata != NULL);
   assert( *consdata != NULL);
   assert( (*consdata)->vars != NULL );

   SCIPdebugMsg(scip, "deleting partial ordering constraint <%s>.\n", SCIPconsGetName(cons));

   n = (*consdata)->n;
   SCIP_CALL( PS_deallocateParentSetData(scip, &((*consdata)->psd), FALSE) );

   for (i = 0; i < n; ++i)
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->vars[i]), n); /*lint !e866*/
   SCIPfreeBlockMemoryArray(scip, &((*consdata)->vars), n);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransPartialOrdering)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_CONSDATA* sourcedata;
   int i;
   int j;
   int k;
   int n;
   char s[SCIP_MAXSTRLEN];

   SCIP_Bool minimal;
   ParentSetData* targetpsd;

   SCIP_VAR* arrow_var;
   SCIP_VAR* new_arrow_var;
   
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   SCIPdebugMsg(scip, "transforming partial ordering constraint <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   n = sourcedata->n;
   consdata->n = n;

   minimal = sourcedata->n;
   consdata->minimal = minimal;

   /* copy parent set data, but does not create new variables */
   SCIP_CALL( PS_copyParentSetData(scip, sourcedata->psd, &targetpsd) );
   /* have to reinitiaise hashtable for edges and arrows since new variables
      will be inserted using put_arrow (see below) */
   SCIP_CALL( hashtablefreeArrow(scip, targetpsd) );
   SCIP_CALL( hashtableCreateArrow(scip, targetpsd) );

   consdata->psd = targetpsd;

   /* transform variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, n) );
   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->vars[i]), n) ); /*lint !e866*/
      for (j = 0; j < n; ++j)
      {
	 if (j != i)
	 {
	    assert( sourcedata->vars[i][j] != NULL );
	    SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->vars[i][j], &(consdata->vars[i][j])) );

            arrow_var = get_arrow(sourcedata->psd,i,j);
            if( arrow_var != NULL )
            {
               SCIP_CALL( SCIPgetTransformedVar(scip, arrow_var, &new_arrow_var) );
               put_arrow(scip, consdata->psd, i, j, new_arrow_var);
            }
         }
         for( k = 0; k < sourcedata->psd->nParentSets[i]; k++ )            
            SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->psd->PaVars[i][k], &(consdata->psd->PaVars[i][k])) );
      }
   }

   /* create constraint */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));

   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpPartialOrdering)
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
      int i, j, n;
      SCIP_VAR*** vars;
      ParentSetData*  psd;
      SCIP_VAR* arrow_var;
      
      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMsg(scip, "adding initial rows for partial ordering constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      n = consdata->n;
      vars = consdata->vars;
      psd = consdata->psd;

      /* add symmetry equation */
      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
            SCIP_ROW* row;
            
            if( j == i )
               continue;
            
            if( j > i )
            {
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sym#%d#%d", i, j);
               SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 1.0, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][i], 1.0) );
               SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
               SCIPdebug( SCIProwPrint(row, NULL) );
#endif
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &row));
               ++nGen;
               
               /* cannot handle infeasible case here - just exit */
               if ( *infeasible )
                  return SCIP_OKAY;
            }

            /* parents are ancestors */

            arrow_var = get_arrow(psd,i,j);
            if( arrow_var == NULL )
               continue;
            
            (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "paisanc#%d#%d", i, j);
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
            SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], -1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, row, arrow_var, 1.0) );
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );
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
      }
   }
   SCIPdebugMsg(scip, "added %d equations.\n", nGen);

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpPartialOrdering)
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
      SCIPdebugMsg(scip, "separating LP solution for partial ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( PartialOrderingSeparate(scip, conshdlr, consdata->n, consdata->vars, NULL, &nGen, &cutoff) );
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
SCIP_DECL_CONSSEPASOL(consSepasolPartialOrdering)
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
      SCIPdebugMsg(scip, "separating solution for partial ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( PartialOrderingSeparate(scip, conshdlr, consdata->n, consdata->vars, sol, &nGen, &cutoff) );
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
SCIP_DECL_CONSENFOLP(consEnfolpPartialOrdering)
{  /*lint --e{715}*/
   char s[SCIP_MAXSTRLEN];
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
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;
      ParentSetData* psd;
      SCIP_Bool minimal;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "enforcing lp solution for partial ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      n = consdata->n;
      vars = consdata->vars;
      psd = consdata->psd;
      minimal = consdata->minimal;
      assert( vars != NULL );

      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
	    SCIP_Real valIJ;
            SCIP_VAR* arrow_i_j;

	    if (j == i)
	       continue;

	    valIJ = SCIPgetSolVal(scip, NULL, vars[i][j]);
            arrow_i_j = get_arrow(psd, i, j);

            /* if parents are not ancestors - should not be the case, if they are added in the beginning */
	    if ( arrow_i_j != NULL && !SCIPisFeasGE(scip, valIJ, SCIPgetSolVal(scip, NULL, arrow_i_j)) )
            {
	       SCIP_ROW *row;
               SCIP_Bool infeasible;

	       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "paisanc#%d#%d", i, j);

	       SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
	       SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], -1.0) );
               SCIP_CALL( SCIPaddVarToRow(scip, row, arrow_i_j, 1.0) );
	       SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	       SCIPdebug( SCIProwPrint(row, NULL) );
#endif
	       SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
	       SCIP_CALL( SCIPreleaseRow(scip, &row));
	       ++nGen;

               if ( infeasible )
               {
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
            }

            if( minimal && SCIPisFeasGT(scip, valIJ, 0.0) )
            {
               /* must have e.g. (1- I(1<-{2,3})) + I(2<--4) + I(3<--4) + (1-(I(1<--4)) >= 1 */
               for( k = 0; k < psd->nParentSets[i]; k++ )
               {
                  SCIP_Real acc;
                  int l;
                  int pa = -1;
                  
                  acc = (1.0 - valIJ) + (1 - SCIPgetSolVal(scip, NULL, psd->PaVars[i][k]));

                  if( SCIPisFeasGE(scip, acc, 1.0) )
                     continue;

                  for( l = 0; l < psd->nParents[i][k] && SCIPisFeasLT(scip, acc, 1.0); l++ )
                  {
                     pa = psd->ParentSets[i][k][l];
                     if( pa == j )
                        break;
                     acc += SCIPgetSolVal(scip, NULL, vars[psd->ParentSets[i][k][l]][j]);
                  }

                  if( pa == j)
                     continue;
                  
                  if( SCIPisFeasLT(scip, acc, 1.0) )
                  {
                     SCIP_ROW *row;
                     SCIP_Bool infeasible;

                     (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "minim#%d#%d#%d", i, j, k);

                     SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -1.0, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
                     SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
                     SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[i][k], -1.0) );
                     SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], -1.0) );
                     for( l = 0; l < psd->nParents[i][k]; l++ )
                        SCIP_CALL( SCIPaddVarToRow(scip, row, vars[psd->ParentSets[i][k][l]][j], 1.0) );
                     SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
                     SCIPdebug( SCIProwPrint(row, NULL) );
#endif
                     SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
                     SCIP_CALL( SCIPreleaseRow(scip, &row));
                     ++nGen;

                     if ( infeasible )
                     {
                        *result = SCIP_CUTOFF;
                        return SCIP_OKAY;
                     }
                  }
               }
            }
            
	    /* if symmetry equations are violated - should not be the case, if they are added in the beginning */
	    if ( ! SCIPisFeasGE(scip, 1.0 - valIJ, SCIPgetSolVal(scip, NULL, vars[j][i])) )
	    {
	       SCIP_ROW *row;
               SCIP_Bool infeasible;

	       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sym#%d#%d", i, j);

	       SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
	       SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
	       SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][i], 1.0) );
	       SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
	       SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
	       SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
	       SCIP_CALL( SCIPreleaseRow(scip, &row));
	       ++nGen;

               if ( infeasible )
               {
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
	    }

	    /* enforce triangle inequalities */
	    for (k = 0; k < n; ++k)
	    {
	       SCIP_Real sum = 0.0;
	       if (k == i || k == j)
		  continue;

	       sum = valIJ + SCIPgetSolVal(scip, NULL, vars[j][k]) - SCIPgetSolVal(scip, NULL, vars[i][k]);

	       /* if sum > 1.0, i.e., the cut is violated */
	       if ( SCIPisFeasGT(scip, sum, 1.0) ) /* this is the only difference to the separation call */
	       {
		  SCIP_ROW *row;
                  SCIP_Bool infeasible;

		  (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "triangle#%d#%d#%d", i, j, k);

		  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 1.0, FALSE, FALSE, TRUE) );
		  SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
		  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][j], 1.0) );
		  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j][k], 1.0) );
		  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[i][k], -1.0) );
		  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
		  SCIPdebug( SCIPprintRow(scip, row, NULL) );
#endif
		  SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
		  SCIP_CALL( SCIPreleaseRow(scip, &row));
		  ++nGen;

                  if ( infeasible )
                  {
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
	       }
	    }
	 }
      }
      if (nGen > 0)
      {
	 *result = SCIP_SEPARATED;
	 return SCIP_OKAY;
      }
   }
   SCIPdebugMsg(scip, "all partial ordering constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsPartialOrdering)
{  /*lint --e{715}*/
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
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;
      ParentSetData* psd;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "enforcing pseudo solution for partial ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      vars = consdata->vars;
      n = consdata->n;
      psd = consdata->psd;
      
      /* check triangle inequalities */
      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
	    SCIP_Bool oneIJ;
            SCIP_Bool onearrowIJ;
            SCIP_VAR* arrow_i_j;

	    if (j == i)
	       continue;

	    /* the priorities should ensure that the solution is integral */
	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[i][j])) );
	    oneIJ = SCIPisGT(scip, SCIPgetSolVal(scip, NULL, vars[i][j]), 0.5);
            
            arrow_i_j = get_arrow(psd, i, j);
            
            /* the priorities should ensure that the solution is integral */
            assert( arrow_i_j == NULL || SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, arrow_i_j)) );

            onearrowIJ = ( arrow_i_j != NULL && SCIPisGT(scip, SCIPgetSolVal(scip, NULL, arrow_i_j), 0.5) );
            
            /* check that parents are partially ordered */
            if( onearrowIJ && !oneIJ )
            {
               SCIPdebugMsg(scip, "constraint <%s> infeasible (parent not partially ordered).\n", SCIPconsGetName(cons));
               *result = SCIP_INFEASIBLE;
	       return SCIP_OKAY;
	    }

            /* if need minimality and j ancestor of i then either j is a parent of i or exists k such that k parent of i, j ancestor of k */ 
            if( consdata->minimal && oneIJ && !onearrowIJ )
            {
               SCIP_Bool kfound = FALSE;
               SCIP_VAR* arrow_i_k;
               
               for (k = 0; k < n; ++k)
               {

                  if (k == i || k == j)
                     continue;

                  /* the priorities should ensure that the solution is integral */
                  assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[k][j])) );
                  if( !SCIPisGT(scip, SCIPgetSolVal(scip, NULL, vars[k][j]), 0.5) )
                     continue;

                  arrow_i_k = get_arrow(psd, i, k);
                  
                  /* the priorities should ensure that the solution is integral */
                  assert( arrow_i_k == NULL || SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, arrow_i_k)) );
                  if( arrow_i_k != NULL && SCIPisGT(scip, SCIPgetSolVal(scip, NULL, arrow_i_k), 0.5) )
                  {
                     kfound = TRUE;
                     break;
                  }
               }
               if( !kfound )
               {
                  SCIPdebugMsg(scip, "constraint <%s> infeasible (partial order not minimal).\n", SCIPconsGetName(cons));
                  *result = SCIP_INFEASIBLE;
                  return SCIP_OKAY;
               }
            }

	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[j][i])) );

	    if ( oneIJ == SCIPisGT(scip, SCIPgetSolVal(scip, NULL, vars[j][i]), 0.5) )
	    {
	       SCIPdebugMsg(scip, "constraint <%s> infeasible (violated equation).\n", SCIPconsGetName(cons));
	       *result = SCIP_INFEASIBLE;
	       return SCIP_OKAY;
	    }

	    for (k = 0; k < n; ++k)
	    {
	       SCIP_Bool oneJK, oneIK;
	       if (k == i || k == j)
		  continue;

	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[j][k])) );
	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, vars[k][i])) );
	       oneJK = SCIPisGT(scip, SCIPgetSolVal(scip, NULL, vars[j][k]), 0.5);
	       oneIK = SCIPisGT(scip, SCIPgetSolVal(scip, NULL, vars[i][k]), 0.5);

	       /* if triangle inequality is violated */
	       if ( oneIJ && oneJK && !oneIK )
	       {
		  SCIPdebugMsg(scip, "constraint <%s> infeasible (violated triangle ineq.).\n", SCIPconsGetName(cons));
		  *result = SCIP_INFEASIBLE;
		  return SCIP_OKAY;
	       }
	    }
	 }
      }
   }
   SCIPdebugMsg(scip, "all partial ordering constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckPartialOrdering)
{  /*lint --e{715}*/
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
      SCIP_VAR*** vars;
      int i;
      int j;
      int k;
      int n;
      ParentSetData* psd;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "checking partial ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );
      assert( consdata->psd != NULL );
      vars = consdata->vars;
      n = consdata->n;
      psd = consdata->psd;

      /* check triangle inequalities and symmetry equations */
      /* and that parents are partially ordered */
      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
	    SCIP_Bool oneIJ;
            SCIP_Bool onearrowIJ;
            SCIP_VAR* arrow_i_j;
	    if (j == i)
	       continue;
            
	    /* the priorities should ensure that the solution is integral */
            assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[i][j])) );
	    oneIJ = SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[i][j]), 0.5);

            arrow_i_j = get_arrow(psd, i, j);
            
            /* the priorities should ensure that the solution is integral */
            assert( arrow_i_j == NULL || SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, arrow_i_j)) );

            onearrowIJ = ( arrow_i_j != NULL && SCIPisGT(scip, SCIPgetSolVal(scip, sol, arrow_i_j), 0.5) );

            /* check that parents are partially ordered */
            if( onearrowIJ && !oneIJ )
            {
               SCIPdebugMsg(scip, "constraint <%s> infeasible (parent not partially ordered).\n", SCIPconsGetName(cons));
               *result = SCIP_INFEASIBLE;
               if( printreason )
               {
                  SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                  SCIPinfoMessage(scip, NULL, "violation: parent not partially ordered <%s> = %.15g and <%s> = %.15g\n",
                     SCIPvarGetName(vars[i][j]), SCIPgetSolVal(scip, sol, arrow_i_j), 0.5,
                     SCIPvarGetName(vars[j][i]), SCIPgetSolVal(scip, sol, vars[i][j]), 0.5);
               }
               return SCIP_OKAY;
            }
            
            /* if need minimality and j ancestor of i then either j is a parent of i or exists k such that k parent of i, j ancestor of k */ 
            if( consdata->minimal && oneIJ && !onearrowIJ )
            {
               SCIP_Bool kfound = FALSE;
               SCIP_VAR* arrow_i_k;
               
               for (k = 0; k < n; ++k)
               {

                  if (k == i || k == j)
                     continue;

                  /* the priorities should ensure that the solution is integral */
                  assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[k][j])) );
                  if( !SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[k][j]), 0.5) )
                     continue;

                  arrow_i_k = get_arrow(psd, i, k);
                  
                  /* the priorities should ensure that the solution is integral */
                  assert( arrow_i_k == NULL || SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, arrow_i_k)) );
                  if( arrow_i_k != NULL && SCIPisGT(scip, SCIPgetSolVal(scip, sol, arrow_i_k), 0.5) )
                  {
                     kfound = TRUE;
                     break;
                  }
               }
               if( !kfound )
               {
                  SCIPdebugMsg(scip, "constraint <%s> infeasible (partial order not minimal).\n", SCIPconsGetName(cons));
                  *result = SCIP_INFEASIBLE;
                  if( printreason )
                  {
                     SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                     SCIPinfoMessage(scip, NULL, "violation: minimality violated <%s> = %.15g\n",
                        SCIPvarGetName(vars[i][j]), SCIPgetSolVal(scip, sol, vars[i][j]), 0.5);
                  }
                  return SCIP_OKAY;
               }
            }
            
	    /* the priorities should ensure that the solution is integral */
	    assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[j][i])) );

	    /* check symmetry equations */
	    if ( oneIJ && SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[j][i]), 0.5) )
	    {
	       SCIPdebugMsg(scip, "constraint <%s> infeasible (violated equation).\n", SCIPconsGetName(cons));
	       *result = SCIP_INFEASIBLE;
               if( printreason )
               {
                  SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                  SCIPinfoMessage(scip, NULL, "violation: symmetry equation violated <%s> = %.15g and <%s> = %.15g\n",
                     SCIPvarGetName(vars[i][j]), SCIPgetSolVal(scip, sol, vars[i][j]), 0.5,
                     SCIPvarGetName(vars[j][i]), SCIPgetSolVal(scip, sol, vars[j][i]), 0.5);
               }
	       return SCIP_OKAY;
	    }

	    for (k = 0; k < n; ++k)
	    {
	       SCIP_Bool oneJK, oneIK;
	       if (k == i || k == j)
		  continue;

	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[j][k])) );
	       assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, vars[k][i])) );
	       oneJK = SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[j][k]), 0.5);
	       oneIK = SCIPisGT(scip, SCIPgetSolVal(scip, sol, vars[i][k]), 0.5);

	       /* if triangle inequality is violated */
	       if ( oneIJ && oneJK && !oneIK )
	       {
		  SCIPdebugMsg(scip, "constraint <%s> infeasible (violated triangle ineq.).\n", SCIPconsGetName(cons));
		  *result = SCIP_INFEASIBLE;
                  if( printreason )
                  {
                     SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                     SCIPinfoMessage(scip, NULL,
                        "violation: triangle inequality violated <%s> = %.15g, <%s> = %.15g, <%s> = %.15g\n",
                        SCIPvarGetName(vars[i][j]), SCIPgetSolVal(scip, sol, vars[i][j]), 0.5,
                        SCIPvarGetName(vars[j][k]), SCIPgetSolVal(scip, sol, vars[j][k]), 0.5,
                        SCIPvarGetName(vars[i][k]), SCIPgetSolVal(scip, sol, vars[i][k]), 0.5);
                  }
		  return SCIP_OKAY;
	       }
	    }
	 }
      }
   }
   SCIPdebugMsg(scip, "all partial ordering constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropPartialOrdering)
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
      SCIP_VAR*** vars;
      SCIP_VAR* arrow_var;
      int i;
      int j;
      int k;
      int n;
      ParentSetData*  psd;
      SCIP_Bool minimal;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMsg(scip, "propagating linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->vars != NULL );

      vars = consdata->vars;
      n = consdata->n;
      psd = consdata->psd;
      minimal = consdata->minimal;
      
      /* check triangle inequalities */
      for (i = 0; i < n; ++i)
      {
	 for (j = 0; j < n; ++j)
	 {
            SCIP_Bool oneIJ;
            
	    if (j == i)
	       continue;

            oneIJ = ( SCIPvarGetLbLocal(vars[i][j]) > 0.5 );

	    /* if x[i][j] == 1 then x[j][i] = 0 */
	    if ( oneIJ )
	    {
	       SCIP_Bool infeasible, tightened;
               SCIP_CALL( SCIPtightenVarUb(scip,vars[j][i],0,TRUE,&infeasible,&tightened) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, " -> node infeasible.\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
	       /* SCIP_CALL( SCIPinferBinvarCons(scip, vars[j][i], FALSE, cons, i*n + j, &infeasible, &tightened) ); */
	       /* if ( infeasible ) */
	       /* { */
	       /*    SCIPdebugMessage(" -> node infeasible.\n"); */
               /*    SCIP_CALL( SCIPinitConflictAnalysis(scip) ); */
               /*    SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) ); */
               /*    SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][i]) ); */
               /*    SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) ); */
	       /*    *result = SCIP_CUTOFF; */
	       /*    return SCIP_OKAY; */
	       /* } */
	       if ( tightened )
		  ++nGen;
	    }

	    /* if i<-j == 1 then i<--j = 1 */
            arrow_var = get_arrow(psd,i,j);
	    if ( arrow_var != NULL && SCIPvarGetLbLocal(arrow_var) > 0.5 )
            {
               SCIP_Bool infeasible, tightened;
               SCIP_CALL( SCIPtightenVarLb(scip,vars[i][j],1,TRUE,&infeasible,&tightened) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, " -> node infeasible.\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  SCIPdebugMsg(scip, "Setting %s to 1\n", SCIPvarGetName(vars[i][j]));
                  ++nGen;
               }
            }


            /* if i<--j == 0 then i<-j = 0 */
	    if ( arrow_var != NULL && SCIPvarGetUbLocal(vars[i][j]) < 0.5 )
            {
               SCIP_Bool infeasible, tightened;
               SCIP_CALL( SCIPtightenVarUb(scip,arrow_var,0,TRUE,&infeasible,&tightened) );
               if( infeasible )
               {
                  SCIPdebugMsg(scip, " -> node infeasible.\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               if( tightened )
               {
                  SCIPdebugMsg(scip, "Setting %s to 0\n", SCIPvarGetName(arrow_var));
                  ++nGen;
               }
            }

            /* if minimal ...
               if parent set for i has been chosen and j is not an ancestor of any of these parents
               then not an ancestor of i either 
            */
            if( minimal )
            {
               for( k = 0; k < psd->nParentSets[i]; k++ )
               {
                  if ( SCIPvarGetLbLocal(psd->PaVars[i][k]) > 0.5 )
                  {
                     SCIP_Bool allnonancestors = TRUE;
                     int l;
                     int pa;
                     
                     for( l = 0; l < psd->nParents[i][k]; l++ )
                     {
                        pa = psd->ParentSets[i][k][l];
                        if( pa == j || SCIPvarGetUbLocal(vars[psd->ParentSets[i][k][l]][j]) > 0.5 )
                        {
                           allnonancestors = FALSE;
                           break;
                        }
                     }
                     
                     if( allnonancestors )
                     {
                        SCIP_Bool infeasible, tightened;
                        SCIP_CALL( SCIPtightenVarUb(scip,vars[i][j],0,TRUE,&infeasible,&tightened) );
                        if( infeasible )
                        {
                           SCIPdebugMsg(scip, " -> node infeasible.\n");
                           *result = SCIP_CUTOFF;
                           return SCIP_OKAY;
                        }
                        if( tightened )
                        {
                           SCIPdebugMsg(scip, "Setting %s to 0\n", SCIPvarGetName(vars[i][j]));
                           ++nGen;
                        }
                     }
                  }
                  else if( SCIPvarGetUbLocal(psd->PaVars[i][k]) < 0.5 )
                     /* if this one not set (either way) then no parent
                        set is set to 1 */
                     break;
               }
            }

            if( !oneIJ )
               continue;
            
	    for (k = 0; k < n; ++k)
	    {
	       if (k == i || k == j)
		  continue;

	       /* if x[i][j] == 1 and x[j][k] == 1 then x[i][k] = 1 */
               /* to get to here have to have oneIJ== TRUE, so don't test again */
	       if ( SCIPvarGetLbLocal(vars[j][k]) > 0.5 )
	       {
		  SCIP_Bool infeasible, tightened;
                  SCIP_CALL( SCIPtightenVarLb(scip,vars[i][k],1,TRUE,&infeasible,&tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMsg(scip, " -> node infeasible.\n");
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }

		  /* SCIP_CALL( SCIPinferBinvarCons(scip, vars[i][k], TRUE, cons, n*n + i*n*n + j*n + k, &infeasible, &tightened) ); */
		  /* if ( infeasible ) */
		  /* { */
		  /*    SCIPdebugMessage(" -> node infeasible.\n"); */
                  /*    SCIP_CALL( SCIPinitConflictAnalysis(scip) ); */
                  /*    SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][j]) ); */
                  /*    SCIP_CALL( SCIPaddConflictBinvar(scip, vars[j][k]) ); */
                  /*    SCIP_CALL( SCIPaddConflictBinvar(scip, vars[i][k]) ); */
                  /*    SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) ); */
		  /*    *result = SCIP_CUTOFF; */
		  /*    return SCIP_OKAY; */
		  /* } */
		  if ( tightened )
		     ++nGen;
	       }

	       /* all other implications occur with other indices i, j, k */
	    }
	 }
      }
   }
   if (nGen > 0)
      *result = SCIP_REDUCEDDOM;
   SCIPdebugMsg(scip, "propagated %d domains.\n", nGen);

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropPartialOrdering)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( bdchgidx != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Propagation resolution of constraint <%s>.\n", SCIPconsGetName(cons));
   *result = SCIP_DIDNOTFIND;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->vars != NULL );

   n = consdata->n;
   vars = consdata->vars;

   assert( 0 <= inferinfo && inferinfo < n*n + n*n*n );

   /* if the conflict came from an equation */
   if ( inferinfo < (n*n) )
   {
      int index1;
      int index2;

      index1 = inferinfo/n;
      index2 = inferinfo % n;
      assert( 0 <= index1 && index1 < n );
      assert( 0 <= index2 && index2 < n );
      assert( vars[index2][index1] == infervar );

      /* if the variable was fixed to 1 */
      if ( SCIPvarGetLbAtIndex(infervar, bdchgidx, FALSE) < 0.5 && SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > 0.5 )
      {
	 SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 1 was x[%d][%d] = 0.\n", index2, index1, index1, index2);
	 /* the reason was that x[i][j] was fixed to 0 */
	 SCIP_CALL( SCIPaddConflictUb(scip, vars[index1][index2], bdchgidx) );
	 *result = SCIP_SUCCESS;
	 return SCIP_OKAY;
      }
   }
   else
   {
      /* otherwise the conflict came from a triangle inequality */
      int index1;
      int index2;
      int index3;

      index1 = (inferinfo - n*n)/(n*n);
      index2 = (inferinfo - n*n - index1 * n*n)/n;
      index3 = (inferinfo - n*n) % n;

      assert( 0 <= index1 && index1 < n );
      assert( 0 <= index2 && index2 < n );
      assert( 0 <= index3 && index3 < n );
      assert( index1 != index2 && index2 != index3 && index1 != index3 );
      assert( vars[index1][index3] == infervar );

      /* the variable should have been fixed to 1 */
      assert( SCIPvarGetUbAtIndex(infervar, bdchgidx, FALSE) < 0.5 && SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) > 0.5 );

      /* the reason was that x[index1][index2] and x[index2][index3] were fixed to 1 */
      SCIPdebugMsg(scip, " -> reason for x[%d][%d] == 1 was x[%d][%d] = x[%d][%d] = 1.\n", index1, index3, index1, index2, index2, index3);
      SCIP_CALL( SCIPaddConflictLb(scip, vars[index1][index2], bdchgidx) );
      SCIP_CALL( SCIPaddConflictLb(scip, vars[index2][index3], bdchgidx) );
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}
#else
#define consRespropPartialOrdering NULL
#endif

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockPartialOrdering)
{  /*lint --e{715}*/
   int i;
   int j;
   int k;
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int n;
   ParentSetData* psd;
   SCIP_VAR* arrow_var;
   
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMsg(scip, "Locking linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);
   assert( consdata->vars != NULL );
   assert( consdata->psd != NULL );
   n = consdata->n;
   vars = consdata->vars;
   psd = consdata->psd;
   
   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
      {
	 if (i != j)
	 {
	    /* the constraint may be violated in any way */
	    SCIP_CALL( SCIPaddVarLocksType(scip, vars[i][j], SCIP_LOCKTYPE_MODEL, nlockspos + nlocksneg, nlockspos + nlocksneg) );

            arrow_var = get_arrow(psd,i,j);
            if( arrow_var != NULL )
               /* the constraint may be violated in any way */
               SCIP_CALL( SCIPaddVarLocks(scip, arrow_var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
	 }
      }
      for( k = 0; k < psd->nParentSets[i]; k++ )
         /* the constraint may be violated in any way */
         SCIP_CALL( SCIPaddVarLocks(scip, psd->PaVars[i][k], nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintPartialOrdering)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int i;
   int j;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->vars != NULL );
   n = consdata->n;
   vars = consdata->vars;

   SCIPinfoMessage(scip, file, "partialordering[");
   for (i = 0; i < n; ++i)
   {
      if ( i > 0 )
	 SCIPinfoMessage(scip, file, ", ");
      SCIPinfoMessage(scip, file, "(");
      for (j = 0; j < n; ++j)
      {
	 if (j != i)
	 {
	    if ( j > 0 && (i > 0 || j > 1) )
	       SCIPinfoMessage(scip, file, ",");
	    SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(vars[i][j]));
	 }
      }
      SCIPinfoMessage(scip, file, ")");
   }
   SCIPinfoMessage(scip, file, "]\n");

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyPartialOrdering)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_VAR*** sourcevars;
   SCIP_VAR*** vars;
   int i;
   int j;
   int k;
   int n;

   SCIP_Bool minimal;
   ParentSetData* targetpsd;
   SCIP_VAR* arrow_var;
   SCIP_VAR* new_arrow_var;
   
   assert( scip != 0 );
   assert( sourceconshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != 0 );
   assert( sourcescip != 0 );
   assert( sourcecons != 0 );
   assert( varmap != 0 );

   *valid = TRUE;

   SCIPdebugMsg(scip, "Copying method for linear ordering constraint handler.\n");

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   n = sourcedata->n;
   minimal = sourcedata->minimal;

   /* copy parent set data, but does not create new variables */
   SCIP_CALL( PS_copyParentSetData(sourcescip, sourcedata->psd, &targetpsd) );
   /* have to reinitiaise hashtable for edges and arrows since new variables
      will be inserted using put_arrow (see below) */
   SCIP_CALL( hashtablefreeArrow(sourcescip, targetpsd) );
   SCIP_CALL( hashtableCreateArrow(sourcescip, targetpsd) );
   
   sourcevars = sourcedata->vars;
   assert( sourcevars != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, n) );
   BMSclearMemoryArray(vars, n);

   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(vars[i]), n) ); /*lint !e866*/

      for (j = 0; j < n && *valid; ++j)
      {
         if ( i != j )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcevars[i][j], &vars[i][j], varmap, consmap, global, valid) );
            assert( !(*valid) || vars[i][j] != NULL );
         }

         arrow_var = get_arrow(sourcedata->psd,i,j);
         if( arrow_var != NULL )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, arrow_var, &new_arrow_var, varmap, consmap, global, valid) );
            assert( !(*valid) || new_arrow_var != NULL );
            put_arrow(scip, targetpsd, i, j, new_arrow_var);
         }
      }

      for( k = 0; k < sourcedata->psd->nParentSets[i]; k++ )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->psd->PaVars[i][k], &(targetpsd->PaVars[i][k]),
               varmap, consmap, global, valid) );
         assert( !(*valid) || targetpsd->PaVars[i][k] != NULL );
      }

   }

   if ( *valid )
   {
      /* create copied constraint */
      if ( name == 0 )
         name = SCIPconsGetName(sourcecons);

      SCIP_CALL( SCIPcreateConsPartialOrdering(scip, cons, name, n, vars, targetpsd, minimal,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   for (i = 0; i < n; ++i)
      SCIPfreeBufferArrayNull(scip, &vars[i]);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** creates the handler for partial ordering constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrPartialOrdering(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* include constraint handler */
   conshdlr = NULL;
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpPartialOrdering, consEnfopsPartialOrdering, consCheckPartialOrdering, consLockPartialOrdering,
         NULL) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeletePartialOrdering) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyPartialOrdering, consCopyPartialOrdering) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransPartialOrdering) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpPartialOrdering) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpPartialOrdering, consSepasolPartialOrdering,
         CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropPartialOrdering, CONSHDLR_PROPFREQ,
         CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   /* SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropPartialOrdering) ); */
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintPartialOrdering) );

   return SCIP_OKAY;
}

/** creates and captures a partial ordering constraint */
SCIP_RETCODE SCIPcreateConsPartialOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of binary variables */
   ParentSetData*        psd,                /**< contains arrow and parent set variables */
   SCIP_Bool             minimal,            /**< whether the partial order must minimal while being consistent with arrows */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   assert(vars != NULL);
   assert(psd != NULL);
   
   /* find the partial ordering constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if (conshdlr == NULL)
   {
      SCIPerrorMessage("partial ordering constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->n = n;
   consdata->minimal = minimal;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, n) );
   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->vars[i]), n) ); /*lint !e866*/
      for (j = 0; j < n; ++j)
      {
	 if (j != i)
	 {
	    assert( vars[i][j] != NULL );
	    consdata->vars[i][j] = vars[i][j];
	 }
      }
   }

   SCIP_CALL( PS_copyParentSetData(scip, psd, &(consdata->psd)) );
   assert(consdata->psd != NULL);
   
   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a partial ordering constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicPartialOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of binary variables */
   ParentSetData*        psd,                /**< contains arrow and parent set variables */
   SCIP_Bool             minimal             /**< whether the partial order must minimal while being consistent with arrows */
   )
{
   SCIP_CALL( SCIPcreateConsPartialOrdering(scip, cons, name, n, vars, psd, minimal,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
