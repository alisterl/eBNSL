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

/**@file   cons_vanilla.c
 * @brief  constraint handler for vanilla constraints
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*#define SCIP_DEBUG*/
#include <assert.h>

#include "cons_vanilla.h"
#include "utils.h"
#include <string.h>
#include "parent_set_data.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "vanilla"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
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
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */




/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for vanilla constraints */
struct SCIP_ConsData
{
   ParentSetData* psd;
   SCIP_VAR*** ancestorvars;
   int** n_betterpa;
   int*** betterpa;
   SCIP_Bool*** store;                /**< store[i][k][j] = TRUE if j is in kth parent set for i */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** Returns true if adding parent to any parent set for child increases the score 
 */
static
SCIP_Bool arrow_preferred(
   ParentSetData* psd,   /**< parent sets data structure */
   SCIP_Bool*** store,
   int child,            /**< child BN variable */
   int parent            /**< parent BN variable */
   )
{
   SCIP_Bool ok;
   SCIP_Bool found_better;
   int j;
   int l;
   int k;
   int kk;
   
   
   for( k = 0; k < psd->nParentSets[child]; ++k )
   {
      if( store[child][k][parent] )
         continue;

      found_better = FALSE;
      for( kk = 0; kk < k; ++kk )
      {
         if( !store[child][kk][parent] )
            continue;
         
         ok = TRUE;
         for( l = 0; l < psd->nParents[child][kk]; ++l )
         {
            j = psd->ParentSets[child][kk][l];
            if( j != parent && !store[child][k][j] )
            {
               ok = FALSE;
               break;
            }
         }

         if( ok )
         {
            found_better = TRUE;
            break;
         }
      }
      if( !found_better )
      {
         return FALSE;
      }
   }

   return TRUE;
}



/** Creates the data for a constraint.
 *  @param scip The SCIP instance to which the constraaint belongs.
 *  @param consdata The location to store the new constraint data.
 *  @param psd The parent set data on which the constraint is based.
 *  @return SCIP_OKAY if successful, or an appropriate error otherwise.
 */
static
SCIP_RETCODE createConsData(
   SCIP* scip,
   SCIP_CONSDATA** consdata,
   ParentSetData* psd,
   SCIP_VAR*** ancestorvars
   )
{

   int i;
   int j;
   int k;

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   SCIP_CALL( PS_copyParentSetData(scip, psd, &((*consdata)->psd)) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->n_betterpa), psd->n) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->betterpa), psd->n) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store), psd->n) );
   
   if( ancestorvars != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->ancestorvars), psd->n) );
      for (i = 0; i < psd->n; ++i)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->ancestorvars[i]), psd->n) ); /*lint !e866*/
         for (j = 0; j < psd->n; ++j)
         {
            if (j != i)
            {
               assert( ancestorvars[i][j] != NULL );
               (*consdata)->ancestorvars[i][j] = ancestorvars[i][j];
            }
         }
      }
   }

   for (i = 0; i < psd->n; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->n_betterpa[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->betterpa[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store[i]), (*consdata)->psd->nParentSets[i]) );

      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         int kk;
         int l;
         int ll;
         int extra;
         int nkpa;

         SCIP_CALL( SCIPallocMemoryArray(scip, &((*consdata)->store[i][k]), (*consdata)->psd->n) );
         for( l = 0; l < (*consdata)->psd->n; ++l )
            (*consdata)->store[i][k][l] = FALSE;
         for( l = 0; l < (*consdata)->psd->nParents[i][k]; ++l )
            (*consdata)->store[i][k][(*consdata)->psd->ParentSets[i][k][l]] = TRUE;

         (*consdata)->n_betterpa[i][k] = 0;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*consdata)->betterpa[i][k]), psd->n) );

         /* look for a better parent set with one extra parent */
         nkpa = psd->nParents[i][k];
         for( kk = 0; kk < k; ++kk )
         {
            if( psd->nParents[i][kk] != nkpa+1 )
               continue;

            extra = -1;
            l = 0;
            for( ll = 0; ll < psd->nParents[i][kk]; ++ll)
            {
               if( psd->ParentSets[i][kk][ll] == psd->ParentSets[i][k][l] )
                  l++;
               else if( extra == -1 )
                  extra = psd->ParentSets[i][kk][ll];
               else
               {
                  extra = -1;
                  break;
               }
            }

            if( extra != -1 )
               (*consdata)->betterpa[i][k][(*consdata)->n_betterpa[i][k]++] = extra;
         }
         /* printf("%d %d: %d better parents.\n",i,k,(*consdata)->n_betterpa[i][k]); */
      }
   }
   return SCIP_OKAY;
}



/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteVanilla)
{  /*lint --e{715}*/

   int i;
   int k;
   int n;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( consdata != NULL);
   assert( *consdata != NULL);
   assert( (*consdata)->psd != NULL );

   SCIPdebugMessage("deleting vanilla constraint <%s>.\n", SCIPconsGetName(cons));

   n = (*consdata)->psd->n;

   for( i = 0; i < n; ++i )
   {
      for( k = 0;  k < (*consdata)->psd->nParentSets[i]; ++k )
         SCIPfreeMemoryArray(scip, &((*consdata)->store[i][k]));
      SCIPfreeMemoryArray(scip, &((*consdata)->store[i]));
   }
   SCIPfreeMemoryArray(scip, &((*consdata)->store));

   SCIP_CALL( PS_deallocateParentSetData(scip, &((*consdata)->psd), FALSE) );

   if( (*consdata)->ancestorvars != NULL )
   {
      for (i = 0; i < n; ++i)
         SCIPfreeBlockMemoryArray(scip, &((*consdata)->ancestorvars[i]), n); /*lint !e866*/
      SCIPfreeBlockMemoryArray(scip, &((*consdata)->ancestorvars), n);
   }


   SCIPfreeBlockMemory(scip, consdata);
   
   return SCIP_OKAY;
}

/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpVanilla)
{  /*lint --e{715}*/

   int c;
   int nGen = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {

      SCIP_CONSDATA* consdata;
      ParentSetData* psd; 
      int i;
      int j;
      int j2;
      int k;
      
      SCIP_CONS* sc_cons;
      SCIP_VAR** vars;
      int nvars;
      
      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMessage("adding initial rows for vanilla constraint <%s>.\n", SCIPconsGetName(conss[c]));
      
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      assert(consdata->psd != NULL);
      
      psd = consdata->psd;
      
      SCIP_CALL( SCIPallocMemoryArray(scip, &vars, psd->n) );
   
      for( i = 0; i < psd->n; ++i )
         vars[i] = psd->PaVars[i][0];

      /* add a constraint rather than a row, but how is this deleted? */
      SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &sc_cons, "sink", psd->n, vars) );
      SCIP_CALL( SCIPaddCons(scip, sc_cons) );
      nGen++;
      /* SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); */
      SCIP_CALL( SCIPreleaseCons(scip, &sc_cons) );

      SCIPfreeMemoryArray(scip, &vars);

      /* if no ancestorvars then no further initial rows for this constraint */
      if( consdata->ancestorvars == NULL )
         continue;
      
      for( i = 0; i < psd->n; ++i )
      {
         for( j2 = 0; j2 < psd->n; ++j2 )
         {
            if( j2 == i )
               continue;

            nvars = 0;
            SCIP_CALL( SCIPallocMemoryArray(scip, &vars, psd->nParentSets[i]) );
            for( k = 0; k < psd->nParentSets[i]; ++k )
            {
               for( j = 0; j < consdata->n_betterpa[i][k]; ++j )
                  if( consdata->betterpa[i][k][j] == j2 )
                  {
                     vars[nvars++] = psd->PaVars[i][k];
                     break;
                  }
            }

            {
               char s[SCIP_MAXSTRLEN];
               SCIP_ROW* row;

               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "critic_anc#%d#%d#%d", i, j2);
               SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s, -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
               for( j = 0; j < nvars; ++j )
                  SCIP_CALL( SCIPaddVarToRow(scip, row, vars[j], 1.0) );
               SCIPfreeMemoryArray(scip, &vars);
               SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->ancestorvars[j2][i], -1.0) );
               SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
               SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ));
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
   }
   
   return SCIP_OKAY;
}





/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpVanilla)
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
      *result = SCIP_FEASIBLE;
   }
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsVanilla)
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
      *result = SCIP_FEASIBLE;
   }
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckVanilla)
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
      int i;
      int k;
      int l;
      
      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("checking vanilla constraint <%s>.\n", SCIPconsGetName(cons));


      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );

      /* be very lazy if no ancestor vars! */
      if( consdata->ancestorvars == NULL )
         continue;

      for( i = 0; i < consdata->psd->n; ++i )
      {
         for( k = 0; k < consdata->psd->nParentSets[i]; ++k )
         {
            SCIP_Bool allok = TRUE;
            assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, consdata->psd->PaVars[i][k])) );
            
            for( l = 0; l < consdata->psd->nParents[i][k]; ++l )
            {
               assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, consdata->ancestorvars[consdata->psd->ParentSets[i][k][l]][i])) );
               if( SCIPisGT(scip, SCIPgetSolVal(scip, sol, consdata->ancestorvars[consdata->psd->ParentSets[i][k][l]][i]), 0.5) )
               {
                  allok = FALSE;
                  break;
               }
            }

            if( allok )
            {
               if( SCIPisLT(scip, SCIPgetSolVal(scip, sol, consdata->psd->PaVars[i][k]), 0.5) )
               {
                  SCIPdebugMessage("constraint <%s> infeasible (Could have chosen <%s>).\n",
                     SCIPconsGetName(cons), SCIPvarGetName(consdata->psd->PaVars[i][k]) );
                  *result = SCIP_INFEASIBLE;
                  return SCIP_OKAY;
               }
               else
               {
                  break;
               }
            }
         }
      }
   }

   SCIPdebugMessage("all vanilla constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropVanilla)
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
      SCIP_VAR*** ancestorvars;
      ParentSetData* psd;
      int i;
      int j;
      int k;
      int kk;
      int n;
      SCIP_Bool infeasible, tightened;

      SCIP_Bool inotancanypa;
      int l;
      
      int** n_betterpa;
      int*** betterpa;
      
      cons = conss[c];
      assert( cons != NULL );
      /* SCIPdebugMessage("propagating vanilla constraint <%s>.\n", SCIPconsGetName(cons)); */

      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );

      /* only do propagations if ancestor variables are available */

      if( consdata->ancestorvars == NULL )
         continue;
      
      ancestorvars = consdata->ancestorvars;
      n = consdata->psd->n;
      psd = consdata->psd;
      n_betterpa = consdata->n_betterpa;
      betterpa = consdata->betterpa;
      
      for( i = 0; i < n; ++i )
      {
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            /* if this one selected then certain ancestor relations must obtain,
               (if they didn't some better parent set would be chosen)
            */
            if( SCIPvarGetLbLocal(psd->PaVars[i][k]) > 0.5 )
            {
               for( j = 0; j < n_betterpa[i][k]; ++j )
               {
                  SCIP_CALL( SCIPtightenVarLb(scip,ancestorvars[betterpa[i][k][j]][i],1,TRUE,&infeasible,&tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible.\n");
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                  {
                     SCIPdebugMessage("Setting ancestor relation (prop #1) %s\n", SCIPvarGetName(ancestorvars[betterpa[i][k][j]][i]));
                     ++nGen;
                  }
               }
               continue;
            }

            /* now the other way */
            for( j = 0; j < n_betterpa[i][k]; ++j )
            {
               /* if i cannot be an ancestor of betterpa[i][k][j] then some better scoring
                  parent set can always be chosen in preference to the kth
               */
               if( SCIPvarGetUbLocal(ancestorvars[betterpa[i][k][j]][i]) < 0.5 )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip,psd->PaVars[i][k],0,TRUE,&infeasible,&tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible.\n");
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                  {
                     SCIPdebugMessage("Ruling out parent set (prop #1) %s\n", SCIPvarGetName(psd->PaVars[i][k]));
                     ++nGen;
                  }
                  break;
               }
            }

            inotancanypa = TRUE;
            for (l = 0; l < psd->nParents[i][k]; ++l)
            {
               /* if not anc[pa][i] == 0 then no prop */
               if( SCIPvarGetUbLocal(ancestorvars[psd->ParentSets[i][k][l]][i]) > 0.5 )
               {
                  inotancanypa = FALSE;
                  break;
               }
            }

            /* kth parent set can always be chosen, rule out all that have worse score */
            
            if( inotancanypa )
            {
               SCIPdebugMessage("Choosing parent set %s will not lead to a cycle\n", SCIPvarGetName(psd->PaVars[i][k]));
               
               for( kk = k+1; kk < psd->nParentSets[i]; ++kk )
               {     
                  SCIP_CALL( SCIPtightenVarUb(scip,psd->PaVars[i][kk],0,TRUE,&infeasible,&tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMessage(" -> node infeasible.\n");
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                  {
                     SCIPdebugMessage("Ruling out parent set (prop #1) %s\n", SCIPvarGetName(psd->PaVars[i][kk]));
                     ++nGen;
                  }
               }

               /* no point looking for further propagations */
               break;
            }
         }
      }
            
      /* /\* for each i, for the best parent set, if i is not an ancestor for any of them, */
      /*    then select it */
      /* *\/ */

      /* for (i = 0; i < n; ++i) */
      /* { */
      /*    SCIP_Bool inotancanypa = TRUE; */
      /*    int l; */
         
      /*    for (l = 0; l < psd->nParents[i][0]; ++l) */
      /*    { */
      /*       /\* if not anc[pa][i] == 0 then no prop *\/ */
      /*       if( SCIPvarGetUbLocal(ancestorvars[psd->ParentSets[i][0][l]][i]) > 0.5 ) */
      /*       { */
      /*          inotancanypa = FALSE; */
      /*          break; */
      /*       } */
      /*    } */

      /*    if( inotancanypa ) */
      /*    { */
      /*       SCIP_Bool infeasible, tightened; */
      /*       SCIP_CALL( SCIPtightenVarLb(scip,psd->PaVars[i][0],1,TRUE,&infeasible,&tightened) ); */
      /*       if( infeasible ) */
      /*       { */
      /*          SCIPdebugMessage(" -> node infeasible.\n"); */
      /*          *result = SCIP_CUTOFF; */
      /*          return SCIP_OKAY; */
      /*       } */
      /*       if ( tightened ) */
      /*       { */
      /*          SCIPdebugMessage("Selecting best parent set %s\n", SCIPvarGetName(psd->PaVars[i][0])); */
      /*          ++nGen; */
      /*       } */
      /*    } */
      /* } */
   }

   return SCIP_OKAY;
}



/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolVanilla)
{

   int c;
   
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   for( c = 0; c < nconss && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      ParentSetData* psd;
      
      int i;
      int j;
      int l;

      SCIP_Bool infeasible, fixed;
      SCIP_VAR* edge_i_j;

      assert(*result != SCIP_CUTOFF);
      
      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);
      assert(consdata->psd != NULL);
      
      SCIPdebugMessage("presolving vanilla constraint <%s>\n", SCIPconsGetName(cons));
      
      psd = consdata->psd;
      
      for( i = 0; i < psd->n; ++i )
      {
         for( l = 0; l < psd->nParents[i][0]; ++l )
         {
            j = psd->ParentSets[i][0][l];

            if( j < i )
               continue;

            /* quick check is: i in the best parent set for j? */
            if( !consdata->store[j][0][i] )
               continue;

            if( arrow_preferred(psd, consdata->store, i, j) && arrow_preferred(psd, consdata->store, j, i) )
            {
               edge_i_j = get_edge(psd,i,j);
               if( edge_i_j == NULL )
               {
                  SCIPdebugMessage("vanilla constraint <%s>: infeasible fixing non-existent edge var == 1\n",
                     SCIPconsGetName(cons) );
                  
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
                  
               /* fix edge variable to 1 */
               SCIP_CALL( SCIPfixVar(scip, edge_i_j, 1.0, &infeasible, &fixed) );
               
               if( infeasible )
               {
                  SCIPdebugMessage("vanilla constraint <%s>: infeasible fixing <%s> == 1\n",
                     SCIPconsGetName(cons), SCIPvarGetName(edge_i_j));

                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }

               if( fixed )
               {
                  SCIPdebugMessage("vanilla constraint <%s>: successfully fixing <%s> == 1\n",
                     SCIPconsGetName(cons), SCIPvarGetName(edge_i_j));

                  *result = SCIP_SUCCESS;
               }
            }
         }
      }
   }
   
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockVanilla)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;
   ParentSetData* psd; 
   int i;
   int j;
   int k;
   SCIP_VAR* var;
      

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->psd != NULL);

   psd = consdata->psd;

   /* lock all variables in both directions */
   for( i = 0; i < psd->n; ++i )
   {
      for( k = 0; k < psd->nParentSets[i]; ++k )
         SCIP_CALL( SCIPaddVarLocks(scip, psd->PaVars[i][k], nlockspos + nlocksneg, nlockspos + nlocksneg) );
      
      for( j = 0; j < psd->n; ++j )
      {
         if( i != j )
         {
            var = get_arrow(psd, i, j);
            if( var != NULL )
               SCIP_CALL( SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
         }
         if( i < j )
         {
            var = get_edge(psd, i, j);
            if( var != NULL )
               SCIP_CALL( SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
         }
      }
   }
   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for vanilla constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrVanilla(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create vanilla constraint handler data */
   conshdlrdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */

   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpVanilla, consEnfopsVanilla, consCheckVanilla, consLockVanilla,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpVanilla) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolVanilla, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteVanilla) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropVanilla, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );

   /* add vanilla constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}

/** creates and captures a vanilla constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsVanilla(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,                /**< parent set data */
   SCIP_VAR***           ancestorvars,       /**< ancestor vars (or NULL if not used ) */
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsVanilla() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the vanilla constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("vanilla constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( createConsData(scip, &consdata, psd, ancestorvars) );

   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a vanilla constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicVanilla(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,                /**< parent set data */
   SCIP_VAR***           ancestorvars        /**< ancestor vars (or NULL if not used ) */
   )
{
   SCIP_CALL( SCIPcreateConsVanilla(scip, cons, name, psd, ancestorvars,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
