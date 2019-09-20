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

/**@file   cons_ci.c
 * @brief  constraint handler for conditional independence constraints
 * @author James Cussens
 *
 *  Doxygen documentation for the locally defined structs @c SCIP_ConsData and @c SCIP_ConshdlrData is unfortunately not available
 *  since structs with the same name are defined in other constraint handlers and Doxygen cannot handle the name clash. You will have to consult
 *  the source code to view this documentation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*#define SCIP_DEBUG*/
#include <assert.h>
#include "cons_ci.h"
#include "parent_set_data.h"
#include "string.h"
#include <ctype.h>
#include "scip/scipdefplugins.h"
#include "solution_info.h"
#include "subip_cuts.h"
#include "utils.h"
#define max(A,B) ((A) > (B) ? (A) : (B))

#define DEFAULT_TIMELIMIT   1e+20
#define DEFAULT_GAPLIMIT  0
#define DEFAULT_ABSGAPLIMIT  0
#define DEFAULT_FORCECUTS TRUE

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "ci"
#define CONSHDLR_DESC          "conditional independence constraint handler"
#define CONSHDLR_SEPAPRIORITY      200000000           /**< priority of the constraint handler for separation (above dagcluster) */
#define CONSHDLR_ENFOPRIORITY         -1 /**< -1, only DAGs are checked. priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        -900000 /**< -1, only DAGs are checked. priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_PRESOLTIMING      SCIP_PRESOLTIMING_ALWAYS

/*
 * Data structures
 */

/*
 * constraint data for ci constraints 
*/
struct SCIP_ConsData
{
   /* ci constraint is that a _|_ b | s
      ( a is independent of b given s ),
      where a, b and s are sets of nodes
   */
   ParentSetData* psd; /* parent set data containing the variables in the constraint */
   SCIP_VAR*** ancestorvars; /* ancestorvars (or NULL if absent) */
   int* a;             /* the set a */
   int* b;             /* the set b */
   int* s;             /* the set s */
   int n_a;            /* number of variables in a */
   int n_b;            /* number of variables in b */
   int n_s;            /* number of variables in s */
   SCIP_Bool* in_abs;  /* in_abs[i]=TRUE iff i is in a \cup b \cup s */
   SCIP_Bool* in_a;    /* in_a[i]=TRUE iff i is in a */
   SCIP_Bool* in_b;    /* in_b[i]=TRUE iff i is in b */
   SCIP_Bool* in_s;    /* in_s[i]=TRUE iff i is in s */
   int*         sol_paset;            /* sol_paset[i] = k if parent set k is chosen in solution under consideration */
   SCIP_Bool*   is_in_an;             /* is_in_an[i]=TRUE iff i is in the smallest ancestral set
                (for the current sol) including a \cup b \cup s */
   int*         an;                   /* integer array representation of the smallest ancestral set
                (for the current sol) including a \cup b \cup s */
   int*         new_an;               /* used to help construct smallest ancestral set
                (for the current sol) including a \cup b \cup s */
   SCIP_Bool**  moral;                /* moral[i][j]=TRUE iff i and j are adjacent in the moral graph of the smallest ancestral set
                (for the current sol) including a \cup b \cup s */
   int**        moral_nbrs;           /* moral_nbrs[i] is the array of neighbours of i in the moral graph of the smallest ancestral set
                (for the current sol) including a \cup b \cup s */
   int*         n_moral_nbrs;         /* n_moral_nbrs[i] is length of moral_nbrs[i] */
   SCIP_Bool*   is_a_reachable;       /* is_a_reachable[i]=TRUE iff i is known to be reachable from some node in a (avoiding s, in moral graph) */
   SCIP_Bool*   is_b_reachable;       /* is_b_reachable[i]=TRUE iff i is known to be reachable from some node in b (avoiding s, in moral graph) */
   int*         a_open;               /* open list for nodes reachable from a (avoiding s, in moral graph) */
   int*         b_open;               /* open list for nodes reachable from b (avoiding s, in moral graph) */
   int*         newnodes;             /* new nodes to add to the frontier */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool forcecuts;               /**< whether to force cuts to be added */
};



/*
 * Local methods
 */

static
SCIP_RETCODE addcicut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< violated ci constraint  */
   SCIP_SOL*             sol                 /**< BN to cut away */
)
{
   SCIP_CONSDATA* consdata;
   int i;
   int k;
   SCIP_ROW* cut;

   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool infeasible;

   ParentSetData* psd;

   /* find the ci constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("ci constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   psd = consdata->psd;

   SCIP_CALL(SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "ci_cut", -SCIPinfinity(scip), (psd->n) - 1,
                                    FALSE, FALSE, TRUE));

   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k])));
         if( SCIPisGT(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k]), 0.5) )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, cut, psd->PaVars[i][k], 1.0) );
            break;
         }
      }

   SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));
   SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );
   return SCIP_OKAY;
}

/** adds an edge to the moral graph of the smallest ancestral set */ 
static
SCIP_RETCODE add_to_moral_graph(
   SCIP_CONSDATA* consdata,        /**< constraint data */ 
   int i,                          /**< first vertex of the edge to add */
   int j                           /**< second vertex of the edge to add */
   )
{
   if( !consdata->moral[i][j] )
   {
      consdata->moral[i][j] = TRUE;
      consdata->moral_nbrs[i][consdata->n_moral_nbrs[i]++] = j;
   }
   if( !consdata->moral[j][i] )
   {
      consdata->moral[j][i] = TRUE;
      consdata->moral_nbrs[j][consdata->n_moral_nbrs[j]++] = i;
   }

   return SCIP_OKAY;
}



/** checks conditional independence constraint for feasibility of given solution */
static
SCIP_RETCODE checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool*            violated            /**< pointer to store whether the constraint is violated */
)
{
   SCIP_CONSDATA* consdata;
   ParentSetData* psd;
   
   int i, j;                      /* node */
   int n_an = 0;                  /* how many in smallest ancestral set containing abs */
   int n_new_an = 0;              /* how many newly found to be in smallest ancestral set containing abs */
   int k;                         /* indexes parent sets */
   int l;                         /* indexes parents in a parent set */
   int i_pa;                      /* node (parent of node i) */
   int ii;                        /* index for sets containing nodes */
   int ll;                        /* indexes parents in a parent set */
   int i_pa2;                     /* node (parent of node i) */
   int n_a_open = 0;              /* length of open list of side a */
   int n_b_open = 0;              /* length of open list of side b */
   int* n_this_open_ptr;          /* pointer to length of 'current' open list */
   int* this_open;                /* 'current' open list */
   SCIP_Bool* is_this_reachable;  /* whether known to be reachable from 'current' side */
   SCIP_Bool* is_other_reachable; /* whether known to be reachable from 'other' side */
   int n_newnodes;                /* number of new nodes to add to the frontier */
   int nbr;                       /* neighbour of node just removed from open list */

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("checking ci constraint <%s> for feasibility of solution %p \n",
                    SCIPconsGetName(cons), (void*)sol);

   psd = consdata->psd;

   /* find which parent set selected in this solution (which is assumed to be integral)
   */

   assert(psd->nParentSets != NULL);
   assert(psd->PaVars != NULL);

   for( i = 0; i < psd->n; ++i )
   {
      consdata->sol_paset[i] = -1;
      assert(psd->PaVars[i] != NULL);
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         assert(psd->PaVars[i][k] != NULL);
         assert(SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k])));
         if( SCIPisGT(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k]), 0.5) )
         {
            consdata->sol_paset[i] = k;
            break;
         }
      }
      assert(consdata->sol_paset[i] > -1);
   }


   /* variables in a \cup b \up s are in the smallest ancestral subgraph of
      sol which contains a and b and s, set others to not be initally
   */

   for( i = 0; i < psd->n; ++i )
      if( consdata->in_abs[i] )
      {
         consdata->is_in_an[i] = TRUE;
         consdata->an[n_an++] = i;
         consdata->new_an[n_new_an++] = i;
      }
      else
         consdata->is_in_an[i] = FALSE;

   /* determine which (additional) variables are in the smallest ancestral subgraph of
      sol which contains a and b and s
   */

   while( n_new_an > 0 )
   {
      i = consdata->new_an[--n_new_an];
      assert(i > -1);
      assert(i < psd->n);
      k = consdata->sol_paset[i];
      assert(k > -1);
      assert(k < psd->nParentSets[i]);
      for( l = 0; l < psd->nParents[i][k]; ++l )
      {
         i_pa = psd->ParentSets[i][k][l];
         if( !consdata->is_in_an[i_pa] )
         {
            consdata->is_in_an[i_pa] = TRUE;
            consdata->an[n_an++] = i_pa;
            consdata->new_an[n_new_an++] = i_pa;
         }
      }
   }

   /* initialise empty moral graph */

   for( i = 0; i < psd->n; ++i )
   {
      consdata->n_moral_nbrs[i] = 0;
      for( j = 0; j < psd->n; ++j )
         consdata->moral[i][j] = FALSE;
   }

   /* construct moral graph */

   for( ii = 0; ii < n_an; ++ii )
   {
      i = consdata->an[ii];
      k = consdata->sol_paset[i];
      for( l = 0; l < psd->nParents[i][k]; ++l )
      {
         /* join parents to the child */
         i_pa = psd->ParentSets[i][k][l];
         SCIP_CALL( add_to_moral_graph(consdata, i, i_pa) );

         /* marry parents */
         for( ll = l + 1; ll < psd->nParents[i][k]; ++ll )
         {
            i_pa2 = psd->ParentSets[i][k][ll];
            SCIP_CALL( add_to_moral_graph(consdata, i_pa, i_pa2) );
         }
      }
   }

   /* initialise bi-directional search */

   for( i = 0; i < psd->n; ++i )
   {
      consdata->is_a_reachable[i] = FALSE;
      consdata->is_b_reachable[i] = FALSE;
   }

   for( ii = 0; ii < consdata->n_a; ++ii )
   {
      i = consdata->a[ii];
      consdata->a_open[n_a_open++] = i;
      consdata->is_a_reachable[i] = TRUE;
   }

   for( ii = 0; ii < consdata->n_b; ++ii )
   {
      i = consdata->b[ii];
      consdata->b_open[n_b_open++] = i;
      consdata->is_b_reachable[i] = TRUE;
   }

   /* start on the 'a' side */

   n_this_open_ptr = &n_a_open;
   this_open = consdata->a_open;
   is_this_reachable = consdata->is_a_reachable;
   is_other_reachable = consdata->is_b_reachable;

   /* do bi-directional search */

   while( (*n_this_open_ptr) > 0 )
   {
      n_newnodes = 0;
      i = this_open[--(*n_this_open_ptr)];
      for( ii = 0; ii < consdata->n_moral_nbrs[i]; ++ii )
      {
         nbr =  consdata->moral_nbrs[i][ii];

         if( is_other_reachable[nbr] )
            /* found a path */
         {
            *violated = TRUE;
            SCIPdebugMessage("solution %p violates ci constraint <%s>\n",
                             (void*)sol, SCIPconsGetName(cons));
            return SCIP_OKAY;
         }

         if( consdata->in_s[nbr] )
            /* nbr in separator, ignore */
            continue;

         if( !is_this_reachable[nbr] )
         {
            /* didn't previously know that nbr was reachable
               so a new node
            */
            consdata->newnodes[n_newnodes++] = nbr;
            is_this_reachable[nbr] = TRUE;
         }
      }
      /* add new nodes to the end of the open list */
      for( ii = 0; ii < n_newnodes; ++ii )
         this_open[(*n_this_open_ptr)++] = consdata->newnodes[ii];

      /* swap sides */
      if( this_open == consdata->a_open )
      {
         n_this_open_ptr = &n_b_open;
         this_open = consdata->b_open;
         is_other_reachable = consdata->is_a_reachable;
         is_this_reachable = consdata->is_b_reachable;
      }
      else
      {
         n_this_open_ptr = &n_a_open;
         this_open = consdata->a_open;
         is_other_reachable = consdata->is_b_reachable;
         is_this_reachable = consdata->is_a_reachable;
      }
   }

   *violated = FALSE;
   SCIPdebugMessage("solution %p satisfies ci constraint <%s>\n",
                    (void*)sol, SCIPconsGetName(cons));

   return SCIP_OKAY;
}


#if 0
static
SCIP_RETCODE parseSet(
   SCIP* scip,
   const char*  str,      /* string to parse */
   int*  set,             /* result set */
   int*  n_set,           /* length of set */
   SCIP_Bool* success     /* success flag */
)
{

   const char* t;
   char tmp[SCIP_MAXSTRLEN];
   int k;

   *n_set = 0;
   t = str;
   while( *t != '\0' )
   {
      k = 0;
      while( *t != '\0' && *t != ',' )
         tmp[k++] = *t++;
      tmp[k] = '\0';
      if( sscanf(tmp, "%d", &(set[(*n_set)++])) != 1 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Expected integer. Got: %s\n", tmp);
         *success = FALSE;
         return SCIP_OKAY;
      }

      if( *t == ',' )
         t++;
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyCi NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeCi)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}



/** initialization method of constraint handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_CONSINIT(consInitCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitCi NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitCi NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreCi NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreCi NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolCi NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolCi NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteCi)
{

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPdebugMessage("deleting ci constraint <%s>.\n", SCIPconsGetName(cons));

   SCIPfreeMemoryArray(scip, &((*consdata)->sol_paset));
   SCIPfreeMemoryArray(scip, &((*consdata)->is_in_an));
   SCIPfreeMemoryArray(scip, &((*consdata)->an));
   SCIPfreeMemoryArray(scip, &((*consdata)->new_an));
   SCIPfreeMemoryArray(scip, &((*consdata)->n_moral_nbrs));
   SCIPfreeMemoryArray(scip, &((*consdata)->is_a_reachable));
   SCIPfreeMemoryArray(scip, &((*consdata)->is_b_reachable));
   SCIPfreeMemoryArray(scip, &((*consdata)->a_open));
   SCIPfreeMemoryArray(scip, &((*consdata)->b_open));
   SCIPfreeMemoryArray(scip, &((*consdata)->newnodes));

   SCIPfreeMemoryArray(scip, &((*consdata)->moral));
   SCIPfreeMemoryArray(scip, &((*consdata)->moral_nbrs));
   SCIPfreeMemoryArray(scip, &((*consdata)->a));
   SCIPfreeMemoryArray(scip, &((*consdata)->b));
   SCIPfreeMemoryArray(scip, &((*consdata)->s));
   SCIPfreeMemoryArray(scip, &((*consdata)->in_abs));
   SCIPfreeMemoryArray(scip, &((*consdata)->in_a));
   SCIPfreeMemoryArray(scip, &((*consdata)->in_b));
   SCIPfreeMemoryArray(scip, &((*consdata)->in_s));

   SCIP_CALL( PS_deallocateParentSetData(scip, &((*consdata)->psd), FALSE) );


   SCIPfreeBlockMemory(scip, consdata);


   return SCIP_OKAY;
}



/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransCi NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpCi)
{

   int c;
   int nGen = 0;

   int i;
   int ii;
   int j;
   int jj;
   int k;
   int l;
   char name[SCIP_MAXSTRLEN];
   SCIP_ROW* row;
   int ia;
   int ib;

   SCIP_Bool found_i;
   SCIP_Bool found_j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* loop through all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CONSDATA* consdata;
      ParentSetData*        psd;
      
      assert(conss != NULL);
      assert(conss[c] != NULL);
      SCIPdebugMessage("adding initial rows for ci constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      psd = consdata->psd;

      /* rule out paths of length two between elements of A and B
       via j not in S */

      for( j = 0; j < psd->n; ++j )
      {
         if( consdata->in_abs[j] )
            continue;

         /* from A to B */
         for( ii = 0; ii < consdata->n_a; ++ii )
         {
            i = consdata->a[ii];
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pathcut(%d,%d)", i, j);


            SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, -SCIPinfinity(scip), 1,
                                             FALSE, FALSE, TRUE));

            found_i = FALSE;
            found_j = FALSE;

            for( k = 0; k < psd->nParentSets[i]; ++k )
               for( l = 0; l < psd->nParents[i][k]; ++l )
                  if( psd->ParentSets[i][k][l] == j )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[i][k], 1.0) );
                     found_i = TRUE;
                     break;
                  }

            for( k = 0; k < psd->nParentSets[j]; ++k )
               for( l = 0; l < psd->nParents[j][k]; ++l )
                  if( consdata->in_b[psd->ParentSets[j][k][l]] )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[j][k], 1.0) );
                     found_j = TRUE;
                     break;
                  }
            if( found_i && found_j )
            {
               SCIPdebug(SCIP_CALL( SCIPprintRow(scip, row, NULL) ));
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               ++nGen;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }


         /* from B to A */
         for( ii = 0; ii < consdata->n_b; ++ii )
         {
            i = consdata->b[ii];
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pathcut2(%d,%d)", i, j);


            SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, -SCIPinfinity(scip), 1,
                                             FALSE, FALSE, TRUE));
            found_i = FALSE;
            found_j = FALSE;


            for( k = 0; k < psd->nParentSets[i]; ++k )
               for( l = 0; l < psd->nParents[i][k]; ++l )
                  if( psd->ParentSets[i][k][l] == j )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[i][k], 1.0) );
                     found_i = TRUE;
                     break;
                  }

            for( k = 0; k < psd->nParentSets[j]; ++k )
               for( l = 0; l < psd->nParents[j][k]; ++l )
                  if( consdata->in_a[psd->ParentSets[j][k][l]] )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[j][k], 1.0) );
                     found_j = TRUE;
                     break;
                  }
            if( found_i && found_j )
            {
               SCIPdebug(SCIP_CALL( SCIPprintRow(scip, row, NULL) ));
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
               ++nGen;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }

         /* A and B with parent not in in S */
         for( ii = 0; ii < consdata->n_a; ++ii )
         {
            ia = consdata->a[ii];
            for( jj = 0; jj < consdata->n_b; ++jj )
            {
               ib = consdata->b[jj];
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "parentcut(%d,%d,%d)", ia, ib, j);
               
               SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, -SCIPinfinity(scip), 1,
                     FALSE, FALSE, TRUE));

               found_i = FALSE;
               found_j = FALSE;

               for( k = 0; k < psd->nParentSets[ia]; ++k )
                  for( l = 0; l < psd->nParents[ia][k]; ++l )
                     if( psd->ParentSets[ia][k][l] == j )
                     {
                        SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[ia][k], 1.0) );
                        found_i = TRUE;
                        break;
                     }

               for( k = 0; k < psd->nParentSets[ib]; ++k )
                  for( l = 0; l < psd->nParents[ib][k]; ++l )
                     if( psd->ParentSets[ib][k][l] == j )
                     {
                        SCIP_CALL( SCIPaddVarToRow(scip, row, psd->PaVars[ib][k], 1.0) );
                        found_j = TRUE;
                        break;
                     }
               if( found_i && found_j )
               {
                  SCIPdebug(SCIP_CALL( SCIPprintRow(scip, row, NULL) ));
                  SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
                  ++nGen;
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
         }
      }
   }

   SCIPdebugMessage("added %d inequalities.\n", nGen);

   return SCIP_OKAY;
}




/** separation method of constraint handler for LP solutions */
static SCIP_DECL_CONSSEPALP(consSepalpCi)
{
   int c;
   int nGen = 0;

   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;

   SCIP_Bool found_efficacious;

   SolutionInfo solinfo;

   int* s;
   int i;
   int n_s = 0;
   int a_i;
   int b_i;
   int must_be_included[2];

   SCIP_Bool cutoff = FALSE;

   ParentSetData* psd;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(conss != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   for( c = 0; c < nconss; ++c )
   {
      cons = conss[c];
      assert(cons != NULL);
      SCIPdebugMessage("separating LP solution for CI constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      psd = consdata->psd;
      
      assert(psd != NULL);
      assert(psd->PaVars != NULL);
      
      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);


      /* store useful information (in consdata members ) on NULL = the LP solution */
      SCIP_CALL( SI_setsolinfo(scip, &solinfo, psd, NULL, TRUE, FALSE) );

      SCIP_CALL( SCIPallocMemoryArray(scip, &s, psd->n) );
      for( i = 0; i < psd->n; ++i )
         if( consdata->in_s[i] )
            s[n_s++] = i;

      *result = SCIP_DIDNOTFIND;

      for( a_i = 0; a_i < consdata->n_a; ++a_i )
      {
         must_be_included[0] = consdata->a[a_i];
         for( b_i = 0; b_i < consdata->n_b; ++b_i )
         {
            must_be_included[1] = consdata->b[b_i];
            SCIP_CALL( IP_findCuts(scip, psd, &solinfo, NULL, &nGen, 1, 1, conshdlr, FALSE,
                  conshdlrdata->forcecuts, &found_efficacious,
                  DEFAULT_TIMELIMIT, DEFAULT_GAPLIMIT, DEFAULT_ABSGAPLIMIT, FALSE,
                  must_be_included, 2, s, n_s, TRUE, FALSE, 0, 0, &cutoff, NULL,
                  FALSE,   /* do not store ci cuts */
                  TRUE    /* only allow knapsack cuts (in case we are pricing) */
                  ) );
            
            if( cutoff )
            {
               SCIPdebugMessage("CI constraint <%s> generated cutoff in LP separator.\n", SCIPconsGetName(cons));
               *result = SCIP_CUTOFF;
               SCIPfreeMemoryArray(scip, &s);
               SI_freesolinfo(&solinfo, psd->n);
               return SCIP_OKAY;
            }
         }
      }

      SCIPfreeMemoryArray(scip, &s);

      SI_freesolinfo(&solinfo, psd->n);
   }

   if( nGen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}



/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolCi NULL
#endif



/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpCi)
{
   /*lint --e{715}*/
   SCIP_Bool violated;
   int i;
   int ncuts;

   *result = SCIP_FEASIBLE;

   SCIPdebugMessage("ci enforcement of %d constraints\n", nconss);

   ncuts = 0;
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, &violated) );
      if( violated )
      {
         SCIP_CALL( addcicut(scip, conss[i], NULL) );
         ncuts++;
      }
   }

   /* adjust the result code */
   if( ncuts > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsCi)
{
   /*lint --e{715}*/
   SCIP_Bool violated;
   int i;

   *result = SCIP_FEASIBLE;


   SCIPdebugMessage("ci pseudo solution enforcement of %d constraints\n", nconss);
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], NULL, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckCi)
{
   /*lint --e{715}*/

   SCIP_Bool violated;
   int i;

   SCIPdebugMessage("ci feasibility check of %d constraints\n", nconss);
   for( i = 0; i < nconss; i++ )
   {
      SCIP_CALL( checkCons(scip, conss[i], sol, &violated) );
      if( violated )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


#if 0
/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropCi NULL
#endif


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolCi)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   int c;
   int i;
   int k;
   int l;
   int ii;
   int jj;
   int j;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   SCIP_Bool found_a_parent;
   SCIP_Bool found_b_parent;

   ParentSetData* psd;
   SCIP_VAR* vars[5];
   SCIP_VAR* var;

   SCIP_Bool use_ancestors;
   int nvars;
   
   *result = SCIP_DIDNOTFIND;

   for( c = 0; c < nconss && !SCIPisStopped(scip); c++ )
   {
      cons = conss[c];
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      psd = consdata->psd;

      /* if non-conditional independence then can rule out ancestor relations */
      
      use_ancestors = ( consdata->n_s == 0 && consdata->ancestorvars != NULL );
      
      /* rule out edges between nodes in a and nodes in b */
      
      for( ii = 0; ii < consdata->n_a; ++ii )
      {
         i = consdata->a[ii];

         for( jj = 0; jj < consdata->n_b; ++jj )
         {
            j = consdata->b[jj];

            vars[0] = get_edge(psd,i,j);
            vars[1] = get_arrow(psd,i,j);
            vars[2] = get_arrow(psd,j,i);
            if( use_ancestors )
            {
               vars[3] = (consdata->ancestorvars)[i][j];
               vars[4] = (consdata->ancestorvars)[j][i];
               nvars = 5;
            }
            else
               nvars = 3;
            
            for( i = 0; i < nvars; ++i)
            {
               var = vars[i];
               if( var != NULL )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, var, 0, FALSE, &infeasible, &tightened) );
                  if( infeasible )
                  {
                     SCIPdebugMessage("CI constraint handler detected an infeasibility");
                     *result = SCIP_CUTOFF;
                     return SCIP_OKAY;
                  }
                  if( tightened )
                  {
                     SCIPdebugMessage("Removing variable <%s>.\n", SCIPvarGetName(var));
                     *result = SCIP_SUCCESS;
                  }
               }
            }
         }

         /* rule out (appropriate) parent set vars for children in a */
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            for( l = 0; l < psd->nParents[i][k]; ++l )
               if( consdata->in_b[psd->ParentSets[i][k][l]] )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, psd->PaVars[i][k], 0, FALSE, &infeasible, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage("Removing variable <%s>.\n", SCIPvarGetName(psd->PaVars[i][k]));
                     *result = SCIP_SUCCESS;
                  }
                  break;
               }
         }
      }

      /* rule out (appropriate) parent set vars for children in b */
      for( ii = 0; ii < consdata->n_b; ++ii )
      {
         i = consdata->b[ii];
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            for( l = 0; l < psd->nParents[i][k]; ++l )
               if( consdata->in_a[psd->ParentSets[i][k][l]] )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, psd->PaVars[i][k], 0, FALSE, &infeasible, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage("Removing variable <%s>.\n", SCIPvarGetName(psd->PaVars[i][k]));
                     *result = SCIP_SUCCESS;
                  }
                  break;
               }
         }
      }
            
      for( i = 0; i < psd->n; ++i )
      {
         if( !consdata->in_s[i] )
            continue;

         /* a child in s cannot have a parent in a and a parent in b */
         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            found_a_parent = FALSE;
            found_b_parent = FALSE;
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               if( consdata->in_a[psd->ParentSets[i][k][l]] )
                  found_a_parent = TRUE;
               else if( consdata->in_b[psd->ParentSets[i][k][l]] )
                  found_b_parent = TRUE;

               if( found_a_parent && found_b_parent )
               {
                  SCIP_CALL( SCIPtightenVarUb(scip, psd->PaVars[i][k], 0, FALSE, &infeasible, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage("Removing variable <%s>.\n", SCIPvarGetName(psd->PaVars[i][k]));
                     *result = SCIP_SUCCESS;
                     break;
                  }
               }
            }
         }
      }
         
   }
   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropCi NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockCi)
{
   /*lint --e{715}*/

   SCIP_CONSDATA* consdata;
   int i;
   int j;
   SCIP_VAR* arrow_i_j;
   
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->psd->arrow != NULL);

   /* only lock arrow variables.
      a conditional independence constraint can only be violated by adding, not removing, an arrow
      hence "nlocksneg, nlockspos" below
   */
   
   for( i = 0; i < consdata->psd->n; ++i )
      for( j = 0; j < consdata->psd->n; ++j )
         if( i != j )
         {
            arrow_i_j = get_arrow(consdata->psd, i, j);
            if( arrow_i_j != NULL )
               SCIP_CALL( SCIPaddVarLocks(scip, arrow_i_j, nlocksneg, nlockspos) );
         }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveCi NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveCi NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableCi NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableCi NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsCi NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintCi)
{

   SCIP_CONSDATA* consdata;
   int i;
   int ii;
   int k;
   int l;
   SCIP_Bool first;
   ParentSetData* psd;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   psd = consdata->psd;
   
   SCIPdebugMessage("Printing method for ci constraint handler\n");

   SCIPinfoMessage(scip, file, "ci(");
   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         if( i > 0 || k > 0 )
            SCIPinfoMessage(scip, file, ",");
         SCIPinfoMessage(scip, file, "%s", SCIPvarGetName(psd->PaVars[i][k]));
      }

   SCIPinfoMessage(scip, file, "|");

   for( i = 0; i < psd->n; ++i )
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {
         if( i > 0 || k > 0 )
            SCIPinfoMessage(scip, file, ",");
         SCIPinfoMessage(scip, file, "%d<-{", i);
         for( l = 0; l < psd->nParents[i][k]; ++l )
         {
            if( l > 0 )
               SCIPinfoMessage(scip, file, ",");
            SCIPinfoMessage(scip, file, "%d", psd->ParentSets[i][k][l]);
         }
         SCIPinfoMessage(scip, file, "}");
      }

   SCIPinfoMessage(scip, file, ",{");
   for( ii = 0; ii < consdata->n_a; ++ii )
   {
      if( ii > 0 )
         SCIPinfoMessage(scip, file, ",");
      SCIPinfoMessage(scip, file, "%d", consdata->a[ii]);
   }

   SCIPinfoMessage(scip, file, "},{");

   for( ii = 0; ii < consdata->n_b; ++ii )
   {
      if( ii > 0 )
         SCIPinfoMessage(scip, file, ",");
      SCIPinfoMessage(scip, file, "%d", consdata->b[ii]);
   }

   SCIPinfoMessage(scip, file, "},{");

   first = TRUE;
   for( i = 0; i < psd->n; ++i )
   {
      if( !consdata->in_s[i] )
         continue;

      if( !first )
         SCIPinfoMessage(scip, file, ",");
      SCIPinfoMessage(scip, file, "%d", i);
      first = FALSE;
   }
   SCIPinfoMessage(scip, file, "})");

   return SCIP_OKAY;


}



/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyCi NULL
#endif


/* this turned off temporarily */
#if 0
/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseCi)
{
   const char* s;
   const char* start;
   char tmp[SCIP_MAXSTRLEN];
   char varname[SCIP_MAXSTRLEN];
   char varinfo[SCIP_MAXSTRLEN];
   char a_str[SCIP_MAXSTRLEN];
   char b_str[SCIP_MAXSTRLEN];
   char s_str[SCIP_MAXSTRLEN];
   int nvars;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int* ch;
   int* npa;
   int** pa;

   int k;
   int nv;
   int tmp_pa[SCIP_MAXSTRLEN];
   int varindex;

   int* a;
   int* b;
   int* sep;
   int n_a;
   int n_b;
   int n_s;

   /* adapted from CONSPARSE in cons_dagcluster */

   assert(success != NULL);

   *success = TRUE;
   s = str;

   /* skip white space */
   while( *s != '\0' && isspace((unsigned char)*s) )
      ++s;

   if(strncmp(s, "ci(", 3) != 0)
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - expected \"ci\": %s\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }

   s += 3;
   start = s;

   /* check that number of variables and number of parent sets are equal and record this number */
   nvars = 0;
   while( *s != '\0' && ! isspace((unsigned char)*s)  && *s != '|' && *s != '.' && *s != ')' )
   {
      if( *s == ',' )
         nvars++;
      s++;
   }

   if( *s != '|' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - missing \"|\": %s\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
      nvars++;

   s++;
   nv = 0;
   while( *s != '\0' && ! isspace((unsigned char)*s)  && *s != ';' && *s != '.' && *s != ')' )
   {
      if( *s == '}' )
         nv++;
      s++;
   }

   if( *s != ')' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - missing \")\": %s\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* should get an extra 3 '}' due to A, B and S */
   if( nvars != nv - 3 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - number of variables (%d) and items of information on variables (%d) are not equal: %s\n", nvars, nv - 3, s);
      *success = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &ch, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &npa, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pa, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &a, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &b, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &sep, nvars) );

   /* rewind */
   s = start;

   /* get variables */
   for( varindex = 0; varindex < nvars; ++varindex )
   {
      k = 0;
      while( *s != '\0' && ! isspace((unsigned char)*s)  && *s != '|' && *s != ',' && *s != '.' && *s != ')' )
         varname[k++] = *s++;
      varname[k] = '\0';

      /* get variable */
      var = SCIPfindVar(scip, varname);
      if( var == NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable <%s>\n", varname);
         *success = FALSE;
         return SCIP_OKAY;
      }
      vars[varindex] = var;

      if( varindex == nvars - 1 && *s != '|' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - expected \"|\": %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      else if( varindex < nvars - 1 && *s != ',' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - expected \",\": %s\n", s);
         *success = FALSE;
         return SCIP_OKAY;
      }
      else
         s++;
   }

   /* skip white space */
   while( *s != '\0' && isspace((unsigned char)*s) )
      s++;

   /* get variable info */
   for( varindex = 0; varindex < nvars; ++varindex )
   {
      k = 0;
      while( *s != '}' )
         varinfo[k++] = *s++;
      varinfo[k++] = *s++;
      varinfo[k] = '\0';

      npa[varindex] = 0;
      if( strcmp(varinfo + k - 4, "<-{}") == 0 )
      {
         if( sscanf(varinfo, "%d", &(ch[varindex])) != 1 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - could not get child with no parents: \"%s\"\n", varinfo);
            *success = FALSE;
            return SCIP_OKAY;
         }
         SCIP_CALL( SCIPallocMemoryArray(scip, &(pa[varindex]), 0) );
      }
      else
      {
         if( sscanf(varinfo, "%d<-{%[^}]}", &(ch[varindex]), tmp) != 2 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - when parsing varinfo: \"%s\"\n", varinfo);
            *success = FALSE;
            return SCIP_OKAY;
         }

         SCIP_CALL( parseSet(scip, tmp, tmp_pa, &(npa[varindex]), success) );
         if( !(*success) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - when parsing set:%s\n", tmp);
            *success = FALSE;
            return SCIP_OKAY;
         }

         SCIP_CALL( SCIPallocMemoryArray(scip, &(pa[varindex]), npa[varindex]) );
         for( k = 0; k < npa[varindex]; ++k )
            pa[varindex][k] = tmp_pa[k];
      }

      if( varindex < nvars - 1 )
      {
         if( *s != ',' )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - expected \",\": %s\n", s);
            *success = FALSE;
            return SCIP_OKAY;
         }
         else
            s++;
      }
   }

   /* now parse A, B and S sets */

   if( sscanf(s, ",{%[^}]},{%[^}]}", a_str, b_str) != 2 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - could not find \"A\" and \"B\" sets:%s\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( parseSet(scip, a_str, a, &n_a, success) );
   if( !(*success) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - could not parse \"A\" set:%s\n", a_str);
      *success = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( parseSet(scip, b_str, b, &n_b, success) );
   if( !(*success) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - could not parse \"B\" set:%s\n", b_str);
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* get to S set */
   while( *s != '}' )
      s++;
   s++;
   while( *s != '}' )
      s++;
   s++;

   if( strcmp(s, ",{})" ) == 0)
      n_s = 0;
   else if( sscanf(s, ",{%[^}]})", s_str ) != 1)
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - could not find \"S\" set: %s\n", s);
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
   {
      SCIP_CALL( parseSet(scip, s_str, sep, &n_s, success) );
      if( !(*success) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error - could not parse \"S\" set:%s\n", s_str);
         *success = FALSE;
         return SCIP_OKAY;
      }
   }

   /* create and populate parent set data struct */


   SCIP_CALL(SCIPcreateConsCi(scip, cons, name, nvars, vars, ch, npa, pa, a, n_a, b, n_b, sep, n_s, initial, separate,
                              enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode));

   for( varindex = 0; varindex < nvars; ++varindex )
      SCIPfreeMemoryArray(scip, &(pa[varindex]));

   SCIPfreeMemoryArray(scip, &vars);
   SCIPfreeMemoryArray(scip, &ch);
   SCIPfreeMemoryArray(scip, &npa);
   SCIPfreeMemoryArray(scip, &pa);
   SCIPfreeMemoryArray(scip, &a);
   SCIPfreeMemoryArray(scip, &b);
   SCIPfreeMemoryArray(scip, &sep);


   return SCIP_OKAY;
}
#else
#define ConsParseCi NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsCi NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsCi)
{
   /*lint --e{715}*/
   SCIPerrorMessage("method of ci power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsCi NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for ci constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrCi(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   conshdlr = NULL;
   conshdlrdata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */

   SCIP_CALL(SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
                                      CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
                                      consEnfolpCi, consEnfopsCi, consCheckCi, consLockCi,
                                      conshdlrdata));
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   /* SCIP_CALL(  SCIPsetConshdlrActive(scip, conshdlr, consActiveCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyCi, consCopyCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveCi)  ); */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteCi) );
   /* SCIP_CALL(  SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrDisable(scip, conshdlr, consDisableCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrEnable(scip, conshdlr, consEnableCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrExit(scip, conshdlr, consExitCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolCi)  ); */
    SCIP_CALL(  SCIPsetConshdlrFree(scip, conshdlr, consFreeCi)  ); 
   /* SCIP_CALL(  SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrInit(scip, conshdlr, consInitCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreCi)  ); */
   /* SCIP_CALL(  SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolCi)  ); */
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpCi) );
   /*SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseCi) );*/
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolCi, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintCi) );
   /* SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropCi, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, */
   /*       CONSHDLR_PROP_TIMING) ); */
   /* SCIP_CALL(  SCIPsetConshdlrResprop(scip, conshdlr, consRespropCi)  ); */
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpCi, consSepasolCi, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   /* SCIP_CALL(  SCIPsetConshdlrTrans(scip, conshdlr, consTransCi)  ); */

   /* add ci constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

   SCIP_CALL(SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/forcecuts",
         "whether to force all cuts to be added",
         &conshdlrdata->forcecuts, TRUE, DEFAULT_FORCECUTS, NULL, NULL));


   return SCIP_OKAY;
}

/** creates and captures a conditional independence constraint \f$ A \perp B | S\f$
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsCi(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,                /**< parent set data containing family, edge and arrow variables */
   SCIP_VAR***           ancestorvars,       /**< anectorvars (or NULL if absent) */
   int*                  a,                  /**< the set \f$A\f$ */
   int                   n_a,                /**< size of set \f$A\f$ */
   int*                  b,                  /**< the set \f$B\f$ */
   int                   n_b,                /**< size of set \f$B\f$ */
   int*                  s,                  /**< the set \f$S\f$ */
   int                   n_s,                /**< size of set \f$S\f$ */
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsCi() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   int i;
   int ii;
   int n;


   if( n_a == 0 )
   {
      SCIPerrorMessage("ci constraint %s: trying to post conditional independence constraint with an empty \"A\" set\n", name);
      return SCIP_ERROR;
   }

   if( n_b == 0 )
   {
      SCIPerrorMessage("ci constraint %s:trying to post conditional independence constraint with an empty \"B\" set\n", name);
      return SCIP_ERROR;
   }

   /* find the ci constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("ci constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create ci constraint handler data */
   consdata = NULL;
   /* TODO: (optional) create constraint handler specific data here */

   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   SCIP_CALL( PS_copyParentSetData(scip, psd, &(consdata->psd)) );

   n = psd->n;

   consdata->ancestorvars = ancestorvars;
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->sol_paset, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->is_in_an, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->an, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->new_an, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->n_moral_nbrs, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->is_a_reachable, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->is_b_reachable, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->a_open, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->b_open, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->newnodes, n) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->moral, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->moral_nbrs, n) );
   for( i = 0; i < n; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(consdata->moral[i]), n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(consdata->moral_nbrs[i]), n) );
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->a, n_a) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->b, n_b) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->s, n_s) );
   consdata->n_a = n_a;
   consdata->n_b = n_b;
   consdata->n_s = n_s;

   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->in_abs, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->in_a, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->in_b, n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &consdata->in_s, n) );

   for( i = 0; i < n; ++i )
   {
      consdata->in_abs[i] = FALSE;
      consdata->in_a[i] = FALSE;
      consdata->in_b[i] = FALSE;
      consdata->in_s[i] = FALSE;
   }

   for( ii = 0; ii < n_a; ++ii )
   {
      i = a[ii];
      consdata->a[ii] = i;
      if( consdata->in_a[i] )
      {
         SCIPerrorMessage("ci constraint %s: %d occurs twice in \"A\" set\n", name, i);
         return SCIP_ERROR;
      }
      consdata->in_a[i] = TRUE;
      consdata->in_abs[i] = TRUE;
   }

   for( ii = 0; ii < n_b; ++ii )
   {
      i = b[ii];
      if( consdata->in_b[i] )
      {
         SCIPerrorMessage("ci constraint %s: %d occurs twice in \"B\" set\n", name, i);
         return SCIP_ERROR;
      }
      consdata->b[ii] = i;
      consdata->in_b[i] = TRUE;
      consdata->in_abs[i] = TRUE;
   }
   for( ii = 0; ii < n_s; ++ii )
   {
      i = s[ii];
      if( consdata->in_s[i] )
      {
         SCIPerrorMessage("ci constraint %s: %d occurs twice in \"S\" set\n", name, i);
         return SCIP_ERROR;
      }
      consdata->s[ii] = i;
      consdata->in_s[i] = TRUE;
      consdata->in_abs[i] = TRUE;
   }

   for( i = 0; i < n; ++i )
   {
      if( consdata->in_a[i] && consdata->in_b[i] )
      {
         SCIPerrorMessage("ci constraint %s: Can't have %d in both \"A\" and \"B\" set\n", name, i);
         return SCIP_ERROR;
      }
      if( consdata->in_a[i] && consdata->in_s[i] )
      {
         SCIPerrorMessage("ci constraint %s: Can't have %d in both \"A\" and \"S\" set\n", name, i);
         return SCIP_ERROR;
      }
      if( consdata->in_b[i] && consdata->in_s[i] )
      {
         SCIPerrorMessage("ci constraint %s: Can't have %d in both \"B\" and \"S\" set\n", name, i);
         return SCIP_ERROR;
      }
   }


   /* create constraint */
   SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
                            local, modifiable, dynamic, removable, stickingatnode));

   return SCIP_OKAY;
}

/** creates and captures a ci constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
/* SCIP_RETCODE SCIPcreateConsBasicCi( */
/*    SCIP*                 scip,               /\**< SCIP data structure *\/ */
/*    SCIP_CONS**           cons,               /\**< pointer to hold the created constraint *\/ */
/*    const char*           name,               /\**< name of constraint *\/ */
/*    int                   nvars,              /\**< number of variables in the constraint *\/ */
/*    SCIP_VAR**            vars,               /\**< array with variables of constraint entries *\/ */
/*    SCIP_Real*            coefs,              /\**< array with coefficients of constraint entries *\/ */
/*    SCIP_Real             lhs,                /\**< left hand side of constraint *\/ */
/*    SCIP_Real             rhs                 /\**< right hand side of constraint *\/ */
/*    ) */
/* { */
/*    SCIP_CALL( SCIPcreateConsCi(scip, cons, name, nvars, vars, coefs, lhs, rhs, */
/*          TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) ); */

/*    return SCIP_OKAY; */
/* } */
