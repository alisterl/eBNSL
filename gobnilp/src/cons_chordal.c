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

/**@file   cons_chordal.c
 * @brief  constraint handler for chordal graph constraints
 * @author James Cussens
 *
 * Implements a constraint that the graph has no immoralities (v-structures).
 * This is equivalent to demanding that the acyclic digraph is equivalent to a 
 * chordal undirected graph (aka decomposable model)
 *
 *  Doxygen documentation for the locally defined structs @c SCIP_ConsData and @c SCIP_ConshdlrData is unfortunately not available
 *  since structs with the same name are defined in other constraint handlers and Doxygen cannot handle the name clash. You will have to consult
 *  the source code to view this documentation.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*#define SCIP_DEBUG*/
#include <assert.h>
#include <string.h>
#include "cons_chordal.h"
#include "utils.h"
#include "scip/pub_misc.h"
#include "scip/scipdefplugins.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "chordal" /**< constraint handler name */
#define CONSHDLR_DESC          "constraint handler for chordal graph constraints" /**< constraint handler description */
#define CONSHDLR_ENFOPRIORITY         -10 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        -10 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        TRUE /**< should separation method be delayed, if other separators found cuts? */
 /* current propagator too slow, so turned off, by default */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP/**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_MEDIUM /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */

#define DEFAULT_MONOLIM 3           /**< default value for limit on size of subset corresponding to the lower bound in a simple monotonicity constraint */
#define DEFAULT_VANILLA FALSE       /**< default value for whether all chordal constraints are for 'vanilla' learning */
#define DEFAULT_ONESINGLETON TRUE   /**< default value for whether clutters for cuts should contain only one singleton */
#define DEFAULT_MAXCUTSROUND 500    /**< default value for maximum number of cuts generated in a separation round (-1 is no limit) */
#define DEFAULT_MAXCUTSROUNDSINGLETON 200 /**< default value for maximum number of cuts generated in a separation round for any singleton (-1 is no limit) */
#define DEFAULT_MINCLUTTERSIZE 3    /**< default value for minimal size of clutter for a cut */
#define DEFAULT_MAXCLUTTERSIZE 4    /**< default value for maximal size of clutter for a cut */
#define DEFAULT_MAXSIZEINCLUTTER -1    /**< default value for maximal size of a subset in a clutter for a cut (-1 is no limit) */

#define min(A,B) ((A) > (B) ? (B) : (A))
#define max(A,B) ((A) < (B) ? (B) : (A))

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
/** tree of integers 
 * used for mapping sequences of integers (ie sets) to integers 
 */
struct tree
{
   int idx;                /**< integer associated with this node */
   int nchildren;          /**< number of children of this node */
   struct tree* children;  /**< children of this node */
};
typedef struct tree TREE;


/** constraint data for chordal constraints 
 */
struct SCIP_ConsData
{
   ParentSetData* psd;      /**< input parent set data */
   int* n_imsetvars;        /**< @c n_imsetvars[size] is the number of imsetvars for components of size @c size */
   SCIP_VAR*** imsetvars;   /**< @c imsetvars[size][i] is the ith imsetvar of size @c size */
   int*** subsets;          /**< @c subsets[size][i] is the subset associated with the ith imsetvar of size @c size */
   SCIP_Bool**** is_subset; /**< @c is_subset[size][ci][size2][cj] if first is a subset of the other (where @c size < @c size2) */
   TREE* tree;              /**< reverse of subsets */
   SCIP_VAR**** pavars;     /**< pavars[size][ci] is the set of pavars whose sum = imsetvars[size][ci] */
   int** n_pavars;          /**< n_pavars[size][ci] is the number of pavars whose sum = imsetvars[size][ci] */
};

/** constraint handler data for chordal constraints */
struct SCIP_ConshdlrData
{
   int monolim;                     /**< limit on size of subset corresponding to the lower bound in a simple monotonicity constraint */
   SCIP_Bool vanilla;               /**< whether all chordal constraints are for 'vanilla' learning */
   SCIP_Bool onesingleton;          /**< whether clutters for cuts should contain only one singleton */
   int maxcutsround;                /**< maximum number of cuts generated in a separation round (-1 is no limit) */
   int maxcutsroundsingleton;       /**< maximum number of cuts generated in a separation round for any singleton (-1 is no limit) */
   int mincluttersize;              /**< minimal size of clutter for a cut */
   int maxcluttersize;              /**< maximal size of clutter for a cut */
   int maxsizeinclutter;           /**< maximal size of subset in a clutter for a cut */
};



/*
 * Local methods
 */

/** free memory taken up by children of a node */
static
void free_children(
   TREE* tree         /**< (pointer to) node whose children are to be deleted */
   )
{
   int i;
   TREE* child;

   assert( tree != NULL );

   for( i = 0; i < tree->nchildren; i++ )
   {
      child = (tree->children) + i;
      if( child->children != NULL )
         free_children(child);
   }
   free(tree->children);
}

/** return the integer ('index') associated with a set of integers  
 * @return the index associated with the given set 
*/
static
int get_idx(
   TREE* tree,       /**< tree containing mapping from sets to integers */ 
   int size,         /**< size of the set */
   int* subset       /**< set whose index is sought */ 
   )
{
   int i;
   int elt;

   assert( tree != NULL );

   for( i = 0; i < size; i++ )
   {
      elt = subset[i];

      if( elt >= tree->nchildren )
         return -1;

      tree = (tree->children) + elt;

      if( tree == NULL )
         return -1;
   }
   return tree->idx;
}
/** associate an integer ('index') with a set of integers  
*/
static
int put_index(
   TREE* tree,    /**< tree containing mapping from sets to integers */ 
   int size,      /**< size of the set */
   int* subset,   /**< set whose index is sought */ 
   int idx        /**< index to associate with the set */
   )
{
   int i;
   int j;
   int elt;
   int nelts;
   TREE* child;

   assert( tree != NULL );

   for( i = 0; i < size; i++ )
   {
      elt = subset[i];
      nelts = elt + 1;

      /* make new child trees if necessary */
      if ( nelts > tree->nchildren )
      {
         if( tree->nchildren == 0 )
            tree->children = (TREE *) malloc(nelts * sizeof(TREE));
         else
            tree->children = (TREE *) realloc(tree->children, nelts * sizeof(TREE));
   
         if( tree->children == NULL )
            /* could not get memory, return error indicator */
            return -1;

         /* initialise new children */
         for( j = tree->nchildren; j < nelts; j++ )
         {
            child = (tree->children)+j;
            assert( child != NULL );
            child->idx = -1;
            child->nchildren = 0;
            child->children = NULL;
         }

         /* update number of children for this tree */
         tree->nchildren = nelts;
      }

      /* move down the tree */
      tree = (tree->children) + elt;
      
      assert( tree != NULL );
   }
   tree->idx = idx;
   return idx;
}

/** add an imset variable to an existing clutter cut 
 * @return value (ie LHS) of cut for solution being separated (has to exceed RHS for cut to be useful) 
*/
static
SCIP_Real update_cut(
   int size,               /**< @c size where @c imsetvars[size][ci] is the variable to add to the cut */
   int ci,                 /**< @c ci where @c imsetvars[size][ci] is the variable to add to the cut */
   SCIP_Real** lpvals,     /**< values of the @c imsetvars[size][ci] in solution to be separated */
   SCIP_VAR*** imsetvars,  /**< @c imsetvars[size][ci] is the cith imsetvar of size @c size */
   SCIP_VAR** cutvars,     /**< output: variables in the cut */
   SCIP_Real* cutvals,     /**< output: values for variables in the cut */
   int* n_cutvars_ptr,     /**< output: number of variables in the cut */
   int* cut_rhs,           /**< output: cut RHS (only changes if imset is of size 1) */ 
   int val                 /**< coefficient for added imset variable in the cut */
   )
   
{
   assert( lpvals != NULL);
   assert( lpvals[size] != NULL);
   assert( imsetvars != NULL);
   assert( cutvars != NULL);
   assert( cutvals != NULL);
   assert( n_cutvars_ptr != NULL);
   assert( cut_rhs != NULL);

   if( ci == -1 )
      return 0;

   if( size == 1)
   {
      (*cut_rhs)++;
      return 0;
   }

   /* printf("f %d %d \n", size, ci); */

   assert( imsetvars[size] != NULL);   
   assert( imsetvars[size][ci] != NULL);   

   cutvars[*n_cutvars_ptr] = imsetvars[size][ci];
   cutvals[(*n_cutvars_ptr)++] = val;
   return val*lpvals[size][ci];
}

/** create a set of integers representation where the size is recorded in element 0 
 * space for the output set must already have been allocated!  
 */
static
void init_union(
   int                   size,  /**< size of input set */
   int*                  set,   /**< input set */
   int*                  uni    /**< output set where @c size is recorded at @c uni[0] */
   )
{
   int j;

   assert( size > -1 );
   assert( set != NULL );

   for(j = 0; j < size; j++)
      uni[j+1] = set[j];
   uni[0] = size;
}

/** compute union of two sets 
 * a set is an ordered sequence of distinct integers
 * space for the output set must already have been allocated!  
 */
static
void get_union(
   int*                  set1,  /**< first input set */
   int*                  set2,  /**< second input set */
   int*                  uni    /**< union of @c set1 and @c set2 */
   )
{
   int i1 = 1;
   int i2 = 1;
   int iuni = 1;
   int size1;
   int size2;

   assert(set1 != NULL);
   assert(set2 != NULL);

   size1 = set1[0];
   size2 = set2[0];
   
   assert( size1 > -1 );
   assert( size2 > -1 );

   while( i1 <= size1 || i2 <= size2 )
   {
      if( i1 > size1 )
      {
         uni[iuni++] = set2[i2++];
         continue;
      }

      if( i2 > size2 )
      {
         uni[iuni++] = set1[i1++];
         continue;
      }

      if( set1[i1] <= set2[i2] )
      {
         if( set1[i1] == set2[i2] )
            i2++;
         uni[iuni++] = set1[i1++];
      }
      else 
      {
         uni[iuni++] = set2[i2++];
      }
   }
   uni[0] = iuni - 1;
}

/* static */
/* void getName( */
/*    char* s, */
/*    int size, */
/*    int ci, */
/*    SCIP_VAR*** imsetvars */
/*    ) */
/* { */
/*    if( size == 1 ) */
/*       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%d",  ci); */
/*    else */
/*       (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "%s",  SCIPvarGetName(imsetvars[size][ci])); */
/* } */

/** brute force approach to separate a solution using clutter cuts 
 *  @return SCIP_OKAY if successful, or an appropriate error otherwise.
 */ 
static
SCIP_RETCODE ChordalSeparate(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_SOL*             sol,                /**< solution to be separated */
   int*                  nGen,               /**< output: pointer to store number of added rows */
   SCIP_Bool*            cutoff              /**< output: pointer to store whether we detected a cutoff */
   )
{

   SCIP_CONSHDLRDATA* conshdlrdata;
   int size;
   int n;
   int nposvars = 0;

   int* n_imsetvars;
   int*** subsets;
   SCIP_Bool**** is_subset;
   SCIP_VAR*** imsetvars;
   TREE* tree;

   int i0;
   int i1;
   int i2;
   int i3;

   int ci0;
   int ci1;
   int ci2;
   int ci3;
   int ci01;
   int ci02;
   int ci03;
   int ci12;
   int ci13;
   int ci23;
   int ci012;
   int ci013;
   int ci023;
   int ci123;
   int ci0123;

   int size0;
   int size1;
   int size2;
   int size3;
   int size01;
   int size02;
   int size03;
   int size12;
   int size13;
   int size23;
   int size012;
   int size013;
   int size023;
   int size123;
   int size0123;

   int union0[2];
   int* union1;
   int* union2;
   int* union3;
   int* union01;
   int* union02;
   int* union03;
   int* union12;
   int* union13;
   int* union23;
   int* union012;
   int* union013;
   int* union023;
   int* union123;
   int* union0123;


   SCIP_VAR* cutvars[14];
   SCIP_Real cutvals[14];
   int n_cutvars;

   int cut_rhs;

   SCIP_Real lhs1;
   SCIP_Real lhs2;
   SCIP_Real lhs3;

   int n_cutvars1;
   int cut_rhs1;
   int n_cutvars2;
   int cut_rhs2;

   int* sizes;
   int* cis;
   SCIP_Real** lpvals;
   int ci;

   SCIP_ROW* cut;   

   /* char s0[SCIP_MAXSTRLEN]; */
   /* char s1[SCIP_MAXSTRLEN]; */
   /* char s2[SCIP_MAXSTRLEN]; */

   /* int t; */
   int total;

   int maxcutsround;
   int maxcutsroundsingleton;
   SCIP_Bool onesingleton;
   int mincluttersize;
   int maxcluttersize;
   int maxsizeinclutter;

   int n_cutsi0;
   
   int biggest_size;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( consdata->psd != NULL );
   assert( nGen != NULL );
   assert( cutoff != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   onesingleton = conshdlrdata->onesingleton;
   maxcutsround = conshdlrdata->maxcutsround;
   maxcutsroundsingleton = conshdlrdata->maxcutsroundsingleton;
   mincluttersize = conshdlrdata->mincluttersize;
   maxcluttersize = conshdlrdata->maxcluttersize;
   maxsizeinclutter = conshdlrdata->maxsizeinclutter;

   if( maxcluttersize < min(2,mincluttersize) )
      return SCIP_OKAY;

   n = consdata->psd->n;
   is_subset = consdata->is_subset;
   n_imsetvars = consdata->n_imsetvars;
   subsets = consdata->subsets;
   imsetvars = consdata->imsetvars;
   tree = consdata->tree;

   if( maxsizeinclutter == -1 )
      maxsizeinclutter = n;

   total = 0;
   for( size = 1; size <= n; size++)
      total += n_imsetvars[size];

   SCIP_CALL( SCIPallocMemoryArray(scip, &sizes, total) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cis, total) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &lpvals, n+1) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &(lpvals[1]),n_imsetvars[1]) );
   for( ci = 0; ci < n_imsetvars[1]; ci++ )
   {
      lpvals[1][ci] = 1;
      sizes[nposvars] = 1;
      cis[nposvars++] = ci;
   }


   biggest_size = 0;
   for( size = 2; size <= n; size++)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(lpvals[size]),n_imsetvars[size]) );
      for( ci = 0; ci < n_imsetvars[size]; ci++ )
      {
         lpvals[size][ci] = SCIPgetSolVal(scip, sol, imsetvars[size][ci]);
         if( SCIPisPositive(scip, lpvals[size][ci]) )
         {
            sizes[nposvars] = size;
            cis[nposvars++] = ci;
            biggest_size = max(biggest_size,size);
         }
      }
   }     

   SCIP_CALL( SCIPallocMemoryArray(scip, &union1, biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union2, biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union3, biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union01, biggest_size+2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union02, biggest_size+2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union03, biggest_size+2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union12, 2*biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union13, 2*biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union23, 2*biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union012, 2*biggest_size+2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union013, 2*biggest_size+2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union023, 2*biggest_size+2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union123, 3*biggest_size+1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &union0123, 3*biggest_size+2) );


   /* brute force search for clutters up to size 4 */ 

   assert( n <= nposvars );
   for(i0 = 0; i0 < n; i0++)
   {

      /* i0 is the singleton in the clutter */
      size0 = 1;
      ci0 = i0;

      union0[0] = 1;
      union0[1] = i0;

      /* initialise cut */
      cut_rhs = 0;
      n_cutvars = 0;
      n_cutsi0 = 0;
      /* don't need lhs0 = 0; */

      /* if only one singleton allowed 
         skip past other singletons
      */
      for(i1 = ( onesingleton ? n : i0 +1 ); i1 < nposvars; i1++)
      {
         size1 = sizes[i1];
         
         if( size1 > maxsizeinclutter )
            continue;

         ci1 = cis[i1];
         
         if( size0 < size1 && is_subset[size0][ci0][size1][ci1] )
            continue;

         lhs1 = update_cut(size1, ci1, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);
  
         /* lhs1 = update_bounds(size1, ci1, lpvals, &cut_rhs, -1); */
         
         init_union(size1, subsets[size1][ci1], union1);
         
         /* now deal with 'interaction terms (only 1) */
         get_union(union0,union1,union01);

         /* printf("union0: "); */
         /* for( t = 0; t < union0[0]; t++) */
         /*    printf("%d,",union0[t+1]); */
         /* printf("\n"); */

         
         /* printf("union1: "); */
         /* for( t = 0; t < union1[0]; t++) */
         /*    printf("%d,",union1[t+1]); */
         /* printf("\n"); */


         /* printf("union01: "); */
         /* for( t = 0; t < union01[0]; t++) */
         /*    printf("%d,",union01[t+1]); */
         /* printf("\n"); */
         
         size01 = union01[0];
         ci01 = get_idx(tree, size01, union01+1);

         lhs1 += update_cut(size01, ci01, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, 1);

         if( lhs1 > cut_rhs && mincluttersize >= 2  )
         {
            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "cluttercut", -SCIPinfinity(scip), cut_rhs, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPaddVarsToRow(scip, cut, n_cutvars, cutvars, cutvals) );

            /* SCIPdebugMessage(" -> Clutter-cut <cluttercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n", */
            /*    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut), */
            /*    SCIPgetCutEfficacy(scip, NULL, cut), */
            /*    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut), */
            /*    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut)); */
            /* SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));   */
               

            SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
            if( *cutoff )
            {
               SCIPdebugMessage("Clutter cut led to cutoff\n");
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               goto TERMINATE;
            }
            SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            (*nGen)++;
            n_cutsi0++;
            if( n_cutsi0 > maxcutsroundsingleton && maxcutsroundsingleton != -1 )
               goto TERMINATE;
            if( *nGen > maxcutsround && maxcutsround != -1 )
               goto TERMINATE;
         }

         n_cutvars1 = n_cutvars;
         cut_rhs1 = cut_rhs;

         if( maxcluttersize < 3 )
         {
            cut_rhs = 0;
            n_cutvars = 0;
            continue;
         }         
         
         for(i2 = i1 + 1; i2 < nposvars; i2++)
         {

            size2 = sizes[i2];

            if( size2 > maxsizeinclutter )
               continue;
         
            ci2 = cis[i2];
            /* only clutters */
            if( (size0 < size2 && is_subset[size0][ci0][size2][ci2]) 
               || (size1 < size2 && is_subset[size1][ci1][size2][ci2]) )
               continue;
            
            lhs2 = update_cut(size2, ci2, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);

            init_union(size2, subsets[size2][ci2], union2);
            
            /* now deal with 'interaction terms (only 3) */
            
            get_union(union0,union2,union02);
            get_union(union1,union2,union12);
            get_union(union01,union2,union012);

            size02 = union02[0];
            size12 = union12[0];
            size012 = union012[0];

            /* printf("union2: "); */
            /* for( t = 0; t < size2; t++) */
            /*    printf("%d,",union2[t+1]); */
            /* printf("\n"); */

            /* printf("union02: "); */
            /* for( t = 0; t < size02; t++) */
            /*    printf("%d,",union02[t+1]); */
            /* printf("\n"); */

            /* printf("union12: "); */
            /* for( t = 0; t < size12; t++) */
            /*    printf("%d,",union12[t+1]); */
            /* printf("\n"); */

            /* printf("union012: "); */
            /* for( t = 0; t < size012; t++) */
            /*    printf("%d,",union012[t+1]); */
            /* printf("\n"); */
         

   
            ci02 = get_idx(tree, size02, union02+1);
            ci12 = get_idx(tree, size12, union12+1);
            ci012 = get_idx(tree, size012, union012+1);
           
            lhs2 += update_cut(size02, ci02, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, 1);
            lhs2 += update_cut(size12, ci12, lpvals, imsetvars, cutvars, cutvals,&n_cutvars, &cut_rhs, 1);
            lhs2 += update_cut(size012, ci012, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);
            
            if( lhs1 + lhs2 > cut_rhs && mincluttersize >= 3 )
            {
               /* got a cut */
               /* getName(s0, size0, ci0, imsetvars); */
               /* getName(s1, size1, ci1, imsetvars); */
               /* getName(s2, size2, ci2, imsetvars); */
               /* printf("clutter: %s %s %s\n", s0, s1, s2); */
               
               SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "cluttercut", -SCIPinfinity(scip), cut_rhs, FALSE, FALSE, TRUE) );
               SCIP_CALL( SCIPaddVarsToRow(scip, cut, n_cutvars, cutvars, cutvals) );

               /* SCIPdebugMessage(" -> Clutter-cut <cluttercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n", */
               /*    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut), */
               /*    SCIPgetCutEfficacy(scip, NULL, cut), */
               /*    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut), */
               /*    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut)); */
               /* SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));   */
               

               SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
               if( *cutoff )
               {
                  SCIPdebugMessage("Clutter cut led to cutoff\n");
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  goto TERMINATE;
               }
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               (*nGen)++;
               n_cutsi0++;
               if( n_cutsi0 > maxcutsroundsingleton && maxcutsroundsingleton != -1 )
                  goto TERMINATE;
               if( *nGen > maxcutsround && maxcutsround != -1 )
                  goto TERMINATE;
            }

            n_cutvars2 = n_cutvars;
            cut_rhs2 = cut_rhs;

            if( maxcluttersize < 4 || mincluttersize > 4 )
            {
               cut_rhs = cut_rhs1;
               n_cutvars = n_cutvars1;
               continue;
            }

            for(i3 = i2+1; i3 < nposvars; i3++)
            {

               size3 = sizes[i3];
               if( size3 > maxsizeinclutter )
                  continue;


               ci3 = cis[i3];
               /* only clutters */
               if( (size0 < size3 && is_subset[size0][ci0][size3][ci3]) 
                  || (size1 < size3 && is_subset[size1][ci1][size3][ci3])
                  || (size2 < size3 && is_subset[size2][ci2][size3][ci3])
                  )
                  continue;
               
               lhs3 = update_cut(size3, ci3, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);
               
               init_union(size3, subsets[size3][ci3], union3);
               
               get_union(union0,union3,union03);
               get_union(union1,union3,union13);
               get_union(union2,union3,union23);
               get_union(union01,union3,union013);
               get_union(union02,union3,union023);
               get_union(union12,union3,union123);
               get_union(union012,union3,union0123);

               size03 = union03[0];
               size13 = union13[0];
               size23 = union23[0];
               size013 = union013[0];
               size023 = union023[0];
               size123 = union123[0];
               size0123 = union0123[0];

               ci03 = get_idx(tree, size03, union03+1);
               ci13 = get_idx(tree, size13, union13+1);
               ci23 = get_idx(tree, size23, union23+1);
               ci013 = get_idx(tree, size013, union013+1);
               ci023 = get_idx(tree, size023, union023+1);
               ci123 = get_idx(tree, size123, union123+1);
               ci0123 = get_idx(tree, size0123, union0123+1);
               
               lhs3 += update_cut(size03, ci03, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, 1);
               lhs3 += update_cut(size13, ci13, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, 1);
               lhs3 += update_cut(size23, ci23, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, 1);
               lhs3 += update_cut(size013, ci013, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);
               lhs3 += update_cut(size023, ci023, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);
               lhs3 += update_cut(size123, ci123, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, -1);
               lhs3 += update_cut(size0123, ci0123, lpvals, imsetvars, cutvars, cutvals, &n_cutvars, &cut_rhs, 1);
            
               if( lhs1 + lhs2 + lhs3 > cut_rhs )
               {
                  /* got a cut */
                  /* getName(s0, size0, ci0, imsetvars); */
                  /* getName(s1, size1, ci1, imsetvars); */
                  /* getName(s2, size2, ci2, imsetvars); */
                  /* printf("clutter: %s %s %s\n", s0, s1, s2); */
               
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "cluttercut", -SCIPinfinity(scip), cut_rhs, FALSE, FALSE, TRUE) );
                  SCIP_CALL( SCIPaddVarsToRow(scip, cut, n_cutvars, cutvars, cutvals) );

                  /* SCIPdebugMessage(" -> Clutter-cut <cluttercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n", */
                  /*    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut), */
                  /*    SCIPgetCutEfficacy(scip, NULL, cut), */
                  /*    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut), */
                  /*    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut)); */
                  /* SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));   */
                  

                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
                  if( *cutoff )
                  {
                     SCIPdebugMessage("Clutter cut led to cutoff\n");
                     SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                     goto TERMINATE;
                  }
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  (*nGen)++;
                  n_cutsi0++;
                  if( n_cutsi0 > maxcutsroundsingleton && maxcutsroundsingleton != -1 )
                     goto TERMINATE;
                  if( *nGen > maxcutsround && maxcutsround != -1 )
                     goto TERMINATE;

               }
               cut_rhs = cut_rhs2;
               n_cutvars = n_cutvars2;
            }
            cut_rhs = cut_rhs1;
            n_cutvars = n_cutvars1;
         }
         cut_rhs = 0;
         n_cutvars = 0;
      }
   }
TERMINATE:
   for( size = 1; size <= n; size++)
      SCIPfreeMemoryArray(scip, &(lpvals[size]));

   SCIPfreeMemoryArray(scip, &lpvals);
   SCIPfreeMemoryArray(scip, &sizes);
   SCIPfreeMemoryArray(scip, &cis);

   SCIPfreeMemoryArray(scip, &union1);
   SCIPfreeMemoryArray(scip, &union2);
   SCIPfreeMemoryArray(scip, &union3);
   SCIPfreeMemoryArray(scip, &union01);
   SCIPfreeMemoryArray(scip, &union02);
   SCIPfreeMemoryArray(scip, &union03);
   SCIPfreeMemoryArray(scip, &union12);
   SCIPfreeMemoryArray(scip, &union13);
   SCIPfreeMemoryArray(scip, &union23);
   SCIPfreeMemoryArray(scip, &union012);
   SCIPfreeMemoryArray(scip, &union013);
   SCIPfreeMemoryArray(scip, &union023);
   SCIPfreeMemoryArray(scip, &union123);
   SCIPfreeMemoryArray(scip, &union0123);
   return SCIP_OKAY;
}




/** Checks whether an acyclic digraph has any immoralities 
 *  @param scip The SCIP instance to which the constraint belongs.
 *  @param psd The parent set data on which the constraint is based.
 *  @param sol The acyclic digraph to check
 *  @return TRUE if digraph has no immoralities else FALSE
*/
static
SCIP_Bool is_feasible(
   SCIP* scip,
   ParentSetData* psd,
   SCIP_SOL* sol
   )
{
   int n = psd->n;
   int i;
   int k;
   int l;
   int npa;
   int l2;
   int pa1;
   int pa2;
   SCIP_VAR* edgevar;

   for( i = 0; i < n; ++i )
   {
      for( k = 0; k < psd->nParentSets[i]; k++ )
         if( SCIPisGT(scip, SCIPgetSolVal(scip, sol, psd->PaVars[i][k]), 0.5) )
         {
            npa = psd->nParents[i][k];
            if( npa > 1 )
            {
               /* any immoralities? */
               for( l = 0; l < npa - 1; l++ )
               {
                  pa1 = psd->ParentSets[i][k][l];
                  for( l2 = l+1; l2 < npa; l2++ )
                  {
                     pa2 = psd->ParentSets[i][k][l2];
                     edgevar = get_edge(psd,pa1,pa2);
                     if( edgevar == NULL || 
                        !SCIPisGT(scip, SCIPgetSolVal(scip, sol, edgevar), 0.5) )
                        return FALSE;
                  }
               }
            }
            /* no need to look for another selected parent set for variable i */
            break;
         }
   }

   return TRUE;
}


/** create imset variables from parent set ('family') variables 
 *  @return SCIP_OKAY if successful, or an appropriate error otherwise.
 */
static
SCIP_RETCODE makeimsetvars(
   SCIP* scip,               /**< The SCIP instance to which the constraint belongs */
   SCIP_CONSDATA** consdata  /**< Constraint data to be updated with imset variables */
   )
{
   ParentSetData* psd;
   int n;
   int nsizes;
   int i;
   int j;
   int k;
   int** powerset;
   int* set;
   int size;
   int size2;
   
   int ci;
   int cj;
   int* n_imsetvars      = NULL;               /* n_imsetvars[size] is the number of imsetvars for components of size size */
   int*** subsets        = NULL;               /* subsets[size][ci] is the subset associated with the cith imsetvar of size size */
   SCIP_VAR*** imsetvars = NULL;               /* imsetvars[size][ci] is the cith imsetvar of size size */
   SCIP_VAR**** pavars   = NULL;               /* pavars[size][ci] is the set of pavars whose sum = imsetvars[size][ci] */
   int** n_pavars        = NULL;               /* n_pavars[size][ci] is the number of pavars whose sum = imsetvars[size][ci] */
   TREE* tree            = NULL;               /* reverse of subsets: retrieves ci for a subset of size size (or -1 if there is none */

   SCIP_Bool found;

   char s[SCIP_MAXSTRLEN];
   char tmp[SCIP_MAXSTRLEN];

   SCIP_CONS* cons;


   SCIP_Bool**** is_subset;

   int* small;
   int* big;
   SCIP_Bool isa_subset;
   int errval;


   psd = (*consdata)->psd;
   n = (*consdata)->psd->n;
   nsizes = n+1;

   tree = (TREE *) malloc(sizeof(TREE));
   tree->idx = -1;
   tree->nchildren = 0;
   tree->children = NULL;

   /* allocate enough space for all sizes from 1 to n (index 0 is not used) */
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_imsetvars, nsizes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &subsets, nsizes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &imsetvars, nsizes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &pavars, nsizes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_pavars, nsizes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &is_subset, nsizes) );

   /* ensure all pointers are initially NULL (should be anyway)
      and initialise number of imsetvars for all sizes to be 0
   */

   for( size = 0; size < nsizes; size++ )
   {
      n_imsetvars[size] = 0;
      subsets[size] = NULL;
      imsetvars[size] = NULL;
      pavars[size] = NULL;
      n_pavars[size] = NULL;
   }

   for( i = 0; i < n; i++ )
   {
      for( k = 0; k < psd->nParentSets[i]; k++ )
      {
         SCIP_CALL( allsubsets(scip, psd->ParentSets[i][k], 
               psd->nParents[i][k], &powerset) );

         /* for each subset of the parents for i in kth parent set 
            there is some imsetvar */
         for( j = 1; j <= powerset[0][0]; j++ )
         {
            /* once size of set has been saved 
               can overwrite set[0] with i */
            set = powerset[j];
            size = set[0]+1;
            set[0] = i;
            SCIPsortInt(set, size);
            
            /* do we have an imsetvar for set? */
            ci = get_idx(tree, size, set);
            if( ci == -1  )
            {
               /* no imsetvar for this subset yet,
                  so assign appropriate index */
               ci = n_imsetvars[size];

               /* and store it */
               errval = put_index(tree, size, set, ci);
               if( errval == -1 )
               {
                  SCIPerrorMessage("Could not insert index in tree of indices.\n");
                  return SCIP_ERROR;
               }


               /* create name for new imset variable 
                  and copy subset */
               if( ci == 0 )
               {
                  /* first imsetvar of this size */
                  SCIP_CALL( SCIPallocMemoryArray(scip, &(subsets[size]), 1) );
                  SCIP_CALL( SCIPallocMemoryArray(scip, &(imsetvars[size]), 1) );
                  if( size > 2 )
                  {
                     SCIP_CALL( SCIPallocMemoryArray(scip, &(pavars[size]), 1) );
                     SCIP_CALL( SCIPallocMemoryArray(scip, &(n_pavars[size]), 1) );
                  }


               }
               else
               {
                  /* one more of this size */
                  SCIP_CALL( SCIPreallocMemoryArray(scip, &(subsets[size]), ci+1) );                  
                  SCIP_CALL( SCIPreallocMemoryArray(scip, &(imsetvars[size]), ci+1) );
                  if( size > 2 )
                  {
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &(pavars[size]), ci+1) );
                     SCIP_CALL( SCIPreallocMemoryArray(scip, &(n_pavars[size]), ci+1) );
                  }
               }
                  
               /* a new imsetvar (of size > 2 )initially has no family variables associated with it */
               if( size > 2 )
                  n_pavars[size][ci] = 0;
               
               /* make the name */
               SCIP_CALL( SCIPallocMemoryArray(scip, &(subsets[size][ci]), size) );
               SCIPsnprintf(s, SCIP_MAXSTRLEN, "c(");
               for( cj = 0; cj < size; cj++ ) 
               {
                  subsets[size][ci][cj] = set[cj];
                  
                  SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "#%s", psd->nodeNames[set[cj]]);
                  strcat(s, tmp);
               }
               strcat(s, ")");
               

               if( size == 2 )
               {
                  /*c-imset components of size 2 already exist as edge variables */
                  imsetvars[size][ci] = get_edge(psd,set[0],set[1]);
                  assert(imsetvars[size][ci] != NULL);
               }
               else
               {
                  /* imset variables for subsets of size 1 are fixed at 1 */

                  SCIP_CALL(SCIPcreateVarBasic(scip, &(imsetvars[size][ci]), s, ((size == 1) ? 1.0 : 0.0), 1.0,
                        0.0,
                        SCIP_VARTYPE_BINARY));
                  SCIP_CALL( SCIPaddVar(scip, imsetvars[size][ci]) );
               }
               /* SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, imsetvars[size][ci]) ); */
               /* SCIP_CALL( SCIPchgVarBranchPriority(scip, imsetvars[size][ci], imsetvarpriority) ); */
               n_imsetvars[size]++;
            }
            
            /* only build up list of associated family vars if size above 2 */
            if( size > 2 )
            {
               assert(n_pavars[size] != NULL );

               if( n_pavars[size][ci] == 0 )
                  SCIP_CALL( SCIPallocMemoryArray(scip, &(pavars[size][ci]), 1) );
               else
                  SCIP_CALL( SCIPreallocMemoryArray(scip, &(pavars[size][ci]), n_pavars[size][ci]+1) );
               
               pavars[size][ci][n_pavars[size][ci]++] = psd->PaVars[i][k];
            }
         }
         
         SCIP_CALL( free_allsubsets(scip, &powerset) );
      }
   }

   /* post all (necessary) linking constraints */
   for( size = 3; size < nsizes; size++ )
   {
      for( ci = 0; ci < n_imsetvars[size]; ci++ )
      {
         assert(imsetvars[size][ci] != NULL);

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, "link", 0, NULL, NULL, 0, 0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, imsetvars[size][ci], -1) );
         for( i = 0; i < n_pavars[size][ci]; i++)
         {
            assert(pavars[size][ci][i] != NULL);
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, pavars[size][ci][i], 1) );
         }
         SCIP_CALL( SCIPaddCons(scip, cons) );
         /* SCIP_CALL( SCIPprintCons(scip, cons, NULL) );  */
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      }
   }

   /* 
      store subset relations 
      consider pairs of subsets with indices (size,ci) and (size2,cj) 
      where size2 is bigger than size
   */

   for( size = 1; size <= n; size++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(is_subset[size]), n_imsetvars[size]) );
      for( ci = 0; ci < n_imsetvars[size]; ci++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &(is_subset[size][ci]), n+1) );
         small = subsets[size][ci];
         
         for( size2 = size+1; size2 <= n; size2++ )
         {
            SCIP_CALL( SCIPallocMemoryArray(scip, &(is_subset[size][ci][size2]), n_imsetvars[size2]) );
            for( cj = 0; cj < n_imsetvars[size2]; cj++ )
            {
               big = subsets[size2][cj];

               /* is small a subset of big? */
               isa_subset = TRUE;
              for( i = 0; i < size; i++ )
               {
                  found = FALSE;
                  for( j = 0; j < size2; j++ )
                  {
                     if( small[i] == big[j] )
                     {
                        found = TRUE;
                        break;
                     }
                  }
                  if( !found )
                  {
                     isa_subset = FALSE;
                     break;
                  }
               }
               is_subset[size][ci][size2][cj] = isa_subset;
            }
         }
      }
   }                  
   
   (*consdata)->n_imsetvars = n_imsetvars;
   (*consdata)->imsetvars = imsetvars;
   (*consdata)->subsets = subsets;
   (*consdata)->is_subset = is_subset;
   (*consdata)->tree = tree;
   (*consdata)->pavars = pavars;
   (*consdata)->n_pavars = n_pavars;
   

   return SCIP_OKAY;
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
   ParentSetData* psd
   )
{

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   /* no need to make a separate copy */
   /* SCIP_CALL( PS_copyParentSetData(scip, psd, &((*consdata)->psd)) ); */
   (*consdata)->psd = psd;
   (*consdata)->n_imsetvars = NULL;
   (*consdata)->imsetvars = NULL;
   (*consdata)->subsets = NULL;
   (*consdata)->is_subset = NULL;
   (*consdata)->tree = NULL;
   (*consdata)->pavars = NULL;
   (*consdata)->n_pavars = NULL;

   SCIP_CALL( makeimsetvars(scip, consdata) );

   return SCIP_OKAY;
}

/*
 * Linear constraint upgrading
 */

#ifdef LINCONSUPGD_PRIORITY
/** tries to upgrade a linear constraint into a chordal constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdChordal)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);

   /* check, if linear constraint can be upgraded to chordal constraint */
   upgrade = FALSE;
   /* TODO: put the constraint's properties here, in terms of the statistics given by nposbin, nnegbin, ... */

   if( upgrade )
   {
      SCIPdebugMessage("upgrading constraint <%s> to chordal constraint\n", SCIPconsGetName(cons));

      /* create the bin Chordal constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsChordal(scip, upgdcons, SCIPconsGetName(cons), nvars, vars, vals, lhs, rhs,
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
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define conshdlrCopyChordal NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeChordal)
{  /*lint --e{715}*/

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
SCIP_DECL_CONSINIT(consInitChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitChordal NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitChordal NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreChordal NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreChordal NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolChordal NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolChordal NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteChordal)
{  /*lint --e{715}*/

   int size;
   int ci;
   int size2;
   int n;
   int nsizes;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( consdata != NULL);
   assert( *consdata != NULL);
   assert( (*consdata)->n_imsetvars != NULL );
   assert( (*consdata)->imsetvars != NULL );
   assert( (*consdata)->subsets != NULL );

   SCIPdebugMessage("deleting chordal constraint <%s>.\n", SCIPconsGetName(cons));

   /* don't delete psd, this is not a copy! */

   n = (*consdata)->psd->n;
   nsizes = n + 1;

   for( size = 1; size < nsizes; size++ )
   {
      for( ci = 0; ci < (*consdata)->n_imsetvars[size]; ci++ )
      {
         for( size2 = size+1; size2 < nsizes; size2++ )
         {
            SCIPfreeMemoryArray(scip, &((*consdata)->is_subset[size][ci][size2]) );
         }
         SCIPfreeMemoryArray(scip, &((*consdata)->is_subset[size][ci]) );
      }
      SCIPfreeMemoryArray(scip, &((*consdata)->is_subset[size]) );
   }
   SCIPfreeMemoryArray(scip, &((*consdata)->is_subset) );


   for( size = 1; size < nsizes; size++ )
   {
      for( ci = 0; ci < (*consdata)->n_imsetvars[size]; ci++ )
      {
         SCIPfreeMemoryArray(scip, &((*consdata)->subsets[size][ci]));
      }
      if( (*consdata)->n_imsetvars[size] > 0 )
      {
         for( ci = 0; ci < (*consdata)->n_imsetvars[size]; ci++ )
            SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->imsetvars[size][ci])) );
         
         SCIPfreeMemoryArray(scip, &((*consdata)->subsets[size]));
         SCIPfreeMemoryArray(scip, &((*consdata)->imsetvars[size]));
      }
   }

   for( size = 3; size < nsizes; size++ )
   {
      for( ci = 0; ci < (*consdata)->n_imsetvars[size]; ci++ )
      {
         SCIPfreeMemoryArray(scip, &((*consdata)->pavars[size][ci]));
      }
      if( (*consdata)->n_imsetvars[size] > 0 )
      {
         SCIPfreeMemoryArray(scip, &((*consdata)->pavars[size]));
         SCIPfreeMemoryArray(scip, &((*consdata)->n_pavars[size]));
      }
   }

   

   SCIPfreeMemoryArray(scip, &((*consdata)->subsets));
   SCIPfreeMemoryArray(scip, &((*consdata)->imsetvars));
   SCIPfreeMemoryArray(scip, &((*consdata)->n_imsetvars));
   SCIPfreeMemoryArray(scip, &((*consdata)->pavars));
   SCIPfreeMemoryArray(scip, &((*consdata)->n_pavars));

   free_children((*consdata)->tree);
   free((*consdata)->tree);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}



/** transforms constraint data into data belonging to the transformed problem */
#if 0
static
SCIP_DECL_CONSTRANS(consTransChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consTransChordal NULL
#endif


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpChordal)
{  /*lint --e{715}*/
   int c;
   int nGen = 0;

   SCIP_VAR* smallvar;
   SCIP_VAR* bigvar;

   SCIP_CONSHDLRDATA* conshdlrdata;
   
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      int sizelim;
      int size;
      int ci;
      int cj;

      assert( conss != NULL );
      assert( conss[c] != NULL );
      SCIPdebugMessage("adding initial rows for chordal constraint <%s>.\n", SCIPconsGetName(conss[c]));

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      sizelim = ( conshdlrdata->monolim < consdata->psd->n ? conshdlrdata->monolim : consdata->psd->n );
      
      for( size = 2; size < sizelim; size++)
      {
         for( ci = 0; ci < consdata->n_imsetvars[size]; ci++ )
         {
            for( cj = 0; cj < consdata->n_imsetvars[size+1]; cj++ )
            {
               if( consdata->is_subset[size][ci][size+1][cj] ) 
               {
                  char s[SCIP_MAXSTRLEN];
                  SCIP_ROW* row;

                  smallvar =  consdata->imsetvars[size][ci];
                  bigvar = consdata->imsetvars[size+1][cj];

                  (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "mono#%s#%s",  SCIPvarGetName(smallvar), SCIPvarGetName(bigvar));
                  SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, s,  -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
                  SCIP_CALL( SCIPcacheRowExtensions(scip, row) );
                  SCIP_CALL( SCIPaddVarToRow(scip, row, bigvar, 1.0) );
                  SCIP_CALL( SCIPaddVarToRow(scip, row, smallvar, -1.0) );
                  SCIP_CALL( SCIPflushRowExtensions(scip, row) );
#ifdef SCIP_DEBUG
                  SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL)) );
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
   }
   SCIPdebugMessage("added %d simple monotoncity constraints.\n", nGen);
   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpChordal)
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
      SCIP_Bool cutoff = FALSE;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("separating LP solution for chordal constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( ChordalSeparate(scip, conshdlr, consdata, NULL, &nGen, &cutoff) );
      if ( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }
   if (nGen > 0)
      *result = SCIP_SEPARATED;
   SCIPdebugMessage("separated %d cuts.\n", nGen);

   return SCIP_OKAY;

}



/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolChordal)
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
      SCIP_Bool cutoff = FALSE;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("separating solution for linear ordering constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      
      *result = SCIP_DIDNOTFIND;
      SCIP_CALL( ChordalSeparate(scip, conshdlr, consdata, sol, &nGen, &cutoff) );
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
SCIP_DECL_CONSENFOLP(consEnfolpChordal)
{    

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
      ParentSetData* psd;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("enforcing lp solution for  chordal graph constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );
      psd = consdata->psd;

      if( !is_feasible(scip, psd, NULL) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all chordal graph constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsChordal)
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
      SCIPdebugMessage("enforcing pseudo solution for chordal graph constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );
      psd = consdata->psd;

      if( !is_feasible(scip, psd, NULL) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all chordal graph constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckChordal)
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

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("checking chordal graph constraint <%s>.\n", SCIPconsGetName(cons));

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );

      if( !is_feasible(scip, consdata->psd, sol) )
      {
         *result = SCIP_INFEASIBLE;
         SCIPdebugMessage("solution violates chordal graph constraint <%s>.\n", SCIPconsGetName(cons));
         return SCIP_OKAY;
      }
   }
   SCIPdebugMessage("all chordal graph constraints are feasible.\n");
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropChordal)
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
      int n;
      ParentSetData* psd;

      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      int i;
      int j;
      int j2;

      SCIP_VAR* arrowvar;
      SCIP_VAR* arrowvar2;
      SCIP_VAR* edgevar;

      int child;
      int k;
      int l;
      int pa;
      SCIP_Bool foundone;      
      SCIP_Bool foundboth;

      cons = conss[c];
      assert( cons != NULL );
      SCIPdebugMessage("propagating chordal graph constraint <%s>.\n", SCIPconsGetName(cons));
      
      *result = SCIP_DIDNOTFIND;

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );
      assert( consdata->psd != NULL );

      psd = consdata->psd;
      n = psd->n;
               
      for( i = 0; i < n; ++i )
      {
         for( j = 0; j < n; ++j )
         {
            if( i == j )
               continue;

            arrowvar =  get_arrow(psd,i,j);

            if( arrowvar != NULL && SCIPvarGetLbLocal(arrowvar) > 0.5 )
            {
               /* propagations which obtain when there is an arrow from j to i */

               for( j2 = j+1; j2 < n; ++j2 )
               {
                  if( i == j2 )
                     continue;
               
                  arrowvar2 =  get_arrow(psd,i,j2);
                  edgevar = get_edge(psd,j,j2);

                  if( arrowvar2 != NULL && SCIPvarGetLbLocal(arrowvar2) > 0.5 )
                  {

                     /* j and j2 are parents of i so must be married */

                     SCIPdebugMessage("<%s> and <%s> are set to TRUE so need to prevent an immorality \n", SCIPvarGetName(arrowvar), SCIPvarGetName(arrowvar2));
               
                     if( edgevar == NULL )
                     {                  
                        SCIPdebugMessage(" edge variable not found -> node infeasible.\n");
                        *result = SCIP_CUTOFF;
                        return SCIP_OKAY;
                     }

                     SCIP_CALL( SCIPtightenVarLb(scip,edgevar,1,TRUE,&infeasible,&tightened) );
                  
                     if( infeasible )
                     {
                        SCIPdebugMessage(" -> node infeasible.\n");
                        *result = SCIP_CUTOFF;
                        return SCIP_OKAY;
                     }
                     
                     if( tightened )
                     {
                        SCIPdebugMessage("Setting: %s to TRUE (to prevent an immorality) \n", SCIPvarGetName(edgevar));
                        ++nGen;
                     }
                  }
                  else if ( edgevar == NULL || SCIPvarGetUbLocal(edgevar) < 0.5 )
                  {
                     /* no marriage so cannot have j2 as parent of i */
                  
                     SCIPdebugMessage("<%s> set to TRUE but no possible marriage so need to prevent an immorality \n", SCIPvarGetName(arrowvar));
                  
                     if( arrowvar2 != NULL )
                     {
                        SCIP_CALL( SCIPtightenVarUb(scip,arrowvar2,0,TRUE,&infeasible,&tightened) );
                        
                        if( infeasible )
                        {
                           SCIPdebugMessage(" -> node infeasible.\n");
                           *result = SCIP_CUTOFF;
                           return SCIP_OKAY;
                        }
                        
                        if( tightened )
                        {
                           SCIPdebugMessage("Setting: %s to FALSE (to prevent an immorality) \n", SCIPvarGetName(arrowvar2));
                           ++nGen;
                        }
                     }
                  }
               }
            }
            edgevar = get_edge(psd,i,j);
            if( edgevar == NULL || SCIPvarGetUbLocal(edgevar) < 0.5 )
            {
               /* propagations which obtain when there is no edge from j to i */
               for( child = 0; child < n; ++child )
               {
                  if( child == i || child == j )
                     continue;
                  
                  for( k = 0; k < psd->nParentSets[child]; k++ )
                  {
                     if( psd->nParents[child][k] < 2 )
                        continue;

                     /* don't waste time on family variables already set to 0 */
                     if( SCIPvarGetUbLocal(psd->PaVars[child][k]) < 0.5 )
                        continue;
                     
                     foundone = FALSE;
                     foundboth = FALSE;
                     for( l = 0; l < psd->nParents[child][k]; l++ )
                     {
                        pa = psd->ParentSets[child][k][l];
                        if( pa == i || pa == j )
                        {
                           if( foundone )
                           {
                              foundboth = TRUE;
                              break;
                           }
                           else
                           {
                              foundone = TRUE;
                           }
                        }
                     }

                     if( foundboth )
                     {
                        SCIPdebugMessage("Edge from %s to %s missing so need to rule out %s (to prevent an immorality) \n", 
                           psd->nodeNames[i], psd->nodeNames[j], SCIPvarGetName(psd->PaVars[child][k]));
                                                   
                        SCIP_CALL( SCIPtightenVarUb(scip,psd->PaVars[child][k],0,TRUE,&infeasible,&tightened) );
                        
                        if( infeasible )
                        {
                           SCIPdebugMessage(" -> node infeasible.\n");
                           *result = SCIP_CUTOFF;
                           return SCIP_OKAY;
                        }
                        
                        if( tightened )
                        {
                           SCIPdebugMessage("Setting: %s to FALSE (to prevent an immorality) \n", SCIPvarGetName(psd->PaVars[child][k]));
                           ++nGen;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   if (nGen > 0)
      *result = SCIP_REDUCEDDOM;
   SCIPdebugMessage("propagated %d domains.\n", nGen);
   
   return SCIP_OKAY;
}



/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolChordal)
{  /*lint --e{715}*/

   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;
   SCIP_Bool cutoff;
   SCIP_Bool fixed;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->vanilla )
      return SCIP_OKAY;
   
   
   for( c = 0; c < nconss && !SCIPisStopped(scip); ++c )
   {
      SCIP_CONS* cons;
      SCIP_CONSDATA* consdata;
      int k;

      assert(*result != SCIP_CUTOFF);

      cons = conss[c];
      assert(cons != NULL);
      consdata = SCIPconsGetData(cons);
      assert(consdata != NULL);

      for( k = 0; k < consdata->psd->nParentSets[0]; k++ )
      {
         if(  consdata->psd->nParents[0][k] > 0 )
         {
            SCIPdebugMessage("trying to fix <%s> to 0 so that variable 0 becomes a source.\n", SCIPvarGetName(consdata->psd->PaVars[0][k]));
            SCIP_CALL( SCIPfixVar(scip, consdata->psd->PaVars[0][k], 0, &cutoff, &fixed) ); 
            if( cutoff )
            {
               SCIPdebugMessage("chordal constraint <%s>: infeasible fixing <%s> == 0\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->psd->PaVars[0][k]) );
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            *result = SCIP_SUCCESS;
         }
      }
   }
   return SCIP_OKAY;
}



/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropChordal NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockChordal)
{ 
   SCIP_CONSDATA* consdata;
   ParentSetData* psd;
   
   int i;
   int j;
   SCIP_VAR* edgevar;


   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   SCIPdebugMessage("Locking chordal graph constraint <%s>.\n", SCIPconsGetName(cons));

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->psd != NULL );
   psd = consdata->psd;

   /* enough to down-lock the edge variables 
      since we just rule out immoralities (see CHECK method) */
   for( i = 0; i < psd->n; i++ )
      for( j = i+1; j < psd->n; j++ )
      {
         edgevar = get_edge(psd,i,j);
         if( edgevar != NULL)
            SCIP_CALL( SCIPaddVarLocks(scip, edgevar, nlockspos, nlocksneg) );
      }
   
   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveChordal NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveChordal NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableChordal NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableChordal NULL
#endif

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsChordal NULL
#endif


/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintChordal NULL
#endif


/** constraint copying method of constraint handler */
#if 0
static
SCIP_DECL_CONSCOPY(consCopyChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consCopyChordal NULL
#endif


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParseChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParseChordal NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsChordal NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsChordal NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsChordal)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of chordal constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsChordal NULL
#endif


/*
 * constraint specific interface methods
 */

/** creates the handler for chordal constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrChordal(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create chordal constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   /* TODO: (optional) create constraint handler specific data here */

   conshdlr = NULL;

   /* include constraint handler */

   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpChordal, consEnfopsChordal, consCheckChordal, consLockChordal,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   /* SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyChordal, consCopyChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveChordal) ); */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteChordal) ); 
   /* SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolChordal) ); */
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeChordal) );
   /* SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreChordal) ); */
   /* SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolChordal) ); */
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpChordal) ); 
   /* SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseChordal) ); */
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolChordal, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) ); 
   /* SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintChordal) ); */
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropChordal, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP, CONSHDLR_PROP_TIMING) );
   /* SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropChordal) ); */
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpChordal, consSepasolChordal, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) ); 
   /* SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransChordal) ); */

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/monolim",
         "limit on size of subset corresponding to the lower bound in a simple monotonicity constraint",
         &conshdlrdata->monolim, FALSE, DEFAULT_MONOLIM, 2, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/vanilla",
         "whether to assume 'vanilla' learning of chordal graphs",
         &conshdlrdata->vanilla, FALSE, DEFAULT_VANILLA, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/onesingleton",
         "whether clutters for cuts should contain only one singleton",
         &conshdlrdata->onesingleton, FALSE, DEFAULT_ONESINGLETON, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxcutsround",
         "maximum number of cuts generated in a separation round (-1 is no limit)",
         &conshdlrdata->maxcutsround, FALSE, DEFAULT_MAXCUTSROUND, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxcutsroundsingleton",
         "maximum number of cuts generated in a separation round for any singleton (-1 is no limit)",
         &conshdlrdata->maxcutsroundsingleton, FALSE, DEFAULT_MAXCUTSROUNDSINGLETON, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/mincluttersize",
         "minimal size of clutter for a cut",
         &conshdlrdata->mincluttersize, FALSE, DEFAULT_MINCLUTTERSIZE, 2, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxcluttersize",
         "maximal size of clutter for a cut (-1 is no limit)",
         &conshdlrdata->maxcluttersize, FALSE, DEFAULT_MAXCLUTTERSIZE, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
         "constraints/"CONSHDLR_NAME"/maxsizeinclutter",
         "maximal size of a subset in a clutter for a cut (-1 is no limit)",
         &conshdlrdata->maxsizeinclutter, FALSE, DEFAULT_MAXSIZEINCLUTTER, -1, INT_MAX, NULL, NULL));


   return SCIP_OKAY;
}

/** creates and captures a chordal constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */

SCIP_RETCODE SCIPcreateConsChordal(
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
   /* TODO: (optional) modify the definition of the SCIPcreateConsChordal() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the chordal constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("chordal constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* Initialise constraint data */
   SCIP_CALL( createConsData(scip, &consdata, psd) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a chordal constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicChordal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd
   )
{
   SCIP_CALL( SCIPcreateConsChordal(scip, cons, name, psd,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}
