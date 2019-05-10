/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2017 James Cussens                       *
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

/** @file
 * Functions for finding matroid cuts
 */

/*#define SCIP_DEBUG*/
#include "matroid_cuts.h"
#include "scip/scipdefplugins.h"
#include <string.h>

#define min(A,B) ((A) > (B) ? (B) : (A))

/** Data for subIP for searching for matroid cuts */
typedef struct
{
   SCIP*           subscip;          /**< sub MIP for finding good matroids for cutting planes */
   SCIP_VAR**      in_ground_set;    /**< vars indicating BN variables being in the matroid's ground set */
   SCIP_VAR**      full_graph;       /**< vars indicating BN variables which appear in cut for full graph */
   SCIP_VAR****    circuit_vars;     /**< circuit_vars[size][i][ci] is the var corresponding to the cith
                                        circuit of size 'size' containing BN variable i */
   int** n_circuits_i;     /**< n_circuits_i[size][i] is the number of circuits of size 'size' containing i */
   int**** circuits_i;     /**< circuits_i[size][i][ci] is the cith circuit of size 'size' containing i */

   SCIP_VAR***    family;  /** family variables which are positive in the current LP solution */
   
} MATROID_AUXIPDATA;


/** Frees sub-IP data */
static SCIP_RETCODE AuxIPDataFree(
   SCIP* scip,                      /**< (Main) SCIP instance */
   MATROID_AUXIPDATA* auxipdata,    /**< Data for the sub-IP (to be freed) */
   int n,                           /**< Number of BN variables in the acyclicity (dagcluster) constraint */
   int circuit_size_lim             /**< The limit on the size of circuits */
)
{
   int i;
   int size;
   int ci;
   
   if( auxipdata->subscip != NULL )
      SCIP_CALL( SCIPfree(&(auxipdata->subscip)) );

   SCIPfreeMemoryArray(scip, &(auxipdata->in_ground_set));
   SCIPfreeMemoryArray(scip, &(auxipdata->full_graph));

   for( size = 2; size <= circuit_size_lim; size++ )
   {
      for( i = 0; i < n ; ++i )
      {
         for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ++ci )
         {
            assert((auxipdata->circuits_i)[size][i][ci] != NULL );
            /* printf("%d\n",(auxipdata->circuits_i)[size][i][ci][0]); */
            if( (auxipdata->circuits_i)[size][i][ci][0] == i )
               /* free((auxipdata->circuits_i)[size][i][ci]); */
               SCIPfreeMemoryArray(scip, &((auxipdata->circuits_i)[size][i][ci]));
         }
      }
   }

   for( size = 2; size <= circuit_size_lim; size++ )
   {
      for( i = 0 ; i < n ; ++i )
      {
         SCIPfreeMemoryArray(scip, &((auxipdata->circuits_i)[size][i]));
         SCIPfreeMemoryArray(scip, &((auxipdata->circuit_vars)[size][i]));
      }

      SCIPfreeMemoryArray(scip, &((auxipdata->n_circuits_i)[size]));
      SCIPfreeMemoryArray(scip, &((auxipdata->circuits_i)[size]));
      SCIPfreeMemoryArray(scip, &((auxipdata->circuit_vars)[size]));
   }
   SCIPfreeMemoryArray(scip, &(auxipdata->n_circuits_i));
   SCIPfreeMemoryArray(scip, &(auxipdata->circuits_i));
   SCIPfreeMemoryArray(scip, &(auxipdata->circuit_vars));

   for( i = 0 ; i < n ; ++i )
      SCIPfreeMemoryArray(scip, &((auxipdata->family)[i]));
   SCIPfreeMemoryArray(scip, &(auxipdata->family));
   
   SCIPfreeMemory(scip, &auxipdata);
   auxipdata = NULL;

   return SCIP_OKAY;
}

/** computes 2^i */
static int mypow2(
   int i           /**< Function returns 2^i */
   )
{

   int j;
   int res = 1;

   assert(i > -1);

   for( j = 0; j < i; j++)
      res *= 2;

   return res;
}
         
      
/** is \f$a\f$ a subset of \f$b \cup \{i\}\f$ ? */
static SCIP_Bool subseti(
   int i,        /**< an element */
   int n_a,      /**< size of set a */
   int* a,       /**< set a */
   int n_b,      /**< size of set b */
   int* b        /**< set b */
   )
{
   int la;
   int lb;
   SCIP_Bool found;
   int elt;
   
   for( la = 0; la < n_a; la++ )
   {
      elt = a[la];
      
      if( elt == i )
         continue;
      
      found = FALSE;
      for( lb = 0; lb < n_b; lb++ )
         if( elt == b[lb] )
         {
            found = TRUE;
            break;
         }
      if( !found )
         return FALSE;
   }
   /* for( la = 0; la < n_a; la++ ) */
   /*    printf("%d,",a[la]); */
   /* printf("subset of "); */
   /* printf("%d,",i); */
   /* for( lb = 0; lb < n_b; lb++ ) */
   /*    printf("%d,",b[lb]); */
   /* printf("\n"); */
   return TRUE;
}

/** is 'elt' a member of set 'a'? */
static SCIP_Bool element_of(
   int n_a,   /**< size of set a */
   int* a,    /**< set a */
   int elt    /**< an element */
   )
{
   int ia;
   for( ia = 0; ia < n_a; ia++)
   {
      if( a[ia] == elt )
         return TRUE;

      if( a[ia] > elt )
         return FALSE;
   }
   return FALSE;
}
         
/** is 'a' a subset_of of 'b' ? */
static SCIP_Bool subset_of(
   int n_a,  /**< size of set a */
   int* a,   /**< set a */
   int n_b,  /**< size of set b */
   int* b    /**< set b */
   )
{
   int ia;
   int ib = 0;

   assert( a != NULL );
   assert( b != NULL );
   assert( n_a > 1);
   assert( n_b > 1);
   
   if( n_a > n_b)
      return FALSE;

   /* sets are represented as ordered arrays of ints */
   for( ia = 0; ia < n_a; ia++ )
   {
      assert( ib < n_b );

      while( b[ib] < a[ia] )
      {
         if( ib < n_b)
            ib++;
         else
            return FALSE;
      }
      if( b[ib] == a[ia] )
      {
         /* have found ia+1 elements of a in b
            so n_a-(ia+1) remain to be found
            there are n_b-(ib+1) elements left in b
            Need n_a-(ia+1) =< n_b-(ib+1)
         */
         if( n_a - ia > n_b - ib )
            return FALSE;
         else
            ib++;
      }
      else
         return FALSE;
   }
   return TRUE;
}

/** both = \f$a \cup b \setminus \{i\}\f$.
    Both a and b are assumed to contain i.
    Sufficient space must have already been allocated for 'both'
    before calling this function
*/
static int myunioni(
   int i,    /** an element (of both a and b) */
   int n_a,  /** size of set a */
   int* a,   /** set a */
   int n_b,  /** size of set b */
   int* b,   /** set b  */
   int* both /** set both (output) */ 
   )
{
   int ia=0;
   int ib=0;
   int iboth=0;
   int j;
   
   while(ia < n_a && ib < n_b )
   {

      /* printf("%d,%d,%d,%d,%d\n",ia,ib,a[ia],b[ib],iboth); */
      
      if( a[ia] < b[ib] )
      {
         both[iboth++] = a[ia++];
         /* printf("A\n"); */
      }
      else if( b[ib] < a[ia] )
      {
         both[iboth++] = b[ib++];
         /* printf("B\n"); */
      }
      else
      {
         /* printf("C\n"); */
         if( a[ia] != i )
            both[iboth++] = a[ia];
         ia++;
         ib++;
      }
      /* printf("Now:ia=%d,ib=%d,iboth=%d\n",ia,ib,iboth); */
   }
   for( j = ia; j < n_a; j++ )
      both[iboth++] = a[j];
   for( j = ib; j < n_b; j++ )
      both[iboth++] = b[j];

   /* printf("Done\n"); */
   return iboth;
}

/** Computes n-choose-k */   
static int choose(
   int n,  /**< number of elements to choose from */
   int k   /**< size of chosen sets */
   )
{
   if( k == 0)
      return 1;

   return (n * choose(n-1,k-1)) / k;
}

/** Main function for finding matroid cuts */
extern SCIP_RETCODE Matroid_findCuts(
   SCIP*           scip,                    /**< SCIP data structure */
   ParentSetData*  psd,                     /**< family variable information */
   SolutionInfo*   solinfo,                 /**< information about the solution to be separated */
   SCIP_SOL*       sol,                     /**< solution to be separated */
   int*            nGen,                    /**< *nGen is number of cutting planes added ( even non-efficacious ones are added ) */
   SCIP_CONSHDLR*  conshdlr,                /**< constraint handler */
   SCIP_Bool       addtopool,               /**< whether to add any found cut to the global cut pool */
   SCIP_Bool       forcecuts,               /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,   /**< to return whether an efficacious cutting plane was found */
   SCIP_Real       limits_time,             /**< limit on how long to spend sub-IP solving */
   SCIP_Real       limits_gap,              /**< limit on size of gap in sub-IP */
   SCIP_Real       limits_absgap,           /**< limit on size of the absolute gap in sub-IP */
   int             circuit_size_lim,        /**< upper bound on size of circuits */
   int             ground_set_size_lim,        /**< upper bound on size of ground set */
   int*            must_be_included,        /**< set of nodes which must be included in any found matroid */
   int             n_must_be_included,      /**< size of the set of nodes which must be included in any matroid */
   int*            must_be_excluded,        /**< set of nodes which must be excluded from any matroid */
   int             n_must_be_excluded,      /**< size of the set of nodes which must be excluded from any found matroid */
   SCIP_Bool*      cutoff                   /**< cutoff = TRUE if a cut is added which leads to a cutoff ( set by SCIPaddRow ) */
)
{

   int* n_circuitsx;       /* n_circuitsx[size] is the number of circuits of size 'size' */
   int*** circuitsx;       /* circuitsx[size][ci] is the cith circuit of size 'size' */
   
   MATROID_AUXIPDATA* auxipdata;          /* data for subscip */
   int size;
   int i;
   int n_choices;
   char s[SCIP_MAXSTRLEN];
   int j;
   int* this_circuit;   /* not called 'circuit' to avoid name clash with function in circuit_cuts! */
   SCIP_VAR* var;
   SCIP_VAR* twovars[2];
   SCIP_CONS* cons;
   int previous_size;
   int ci;
   int* small_circuit;
   int cj;
   int elt;
   char tmp[SCIP_MAXSTRLEN];
   int** small_circuits;
   int bigger_size;
   int** big_circuits;
   int big_ci;
   int* big_circuit;
   int ki;
   int k;
   SCIP_Real val;
   int n_pas;
   int* pas;
   SCIP_VAR** vars;
   int nvars;
   int maxnvars;
   int nvars_max;
   int** circuits;
   int size2;
   int start;
   int* circuit2;
   int* both;
   int n_both;
   int ci2;
   int size3;
   int* circuit3;
   int ci3;

   SCIP_STATUS status;
   int si;
   int nsols;
   SCIP_SOL** subscip_sols;
   SCIP_SOL* subscip_sol;

   SCIP_Longint* weights;

   int**** set_circuits;
   int** n_set_circuits;

   int rhs;
   SCIP_ROW* cut;
   SCIP_Bool included;
   
   /* /\* debugging *\/ */
   /* int a[3] = {0,1,3}; */
   /* int b[3] = {0,2,3}; */
   /* int c[3] = {-1,-1,-1}; */

   /* n_both = myunioni(0,3,a,3,b,c); */

   SCIP_CALL( SCIPallocMemory(scip, &auxipdata) );
   auxipdata->subscip = NULL;
   auxipdata->in_ground_set = NULL;
   auxipdata->full_graph = NULL;
   auxipdata->circuit_vars = NULL;
   auxipdata->n_circuits_i = NULL;
   auxipdata->circuits_i = NULL;
   auxipdata->family = NULL;

   /* allocate memory for circuits and variables */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->in_ground_set), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->full_graph), psd->n) );
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_circuitsx, circuit_size_lim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &circuitsx, circuit_size_lim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->n_circuits_i), circuit_size_lim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->circuits_i), circuit_size_lim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->circuit_vars), circuit_size_lim) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &set_circuits, circuit_size_lim) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &n_set_circuits, circuit_size_lim) );

   /* don't consider circuits of size 0 or 1
      so set all to NULL / -1 to indicate that they should not
      be even looked at
   */
   
   for( size = 0; size < 2; size++ )
   {
      (auxipdata->circuit_vars)[size] = NULL;
      n_circuitsx[size] = -1;
      circuitsx[size] = NULL;
      (auxipdata->n_circuits_i)[size] = NULL;
      (auxipdata->circuits_i)[size] = NULL;
   }
      
   for( size = 2; size <= circuit_size_lim; size++ )
   {
      n_circuitsx[size] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->n_circuits_i[size]), psd->n) );
      for( i = 0; i < psd->n; ++i )
         (auxipdata->n_circuits_i)[size][i] = 0;

      SCIP_CALL( SCIPallocMemoryArray(scip, &(set_circuits[size]), psd->n) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(n_set_circuits[size]), psd->n) );
      
      SCIP_CALL( SCIPallocMemoryArray(scip, &(circuitsx[size]), choose(psd->n,size)) );
      /* printf("%d subsets from %d of size %d\n", choose(psd->n,size), psd->n, size); */
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->circuits_i[size]), psd->n));
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->circuit_vars[size]), psd->n));
      for( i = 0; i < psd->n; ++i )
      {
         n_choices = choose((psd->n)-1,size-1);
         /* printf("%d subsets from %d of size %d containing %d\n", n_choices, psd->n, size, i); */
         SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->circuits_i[size][i]), n_choices) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->circuit_vars[size][i]), n_choices) );
         SCIP_CALL( SCIPallocMemoryArray(scip, &set_circuits[size][i], n_choices) );
      }
   }

   /* set up sub-MIP */
   SCIP_CALL( SCIPcreate(&(auxipdata->subscip)) );
   SCIP_CALL( SCIPincludeDefaultPlugins(auxipdata->subscip) );
   SCIP_CALL( SCIPcreateProb(auxipdata->subscip, "Matroid separating MIP", NULL, NULL , NULL , NULL , NULL , NULL , NULL) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/time", limits_time) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/gap", limits_gap) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/absgap", limits_absgap) );

   /* SubMIP variables indicating being in the ground set of the matroid 
      and indicating that full graph parent is in the cut
   */

   SCIP_CALL( SCIPallocMemoryArray(scip, &weights, psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "I_g#%s", psd->nodeNames[i]);
      SCIP_CALL( SCIPcreateVarBasic(auxipdata->subscip, &(auxipdata->in_ground_set[i]), s, 0.0, 1.0, 0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->in_ground_set[i]) );
      /* SCIP_CALL( SCIPreleaseVar(auxipdata->subscip, &(auxipdata->in_ground_set[i])) ); */
      
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "I_f#%s", psd->nodeNames[i]);
      SCIP_CALL( SCIPcreateVarBasic(auxipdata->subscip, &(auxipdata->full_graph[i]), s, 0.0, 1.0, -1, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->full_graph[i]) );
      /* SCIP_CALL( SCIPreleaseVar(auxipdata->subscip, &(auxipdata->full_graph[i])) ); */

      weights[i] = 1;
   }

   SCIP_CALL( SCIPcreateConsBasicKnapsack(
         auxipdata->subscip, &cons, "ground_set_size_lim", psd->n, auxipdata->in_ground_set, weights, (SCIP_Longint) ground_set_size_lim) );
   SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );
   /*SCIP_CALL( SCIPprintCons(scip, cons, NULL)  ); */
   SCIP_CALL( SCIPreleaseCons(auxipdata->subscip, &cons) );

   SCIPfreeMemoryArray(scip, &weights);


/* create circuits and circuit vars and add constraints on I_f and I_g variables along the way */
   
   /* create circuits of size 2 */
   for( i = 0; i < psd->n; i++)
      for( j = i+1; j < psd->n; j++)
      {

         /* printf("2-cluster %d %d\n",i,j); */
         
         /* create circuit */
         this_circuit = NULL;
         SCIP_CALL( SCIPallocMemoryArray(scip, &this_circuit, 2) );
         /* this_circuit = (int *) malloc(2 * sizeof(int)); */
         assert( this_circuit != NULL );
         this_circuit[0] = i;
         this_circuit[1] = j;
         
         /* create corresponding variable */
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "c#%s#%s", psd->nodeNames[i], psd->nodeNames[j]);
         SCIP_CALL( SCIPcreateVarBasic(auxipdata->subscip, &var, s, 0.0, 1.0, 0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(auxipdata->subscip, var) );

         SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, var, &(twovars[0])) );

         /* c(i,j) =< I_g(i) */
         twovars[1] = auxipdata->in_ground_set[i];
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "g(%s,%s)", SCIPvarGetName(var),SCIPvarGetName(twovars[1]));
         SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, 2, twovars) );
         /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );

         /* c(i,j) =< I_g(j) */
         twovars[1] = auxipdata->in_ground_set[j];
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "g(%s,%s)", SCIPvarGetName(var),SCIPvarGetName(twovars[1]));
         SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, 2, twovars) );
         /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );

         /* c(i,j) =< I_f(j) [ NB no constraint on I_f(i) ] */
         twovars[1] = auxipdata->full_graph[j];
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "f(%s,%s)", SCIPvarGetName(var),SCIPvarGetName(twovars[1]));
         SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, 2, twovars) );
         /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );

         /* and store */
         circuitsx[2][n_circuitsx[2]++] = this_circuit;
         (auxipdata->circuits_i)[2][i][(auxipdata->n_circuits_i)[2][i]] = this_circuit;
         (auxipdata->circuits_i)[2][j][(auxipdata->n_circuits_i)[2][j]] = this_circuit;
         (auxipdata->circuit_vars)[2][i][(auxipdata->n_circuits_i)[2][i]++] = var;
         (auxipdata->circuit_vars)[2][j][(auxipdata->n_circuits_i)[2][j]++] = var;

         assert( n_circuitsx[2] <= choose(psd->n,2) );
         /* printf("%d:%d\n",(auxipdata->n_circuits_i)[2][i], choose((psd->n)-1,1) ); */
         assert( (auxipdata->n_circuits_i)[2][i] <= choose((psd->n)-1,1) );
         assert( (auxipdata->n_circuits_i)[2][j] <= choose((psd->n)-1,1) );
         
      }
   
   /* create circuits of size greater than 2 
      by adding elements which are greater than the last element of an existing circuit
   */
   for( size = 3; size <= circuit_size_lim; size++ )
   {
      /* printf("size = %d\n", size); */
      previous_size = size-1;
      for( ci = 0; ci < n_circuitsx[previous_size]; ci++ )
      {
         small_circuit = circuitsx[previous_size][ci];
         assert( small_circuit != NULL );
         
         for( i = small_circuit[previous_size-1]+1; i < psd->n; i++)
         {
            assert( i > -1 && i < psd->n );

            this_circuit = NULL;
            SCIP_CALL( SCIPallocMemoryArray(scip, &this_circuit, size) ); 
            /* this_circuit = (int *) malloc( size * sizeof(int));  */
            assert( this_circuit != NULL );
            
            /* copy elements from existing circuit and build up name for variable  */
            SCIPsnprintf(s, SCIP_MAXSTRLEN, "c");
            for( cj = 0; cj < previous_size; cj++ )
            {
               elt = small_circuit[cj];
               assert( elt > -1 && elt < psd->n);
               
               SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "#%s", psd->nodeNames[elt]);
               strcat(s, tmp);

               this_circuit[cj] = elt;
               (auxipdata->circuits_i)[size][elt][(auxipdata->n_circuits_i)[size][elt]] = this_circuit;

            }
            /* ... and then add new element */
            this_circuit[cj] = i;

            /* create variable */
            SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "#%s", psd->nodeNames[i]);
            strcat(s, tmp);
            SCIP_CALL( SCIPcreateVarBasic(auxipdata->subscip, &var, s, 0.0, 1.0, 0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(auxipdata->subscip, var) );
            /* SCIP_CALL( SCIPreleaseVar(auxipdata->subscip, &var) ); */

            /* add constraints */
            SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, var, &(twovars[0])) );
            for( cj = 0; cj < previous_size; cj++ )
            {

               elt = small_circuit[cj];
               (auxipdata->circuit_vars)[size][elt][(auxipdata->n_circuits_i)[size][elt]++] = var;
               assert( (auxipdata->n_circuits_i)[size][elt] <= choose((psd->n)-1,size-1) );
               
               /* this_circuit =< I_g(elt) */
               twovars[1] = auxipdata->in_ground_set[elt];

               SCIPsnprintf(s, SCIP_MAXSTRLEN, "g(%s,%s)", SCIPvarGetName(var),SCIPvarGetName(twovars[1]));
               SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, 2, twovars) );
               /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
               SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );

            }

            /* store */
            (auxipdata->circuits_i)[size][i][(auxipdata->n_circuits_i)[size][i]] = this_circuit;
            (auxipdata->circuit_vars)[size][i][(auxipdata->n_circuits_i)[size][i]++] = var;
            circuitsx[size][n_circuitsx[size]++] = this_circuit;

            assert( n_circuitsx[size] <= choose(psd->n,size) );
            assert( (auxipdata->n_circuits_i)[size][i] <= choose((psd->n)-1,size-1) );

            /* this_circuit =< I_g(i) */
            twovars[1] = auxipdata->in_ground_set[i];
            SCIPsnprintf(s, SCIP_MAXSTRLEN, "g(%s,%s)", SCIPvarGetName(var),SCIPvarGetName(twovars[1]));
            SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, 2, twovars) );
            /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
            SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );

            /* this_circuit =< I_f(i) */
            twovars[1] = auxipdata->full_graph[i];
            SCIPsnprintf(s, SCIP_MAXSTRLEN, "f(%s,%s)", SCIPvarGetName(var),SCIPvarGetName(twovars[1]));
            SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, 2, twovars) );
            /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
            SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );

         }
      }
   }

   /* finished with circuitsx, n_circuitsx */
   for( size = 2; size <= circuit_size_lim; size++ )
      SCIPfreeMemoryArray(scip,&(circuitsx[size]));
   SCIPfreeMemoryArray(scip,&circuitsx);
   SCIPfreeMemoryArray(scip,&n_circuitsx);

   
   /* /\* for debugging only *\/ */
   /* for( size = 2; size <= circuit_size_lim; size++ ) */
   /* { */
   /*    printf("Circuits of size %d ",size); */
   /*    for( i = 0; i < psd->n; i++) */
   /*    { */
   /*       printf("containing %d\n",i); */
   /*       circuits = (auxipdata->circuits_i)[size][i]; */
   /*       for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ci++ ) */
   /*       { */
   /*          for( j = 0; j < size; j++) */
   /*             printf("%d,",circuits[ci][j]); */
   /*          printf(" %s\n", SCIPvarGetName((auxipdata->circuit_vars)[size][i][ci])); */
   /*       } */
   /*       printf("\n"); */
   /*    } */
   /*    printf("\n"); */
   /* } */

   
   /* simple binary constraints that circuits form a clutter */

   for( size = 2; size <= circuit_size_lim; size++ )
   {
      for( i = 0; i < psd->n; i++)
      {
         small_circuits = (auxipdata->circuits_i)[size][i];
         for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ++ci )
         {
            small_circuit = small_circuits[ci];
            assert( small_circuit != NULL );
            
            /* this prevents posting the same constraint twice */
            if( small_circuit[0] != i )
               continue;

            twovars[0] = (auxipdata->circuit_vars)[size][i][ci];
            assert( twovars[0] != NULL );
            
            for( bigger_size = size+1; bigger_size <= circuit_size_lim; bigger_size++ )
            {
               big_circuits = (auxipdata->circuits_i)[bigger_size][i];
               for( big_ci = 0; big_ci < (auxipdata->n_circuits_i)[bigger_size][i]; ++big_ci )
               {
                  big_circuit = big_circuits[big_ci];
                  assert( big_circuit != NULL );

                  /* for(j = 0; j < size; j++) */
                  /*    printf("%d,",small_circuit[j]); */
                  /* printf("\n"); */

                  /* for(j = 0; j < bigger_size; j++) */
                  /*    printf("%d,",big_circuit[j]); */
                  /* printf("\n"); */

                  
                  if( subset_of(size,small_circuit,bigger_size,big_circuit) )
                  {
                     twovars[1] = (auxipdata->circuit_vars)[bigger_size][i][big_ci];
                     assert( twovars[1] != NULL );
                     SCIPsnprintf(s, SCIP_MAXSTRLEN, "i(%s,%s)", SCIPvarGetName(twovars[0]),SCIPvarGetName(twovars[1]));
                     SCIP_CALL( SCIPcreateConsBasicSetpack(auxipdata->subscip, &cons, s, 2, twovars) );
                     /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
                     SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );
                  }
               }
            }
         }
      }
   }

   /* the matroid must be connected: if two distinct elements are both in the ground
      set, then there is a circuit containing both */

   maxnvars = 2;
   for( size = 0; size <= circuit_size_lim-2; size++ )
      maxnvars += choose((psd->n)-2,size);
   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, maxnvars) );
   for( i = 0; i < psd->n; i++)
   {
      SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, (auxipdata->in_ground_set)[i], &(vars[0])) );
      for( j = i+1; j < psd->n; j++)
      {
         SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, (auxipdata->in_ground_set)[j], &(vars[1])) );
         nvars = 2;
         for( size = 2; size <= circuit_size_lim; size++ )
         {
            for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ci++ )
            {
               this_circuit = (auxipdata->circuits_i)[size][i][ci];
               if( element_of(size,this_circuit,j) )
               {
                  vars[nvars++] = (auxipdata->circuit_vars)[size][i][ci];
                  assert( nvars <= maxnvars );
               }
            }
         }
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "cnct#%d#%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, nvars, vars) );
         /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );
      }
   }
   
   /* circuit elimination constraints */
   maxnvars = mypow2(2*(circuit_size_lim-1))-2*(circuit_size_lim-1)-1+2;
   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, maxnvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &both, 2*circuit_size_lim - 2) ); 
   for( i = 0; i < psd->n; ++i )
   {
      for( size = 2; size <= circuit_size_lim; size++ )
      {
         for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ++ci )
         {
            this_circuit = (auxipdata->circuits_i)[size][i][ci];
            SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, (auxipdata->circuit_vars)[size][i][ci], &(vars[0])) );
            for( size2 = size; size2 <= circuit_size_lim; size2++ )
            {
               /* symmetry breaking .. */
               if( size2 == size )
                  start = ci+1;
               else
                  start = 0;
               
               for( ci2 = start; ci2 < (auxipdata->n_circuits_i)[size2][i]; ++ci2 )
               {

                  circuit2 = (auxipdata->circuits_i)[size2][i][ci2];

                  /* only consider incomparable pairs of subsets */
                  if( subset_of(size,this_circuit,size2,circuit2) )
                     continue;
                  
                  SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, (auxipdata->circuit_vars)[size2][i][ci2], &(vars[1])) );
                 
                  /* this_circuit and circuit2 are distinct, incomparable and both contain i */
                 
                  /* printf("foo %d %d %d %d %d %d\n", i, size, ci, size2, ci2, size+size2-2); */
                  /* for( j = 0; j < size; j++) */
                  /*    printf("%d,",this_circuit[j]); */
                  /* printf("\n"); */
                  /* for( j = 0; j < size2; j++) */
                  /*    printf("%d,",circuit2[j]); */


                  assert( both != NULL );
                  n_both = myunioni(i,size,this_circuit,size2,circuit2,both);
                  /* printf("\n"); */
                  /* for( j = 0; j < n_both; j++) */
                  /*    printf("%d,",both[j]); */
                  /* printf("\n"); */
                  /* printf("****************\n"); */
                  nvars = 2;

                  /* printf("Working on : %d, %s, %s\n",i,SCIPvarGetName((auxipdata->circuit_vars)[size][i][ci]), */
                  /*    SCIPvarGetName((auxipdata->circuit_vars)[size2][i][ci2])); */
                  
                  /* find all circuits which are subsets of 'both' (which does not have i in it ) */
                  for( j = 0; j < psd->n; ++j )
                  {
                     if( j == i )
                        continue;

                     for( size3 = 2; size3 < min(circuit_size_lim,n_both+1); size3++ )
                        for( ci3 = 0; ci3 < (auxipdata->n_circuits_i)[size3][j]; ++ci3 )
                        {
                           circuit3 = (auxipdata->circuits_i)[size3][j][ci3];

                           /* include each circuit only once */
                           if( circuit3[0] != j )
                              continue;

                           /* printf("Considering: %s\n",SCIPvarGetName((auxipdata->circuit_vars)[size3][j][ci3])); */
                           
                           /* exclude those which are subsets of either 'circuit'
                              or 'circuit2'
                           */
                           if(
                              subset_of(size3,circuit3,size,this_circuit) ||
                              subset_of(size3,circuit3,size2,circuit2) )
                           {
                              /* printf("Rejected as a subset\n"); */
                              continue;
                           }
                           
                           if( subset_of(size3,circuit3,n_both,both) )
                           {
                              vars[nvars++] = (auxipdata->circuit_vars)[size3][j][ci3];
                              assert( nvars <= maxnvars );
                           }
                           /* else */
                           /*    printf("Rejected as NOT a subset\n"); */
                        }
                  }

                  SCIPsnprintf(s, SCIP_MAXSTRLEN, "ce(%s,%s,%s)",psd->nodeNames[i],
                     SCIPvarGetName((auxipdata->circuit_vars)[size][i][ci]),
                     SCIPvarGetName((auxipdata->circuit_vars)[size2][i][ci2]));
                  SCIP_CALL( SCIPcreateConsBasicSetcover(auxipdata->subscip, &cons, s, nvars, vars) );
                  /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
                  SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );
               }
            }
         }
      }
   }
   SCIPfreeMemoryArray(scip, &both);
   SCIPfreeMemoryArray(scip, &vars);
   
   /* START of LP solution-specific code */
               
   /* collect all variables with positive value in the current LP solution */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->family), psd->n) );
   for( i = 0; i < psd->n; ++i )
   {
      assert(solinfo->nposvars[i] > -1);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->family[i]), solinfo->nposvars[i]) );
      for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
      {
         k = solinfo->posvars[i][ki];
         assert(k > -1 && k < psd->nParentSets[i]);

         /* val = SCIPgetSolVal(scip, sol, psd->PaVars[i][k]); */
         val = solinfo->lpsolvals[i][k];
         assert( val > 0 );
         assert(SCIPisPositive(scip, val));

         n_pas = psd->nParents[i][k];
         pas = psd->ParentSets[i][k];
         SCIPsnprintf(s, SCIP_MAXSTRLEN, "I#%d",i);
         for( j = 0; j < n_pas; j++)
         {
            SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "#%d", psd->ParentSets[i][k][j]);
            strcat(s, tmp);
         }
         SCIP_CALL( SCIPcreateVarBasic(auxipdata->subscip, &(auxipdata->family[i][ki]),
               SCIPvarGetName(psd->PaVars[i][k]), 0.0, 1.0, val, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->family[i][ki]) );
         /* SCIP_CALL( SCIPreleaseVar(auxipdata->subscip, &(auxipdata->family[i][ki])) ); */

         
         /* only consider circuits containing i which are small enough
            to potentially be subsets of pas \cup \{i\}
         */
         nvars_max = mypow2(n_pas);
         /* printf("foo %d\n",nvars_max); */
         SCIP_CALL( SCIPallocMemoryArray(scip, &vars, nvars_max) );
         nvars = 0;
         for( size = 2; size < min(circuit_size_lim,n_pas + 1); size++ )
         {
            circuits = (auxipdata->circuits_i)[size][i];
            for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ++ci )
            {
               this_circuit = circuits[ci];
               if( subseti(i,size,this_circuit,n_pas,pas) )
               {
                  vars[nvars++] = (auxipdata->circuit_vars)[size][i][ci];
                  /* printf("nvars %d\n",nvars); */
                  assert( nvars <= nvars_max ); 
               }
            }
         }
         /* use same name 's' for cons and variable */
         
         SCIP_CALL( SCIPcreateConsBasicOr(auxipdata->subscip, &cons, s, auxipdata->family[i][ki], nvars, vars) );
         /* SCIP_CALL(  SCIPprintCons(scip, cons, NULL)  ); */
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, cons) );
         
         SCIPfreeMemoryArray(scip, &vars);
      }
   }

   /* set obj limit for subcip */
   SCIP_CALL_ABORT(SCIPsetObjsense(auxipdata->subscip, SCIP_OBJSENSE_MAXIMIZE));
   SCIP_CALL( SCIPsetObjlimit(auxipdata->subscip, 0) );

   /* SCIP_CALL( SCIPwriteOrigProblem(auxipdata->subscip,NULL,NULL,FALSE)  );  */

   /* finally, solve the subMIP! */
   SCIP_CALL( SCIPsolve(auxipdata->subscip) );


   status = SCIPgetStatus(auxipdata->subscip);

   if( status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_INFORUNBD )
   {
      SCIPdebugMessage("Could not find a matroid cut.\n");
      goto TERMINATE;
   }

   nsols = SCIPgetNSols(auxipdata->subscip);
   if( nsols > 0 && SCIPisFeasLE(auxipdata->subscip, SCIPgetSolOrigObj(auxipdata->subscip, SCIPgetBestSol(auxipdata->subscip)), 0.0) )
   {
      SCIPdebugMessage("Could not find a matroid cut: best objective too low\n");
      goto TERMINATE;
   }

   if( status != SCIP_STATUS_SOLLIMIT && status != SCIP_STATUS_GAPLIMIT && status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_NODELIMIT  && status != SCIP_STATUS_TIMELIMIT )
   {
      SCIPerrorMessage("Solution of subscip for matroid separation returned with invalid status %d.\n", status);
      return SCIP_ERROR;
   }

   subscip_sols = SCIPgetSols(auxipdata->subscip);
   for( si = 0; si <  nsols; ++si )
   {
      subscip_sol = subscip_sols[si];

      /* give up looking for cuts as soon as we come across a solution which is not good enough
         since they are ordered
      */
      if( SCIPisFeasLE(auxipdata->subscip, SCIPgetSolOrigObj(auxipdata->subscip, subscip_sol), 0.0) )
         break;
      
      for( size = 2; size <= circuit_size_lim; size++ )
      {
         for( i = 0; i < psd->n; i++)
         {
            n_set_circuits[size][i] = 0;

            if( SCIPgetSolVal(auxipdata->subscip, subscip_sol, auxipdata->in_ground_set[i]) < 0.5 )
               continue;

            for( ci = 0; ci < (auxipdata->n_circuits_i)[size][i]; ++ci )
            {
               this_circuit = (auxipdata->circuits_i)[size][i][ci];

               if( SCIPgetSolVal(auxipdata->subscip, subscip_sol, (auxipdata->circuit_vars)[size][i][ci]) > 0.5 )
               {
                  set_circuits[size][i][n_set_circuits[size][i]++] = this_circuit;
#ifdef SCIP_DEBUG
               if( this_circuit[0] == i )
                  SCIPdebugMessage("Matroid has this circuit: %s\n", SCIPvarGetName(auxipdata->circuit_vars[size][i][ci]));
#endif
               }
            }
         }
      }
      
      rhs = 0;
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "matroidcut", -SCIPinfinity(scip), 0, FALSE, FALSE, TRUE) );

      for( i = 0; i < psd->n; ++i )
      {
         if( SCIPgetSolVal(auxipdata->subscip, subscip_sol, auxipdata->in_ground_set[i]) < 0.5 )
            continue;
         
         if( SCIPgetSolVal(auxipdata->subscip, subscip_sol, auxipdata->full_graph[i]) > 0.5 )
            rhs++;

         for( k = 0; k < psd->nParentSets[i]; ++k )
         {
            included = FALSE;
            
            for( size = 2; size < min(circuit_size_lim,psd->nParents[i][k]+2); size++ )
            {
               if( included )
                  break;
               
               for(ci = 0; ci < n_set_circuits[size][i]; ci++ )
               {
                  this_circuit = set_circuits[size][i][ci];
                  /* need to check that parent sets are sorted */
                  if( subseti(i,size,this_circuit,psd->nParents[i][k],psd->ParentSets[i][k]) )
                  {
                     SCIP_CALL( SCIPaddVarToRow(scip, cut, psd->PaVars[i][k], 1.0) );
                     included = TRUE;
                     break;
                  }
               }
            }
         }
      }
      SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
      assert(SCIPisIntegral(scip, rhs));
      SCIPdebugMessage(" -> Matroid-cut <matroidcut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
         SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
         SCIPgetCutEfficacy(scip, NULL, cut),
         SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
         SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
      SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));
               
      SCIP_CALL( SCIPaddRow(scip, cut, FALSE, cutoff) );
      if( *cutoff )
      {
         SCIPdebugMessage("Matroid cut led to cutoff\n");
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         return SCIP_OKAY;
      }

      if( SCIPisCutEfficacious(scip, sol, cut) )
         *found_efficacious_ptr = TRUE;

      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }

   for( size = 2; size <= circuit_size_lim; size++ )
   {
      for( i = 0 ; i < psd->n ; ++i )
         SCIPfreeMemoryArray(scip, &(set_circuits[size][i]));
      SCIPfreeMemoryArray(scip, &(set_circuits[size]));
      SCIPfreeMemoryArray(scip, &(n_set_circuits[size]));
   }
   SCIPfreeMemoryArray(scip, &set_circuits);
   SCIPfreeMemoryArray(scip, &n_set_circuits);

   /* free auxipdata */

TERMINATE:
   SCIP_CALL( AuxIPDataFree(scip, auxipdata, psd->n, circuit_size_lim) );
   return SCIP_OKAY;
   
}
   
