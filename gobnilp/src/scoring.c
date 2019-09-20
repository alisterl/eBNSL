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

/** @file scoring.c
 *  @brief  Generates local scores from data
 *  @author James Cussens
 *  @author Mark Bartlett
 */

/* uncomment next line for debugging output */
/*#define SCIP_DEBUG*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include "scoring.h"
#include "scorecache.h"
#ifdef BLAS
#include "bge_matrix.h"
#include "bge_vector.h"
#include "bge_posterior.h"
#include "bge_score.h"
#endif
#include "utils.h"
#include "versiongit.h"
#include "pricer_family.h"

#define BLOCKSIZE 10000
#define MAXARITY UCHAR_MAX


typedef unsigned int ROW;        /**< Index of a row (i.e.\ datapoint) in the data */
typedef unsigned char ARITY;     /**< Arity of a variable in the data */
typedef unsigned char VALUE;     /**< Value of a variable in the data */
typedef double SCOREARG;         /**< An argument for a local score function. */
typedef COUNT* FLATCONTAB;       /**< Flat contingency table */

enum score_code
{
   SCORE_BDEU = 0,   /**< BDeu local scores */
   SCORE_BGE =  1,   /**< BGe local scores */
   SCORE_BIC =  2    /**< BIC local scores */
};
typedef enum score_code SCORE_CODE;           /**< Code for different scoring methods */

/** represents data for a 'query', e.g.\ X1=1,X3=2
 *
 * There is always a count of how many datapoints satisfy the query.
 * If, as is most common,  the highest-indexed variable (in the entire dataset) is not mentioned in the query then:
 * if the count < rmin there is a pointer to the datapoint indices for the data
 * otherwise there is a pointer to an 'array' of 'vary nodes' one for each
 * of the remaining variables.
 */
struct adtree
{
   COUNT count;                /**< how many datapoints for this query */
   struct varynode *children;  /**< one for each variable specialising  the query, if any
                                  (NULL if leaflist is used or there are no specialisations)*/
   ROW *leaflist;              /**< leaflist, if used (NULL otherwise ) */
};
typedef struct adtree ADTREE;  /**< An AD tree */


struct varynode                 /** for splitting data on values of a particular variable (variable not stored in varynode) */
{
   struct adtree **children;    /**< children[val] is a pointer to the ADTREE for the specialisation var=val
                                   (or NULL if no data for this, mcv) */
   VALUE mcv;                   /**< most common value for this variable (in appropriate 'conditional' dataset) */
};
typedef struct varynode VARYNODE; /**< A varynode (in an AD tree) */


union treecontab                 /** A tree-shaped contingency table */
{
   union treecontab *children;   /**< when there are variables...
                                  if treecontab.children == NULL then there are only zero counts in the contingency table.
                                  if treecontab.children != NULL then treecontab.children[i] is the treecontab formed by specialising on
                                  the ith value of the first variable */
   COUNT count;                  /**< when there are no variables treecontab.count just stores a count */
};
typedef union treecontab TREECONTAB; /**< A tree-shaped contingency table. The variables for the contingency table are not stored in
                                    this data structure. */


struct forchild
{
   double neglogaritychild;    /**< -log of the arity of the child */
   double penalty;             /**< log(N)/2 * (arity[child]-1) (for BIC) */
};
typedef struct forchild FORCHILD;

struct scoringparams
{
   SCORE_CODE score_type;  /**< scoring function, e.g. BDeu or BGe */
   SCIP_Bool continuous;   /**< if true then data is continuous otherwise discrete */
   int palim;              /**< limit on the number of parents */
   SCIP_Bool pruning;      /**< whether to prune parent sets with a score lower than one of their subsets */
   SCIP_Real prunegap;     /**< value of gobnilp/scoring/prunegap parameter, see its documentation */
   int initpalim;          /**< initial limit on the number of parents */
   SCIP_Bool pricing;      /**< whether pricing is enabled */
};
typedef struct scoringparams SCORINGPARAMS; /**< scoring parameters */

struct prior
{
   double alpha;             /**< The 'effective sample size' governing the BDeu score */
   ARITY* arity;             /**< arity[i] is the arity of variable i (NULL for purely continuous data) */
   int nvars;                /**< number of BN variables */
   char** nodeNames;         /**< BN variable names */
   char*** labels;           /**< Maps integer encoding of a value to the corresponding string (only for discrete) */
   int** is_parent;          /**< is_parent[i][j] is set to 1 (-1) if j can(not) be a parent of i */
   VARIABLE** is_parents;    /**< is_parents[i] is the set of required parents for i */
   int*  n_is_parents;       /**< n_is_parents[i] is the number of required parents for i */
#ifdef BLAS
   SCIP_Real alpha_mu;       /**< \f$ \alpha_{\mu} \f$ (See <a href="https://projecteuclid.org/euclid.aos/1407420013">Kuipers et al</a>)*/
   int alpha_omega;          /**< \f$ \alpha_{\omega} \f$ (See <a href="https://projecteuclid.org/euclid.aos/1407420013">Kuipers et al</a>)*/
   Bge_Matrix* prior_matrix; /**< \f$ T \f$ (See <a href="https://projecteuclid.org/euclid.aos/1407420013">Kuipers et al</a>)*/
   Bge_Vector* nu;           /**< \f$ \mathbf \nu \f$ (See <a href="https://projecteuclid.org/euclid.aos/1407420013">Kuipers et al</a>)
                                A NULL value indicates that \f$ \mathbf{\nu}=\bar{\mathbf{x}}\f$ where \f$ \bar{\mathbf{x}}\f$
                                is the sample mean */
#endif
};
typedef struct prior PRIOR;  /**< Prior values (arity is prior since independent of oberved data values) */


struct posterior
{
   ADTREE* adtree;                /**< An AD tree (NULL for purely continuous data) */
   VALUE** data;                  /**< data[i][j] is the value of variable i in row j (discrete) */
   int nrows;                     /**< number of datapoints */
   int* num_observed;             /**< Number of distinct values observed in the data for each variable */
   double logndiv2;               /**< log(nrows) / 2, which is used by BIC */
#ifdef BLAS
   Bge_Matrix* continuous_data;   /**< Continuous data (NULL for purely discrete data) */
   Bge_Matrix* posterior_matrix;  /**< The posterior matrix R (NULL for purely discrete data) */
   double log_prefactor;          /**< The ratio of logarithms for the score (ignored for purely discrete data) */
   double* log_gamma_ratio_table; /**< A table that contains all the possible gamma_ratios
                                     for each different size parent set (NULL for purely discrete data) */
#endif
};
typedef struct posterior POSTERIOR;  /**< Posterior values, e.g. the data and values which depend on the data */

struct scoring_extras
{
   SCIP* scip;                    /**< SCIP instance */
   SCIP_Bool fast;                /**< If true then C's lgamma function is used for computing BDeu scores,
                                     otherwise a slower, more accurate sum of logs */
   SCORECACHE* scorecache;        /**< score cache */
   int maxflatcontabsize;         /**< maximum size for a flat contingency table (as opposed to an TREECONTAB tree contingency table )*/
   COUNT npos_cells;              /**< number of cells with a positive count in most recent contingency table used for scoring */
};
typedef struct scoring_extras SCORING_EXTRAS;     /**< Additional information to send to a local scoring function, e.g. a cache */

struct trienode
{
   VARIABLE elt;             /**< last element in the set (this value ignored if node corresponds to empty set ) */
   SCIP_Bool keep;           /**< whether to keep the node to make a family variable */
   SCORE score;              /**< if keep==TRUE then score for this parent set
                                else the score of this parent set's best scoring proper subset */
   SCORE ub;                 /**< upper bound on score of supersets of this parent set */
   struct trienode* mother;  /**< parent node: (or NULL if the root node) */
   struct trienode* child;   /**< leftmost child: extensions to this parent set (or NULL if none) */
   struct trienode* next;    /**< sibling to the right: a parent set of same size whose non-final elements are the same and
                                whose final element is greater than elt (or NULL if none) */
};
typedef struct trienode TRIENODE;  /**< Trie structure for storing scored parent sets */

struct data_etc
{
   PRIOR* prior;
   POSTERIOR* posterior;
   SCORING_EXTRAS* scoring_extras;
   SCORINGPARAMS* scoring_params;
   TRIENODE*** openlist;          /**< openlist[child] is an open list of nodes */
   int* nopenlist;                /**< length of the open list */
};

/* struct indexedcontab */
/* { */
/*    int nvariables;     /\**< number of variables in the contingency table *\/ */
/*    int nrows;          /\**< number of rows in the contingency table *\/ */
/*    COUNT* rows;        /\**< The value of variable i in row j is rows[(nvariables+1)*j+i]  */
/*                         The count associated with row j is rows[(nvariables+1)*j+nvariables] *\/ */
/* }; */
/* typedef struct indexedcontab INDEXEDCONTAB; */

/** Checks whether a string represents a valid arity for a discrete variable
 * @return TRUE if the string represents a valid arity for a discrete variable, else FALSE
 */
static
SCIP_Bool isarity(
   char* str   /**< string to check */
   )
{
   unsigned int i;

   if( str[0] == '0' )
      return FALSE;
   for( i = 0; i < strlen(str); i++ )
      if( !isdigit(str[i]) )
         return FALSE;
   return TRUE;
}

/** Checks whether a string represents a float
 * @return TRUE if the string represents a float, else FALSE
 */
static
SCIP_Bool lookslikefloat(
   char* str   /**< string to check */
   )
{
   unsigned int i;
   SCIP_Bool haspoint = FALSE;
   char c;

   c = str[0];
   if( !( isdigit(c) || c == '-' || c == '.') )
      return FALSE;

   for( i = 1; i < strlen(str); i++ )
   {
      c = str[i];
      if( c == '.' )
         haspoint = TRUE;
      else if( !(isdigit(c)) )
         return FALSE;
   }
   return haspoint;
}


#if 0
/** return whether a cluster has an intersection with the given parent set of at least k
 *
 *  sets represented by ordered sets of integers
 *  @return TRUE if intersection of cluster and parent set is of size at least k
 */
static
SCIP_Bool intersect(
   int* ps,            /**< the parent set */
   int npa,            /**< size of the parent set */
   int* cluster,       /**< the cluster */
   int cluster_size,   /**< size of the cluster */
   int k               /**< k-value of cluster */
   )
{
   int intersection_size = 0;
   int i = 0;
   int j = 0;

   assert( ps != NULL );
   assert( cluster != NULL );
   assert( k > 0 );

   while( i < npa && j < cluster_size )
   {
      if( cluster[j] == ps[i] )
      {
         if( ++intersection_size >= k )
            return TRUE;
         i++;
         j++;
      }
      else if( cluster[j] < ps[i] )
      {
         j++;
      }
      else
         i++;
   }

   return FALSE;
}
#endif


#if 0
/** Computer the 'dual penalty' for a given family (= child and parent set).
 *
 *  This just adds up the dual values for LP rows which
 *  the corresponding family variable would participate in if it were to be created
 */
static
SCIP_Real getDualPenalty(
   int child,             /**< child in family */
   int* ps,               /**< parent set */
   int npa,               /**< size of parent set */
   DUALINFO* dualinfo     /**< info on dual values for current LP */
   )
{
   SCIP_Real dualpenalty;
   int i;

   /* initialise to penalty for all family variables with this child
      which is dual value of constraint stating that there is at most one
      non-empty parent set
   */
   dualpenalty = dualinfo->childpenalties[child];

   /* add in penalties from constraints stating that the arrow variable i<-j
      is an upperbound on any family variable with child i and containing j
      as a parent
   */
   for( i = 0; i < npa; i++)
   {
      dualpenalty += dualinfo->arrowpenalties[dualinfo->n*child+ps[i]];
   }

   /* penalties from k-cluster cuts */
   for( i = 0; i < dualinfo->nclusters[child]; i++)
   {
      CLUSTER_CUT* cluster_cut = dualinfo->clusters[child][i];

      if( intersect(ps,npa,cluster_cut->elts,cluster_cut->nelts,cluster_cut->k) )
         dualpenalty += cluster_cut->dualsol;
   }

   return dualpenalty;
}
#endif

/** Parses a string to extract a set from it.
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE parseSet(
   SCIP* scip,            /**< The SCIP instance that this belongs to. */
   int n,                 /**< The number of variables */
   char** nodeNames,      /**< Variable names */
   const char*  str,      /**< string to parse */
   int*  set,             /**< result set */
   int*  n_set,           /**< length of set */
   SCIP_Bool* success     /**< success flag */
   )
{
   const char* t;
   char tmp[SCIP_MAXSTRLEN];
   char nodename[SCIP_MAXSTRLEN];
   int k;
   int nodeindex;

   *n_set = 0;
   t = str;
   while( *t != '\0' )
   {
      k = 0;
      while( *t != '\0' && *t != ',' )
         tmp[k++] = *t++;
      tmp[k] = '\0';
      if( sscanf(tmp, "%d", &nodeindex) != 1 )
      {
         if( sscanf(tmp, "%s", nodename) != 1 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Expected variable name. Got: %s\n", tmp);
            *success = FALSE;
            return SCIP_OKAY;
         }
         nodeindex = get_index_names(nodename, n, nodeNames);
         if( nodeindex == -1 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Expected variable name. Got: %s\n", tmp);
            *success = FALSE;
            return SCIP_OKAY;
         }
      }

      set[(*n_set)++] = nodeindex;

      if( *t == ',' )
         t++;
   }
   return SCIP_OKAY;
}


/** Enforce or prohibit a parent for a child, returning an error if this leads to
 * infeasibilty
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE setcheck(
   int** is_parent,   /**< is_parent[i][j] = -1/1 if j cannot/must be a parent of i */
   int i,             /**< child variable */
   int j,             /**< parent variable */
   int val            /**< -1 for prohibition, 1 for enforcement */
   )
{
   if(
      (is_parent[i][j] == -1 && val == 1)
      ||
      (is_parent[i][j] == 1 && val == -1)
      )
   {
      SCIPerrorMessage("Infeasibility detected: variable %d is both required and prohibited from being a parent of variable %d\n", j, i);
      return SCIP_ERROR;
   }
   else if( is_parent[i][j] == 1 && is_parent[j][i] == 1)
   {
      SCIPerrorMessage("Infeasibility detected: variable %d is required to be both a parent and a child of variable %d\n", j, i);
      return SCIP_ERROR;
   }
   else
   {
      is_parent[i][j] = val;
      return SCIP_OKAY;
   }
}


/** Use a conditional independence (ci) constraint \f$ A \perp B | S \f$ to rule out parents.
 *  A and B are given as comma separated variable names, S not provided since it would not be used.
 *
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param n The number of variables
 *  @param nodeNames Names of the variables
 *  @param a_str A set
 *  @param b_str B set
 *  @param is_parent is_parent[i][j] is set to 1 (-1) if j can(not) be a parent of i
 *
 *  @return SCIP_OKAY if the constraint on parents was added successfully or an error code otherwise.
 */
static
SCIP_RETCODE ci_constraint(
   SCIP* scip,
   int n,
   char** nodeNames,
   const char* a_str,
   const char* b_str,
   int** is_parent
   )
{
   int i;
   int j;

   int a[SCIP_MAXSTRLEN];
   int b[SCIP_MAXSTRLEN];
   int n_a;
   int n_b;

   SCIP_Bool success;

   success = TRUE;
   SCIP_CALL( parseSet(scip, n, nodeNames, a_str, a, &n_a, &success) );
   if( !success )
      return SCIP_ERROR;
   SCIP_CALL( parseSet(scip, n, nodeNames, b_str, b, &n_b, &success) );
   if( !success )
      return SCIP_ERROR;

   /* If A _|_ B | S rule out edges (in either direction) between
      elements of A and B
   */
   for( i = 0; i < n_a; ++i)
      for( j = 0; j < n_b; ++j)
      {
         if( a[i] == b[j] )
            return SCIP_ERROR;
         SCIP_CALL( setcheck(is_parent,a[i],b[j],-1) );
         SCIP_CALL( setcheck(is_parent,b[i],a[j],-1) );
      }

   return SCIP_OKAY;
}



/** Processes a constraint on the DAG structure.
 *
 *  @param scip The SCIP instance in which to add the constraint.
 *  @param n The number of variables
 *  @param nodeNames Names of the variables
 *  @param line The description of the constraint to add.
 *  @param is_parent is_parent[i][j] is set to 1 (-1) if j must(not) be a parent of i
 *
 *  @return SCIP_OKAY if the constraint was added or an error otherwise.
 */
static
SCIP_RETCODE process_constraint(
   SCIP* scip,
   int n,
   char** nodeNames,
   const char* line,
   int** is_parent
   )
{

   int i;
   int j;
   int child;

   char a_str[SCIP_MAXSTRLEN];
   char b_str[SCIP_MAXSTRLEN];
   char s_str[SCIP_MAXSTRLEN];

   if( line[0] == '#' )
      return SCIP_OKAY;

   if( sscanf(line, "%[^~<>-]-%[^~<>-]", a_str, b_str) == 2 )
   {
      /* enforced edge in undirected skeleton */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      /* do nothing at present */
   }
   else if( sscanf(line, "%[^~<>-]<--%[^~<>-]", a_str, b_str) == 2 )
   {
      /* enforced partial order relation */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      /* do nothing at present */
   }
   else if( sscanf(line, "~%[^~<>-]<--%[^~<>-]", a_str, b_str) == 2 )
   {
      /* prohibited partial order relation */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      /* do nothing at present */
   }
   else if( sscanf(line, "~%[^~<>-]-%[^~<>-]", a_str, b_str) == 2 )
   {
      /* prohibited edge in undirected skeleton */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( setcheck(is_parent,i,j,-1) );
      SCIP_CALL( setcheck(is_parent,j,i,-1) );
   }
   else if( sscanf(line, "%[^~<>-]<-%[^~<>-]", a_str, b_str) == 2 )
   {
      /* enforced arrow */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( setcheck(is_parent,i,j,1) );
   }
   else if( sscanf(line, "~%[^~<>-]<-%[^~<>-]", a_str, b_str) == 2 )
   {
      /* prohibited arrow */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      SCIP_CALL( setcheck(is_parent,i,j,-1) );
   }
   else if( sscanf(line, "%[^~<>-]->%[^~<>-]<-%[^~<>-]", a_str, s_str, b_str) == 3 )
   {
      /* enforced immorality */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      child = get_index_names(s_str, n, nodeNames);
      if( i == -1 || j == -1 || child == -1 )
         return SCIP_READERROR;

      SCIP_CALL( setcheck(is_parent,child,i,1) );
      SCIP_CALL( setcheck(is_parent,child,j,1) );

   }
   else if( sscanf(line, "~%[^~<>-]->%[^~<>-]<-%[^~<>-]", a_str, s_str, b_str) == 3 )
   {
      /* prohibited immorality */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      child = get_index_names(s_str, n, nodeNames);
      if( i == -1 || j == -1 || child == -1 )
         return SCIP_READERROR;

      /* do nothing at present */
   }
   else if( sscanf(line, "%[^_~<>-]_|_%[^|~<>-]|%[^~<>-]", a_str, b_str, s_str) == 3 )
      SCIP_CALL( ci_constraint(scip, n, nodeNames, a_str, b_str, is_parent) );
   else if( sscanf(line, "%[^_~<>-]_|_%[^~<>-]", a_str, b_str) == 2 )
      SCIP_CALL( ci_constraint(scip, n, nodeNames, a_str, b_str, is_parent) );
   else if( sscanf(line, "%[^~<>-]<%[^~<>-]", a_str, b_str) == 2 )
   {
      /* enforced total order relation */
      i = get_index_names(a_str, n, nodeNames);
      j = get_index_names(b_str, n, nodeNames);
      if( i == -1 || j == -1 )
         return SCIP_READERROR;

      /* just rule out j from being a parent of i at present */
      SCIP_CALL( setcheck(is_parent,i,j,-1) );
   }
   else if ( sscanf(line, "sink%s", a_str) == 1 )
   {
      /* enforced sink */
      i = get_index_names(a_str, n, nodeNames);
      if( i == -1 )
         return SCIP_READERROR;

      for( j = 0; j < n; ++j )
         if( j != i )
            SCIP_CALL( setcheck(is_parent,j,i,-1) );
   }
   else
   {
      SCIPerrorMessage("Not recognised as a DAG constraint: %s\n", line);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}


/** Reads in structural constraints to restrict which parent sets are considered
 *
 *  (Altered version of function of same name in probdata_bn.c)
 *  @param scip The SCIP instance to add the constraint to.
 *  @param n The number of variables
 *  @param nodeNames Names of the variables
 *  @param is_parent is_parent[i][j] = 1 (-1) if j must (not) be a parent of i
 *  @return SCIP_OKAY if the constraints could be successfully added.
 */
static
SCIP_RETCODE addGeneralDAGConstraints(
   SCIP* scip,
   int n,
   char** nodeNames,
   int** is_parent
   )
{
   int status;
   char s[SCIP_MAXSTRLEN];
   char* dagconstraintsfile;
   FILE* dagconstraints;

   SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/dagconstraintsfile", &dagconstraintsfile) );

   if( strcmp(dagconstraintsfile, "") == 0 )
      return SCIP_OKAY;

   dagconstraints = fopen(dagconstraintsfile, "r");
   if( dagconstraints == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", dagconstraintsfile);
      return SCIP_NOFILE;
   }

   status = fscanf(dagconstraints, "%[^\n]%*c", s);
   while( status == 1 )
   {
      SCIP_CALL( process_constraint(scip, n, nodeNames, s, is_parent) );
      status = fscanf(dagconstraints, "%[^\n]%*c", s);
   }

   status = fclose(dagconstraints);
   if( status != 0 )
   {
      SCIPerrorMessage("Could not close file %s.\n", dagconstraintsfile);
      return SCIP_ERROR;
   }
   return SCIP_OKAY;
}


/** get the node corresponding to a non-empty extension of a given
 * set or NULL if the extension is not stored
 * @return node corresponding to extension or NULL if it does not exist
 */
static
TRIENODE* getnode(
   const TRIENODE* current,  /**< set to be extended */
   const VARIABLE* elts,     /**< elements in the extension */
   int nelts                 /**< number of elements in the extension */
   )
{

   assert( current != NULL );
   assert( elts != NULL );
   assert( nelts > 0);

   current = current->child;
   while( current != NULL && current->elt < elts[0] )
      current = current->next;

   if( current == NULL || current->elt != elts[0] )
      return NULL;
   else if( nelts == 1)
      return (TRIENODE*) current;
   else
      return getnode(current,elts+1,nelts-1);
}

#if 0
/** get the size of the parent set corresponding to a trie node
 * @return the size of the parent set corresponding to the given trie node
 */
static
int getnpa(
   TRIENODE* node /**< the trie node */
   )
{
   if( node->mother == NULL )
      return 0;
   else
      return getnpa(node->mother) + 1;
}
#endif

/** get the set corresponding to a node.
 *  space for set must have already been allocated
 */
static
void getset(
   const TRIENODE* current,  /**< node */
   VARIABLE* elts,           /**< set (to populate) */
   int depth                 /**< depth of node (= size of set) */
   )
{
   int i;

   assert( current != NULL );
   assert( elts != NULL );

   for( i=depth-1; i >= 0; i--)
   {
      elts[i] = current->elt;
      assert(current->mother != NULL);
      current = current->mother;
   }
}


#ifdef SCIP_DEBUG
/** print a TRIENODE in a 'raw' format (for debugging ony ) */
static
void printnode(
   const TRIENODE* node  /**< TRIENODE to be printed */
   )
{
   printf("self=%p,elt=%d,keep=%d,score=%g,mother=%p,child=%p,next=%p\n",
      (void*)node,node->elt,node->keep,node->score,(void*)node->mother,(void*)node->child,(void*)node->next);
}
#endif

/** Find best score and lower upper bound over each smaller-by-1 subset of a new parent set:  current + { newpa }.
 *
 * If any such subsets are not stored (in the trie) then they must have been exponentially pruned
 * @return FALSE if all such subsets are stored else TRUE
 */
static
SCIP_Bool exppruned(
   SCIP* scip,               /**< SCIP instance */
   const TRIENODE* current,  /**< trie node corresponding to the old parent set  */
   int size,                 /**< size of old parent set */
   VARIABLE newpa,           /**< new parent to be added to old parent set */
   SCIP_Real prunegap,
   SCORE* bestscore,         /**< returns the best score over smaller-by-1 subsets of a new parent set current + newpa */
   SCORE* lowest_ub          /**< returns the lowest upper bound over smaller-by-1 subsets of a new parent set current + newpa */
   )
{
   int i;
   VARIABLE* subset;
   TRIENODE* node;

   SCIP_Bool result = FALSE;

   /* initialise to the 'score' of current
      N.B. this could actually be the score of one of its subsets
   */
   *bestscore = current->score;

   /* initialise lowest upper bound to be the upper bound associated with 'current'
      so initially, at least, the lowest upper bound is not lower than the best score
   */
   *lowest_ub = current->ub;
   assert(*lowest_ub >= *bestscore + prunegap);

   /* suppose current = {1,3,4} and npa = 5
      make subset = {1,3,4,5}
   */
   SCIP_CALL( SCIPallocBufferArray(scip, &subset, size+1) );
   getset(current,subset,size);
   subset[size] = newpa;

   /* find all subsets of size=size by going up trie */
   for(i = size; i > 0; i--)
   {
      current = current->mother;
      assert(current != NULL);
      node = getnode(current,subset+i,(size-i+1));
      if( node == NULL )
      {
         /* node missing so none of its supersets
            (including: current \cup \{ newpa \})
            should be considered
         */
         result = TRUE;
         break;
      }
      *bestscore = MAX(*bestscore,node->score);
      *lowest_ub = MIN(*lowest_ub,node->ub);
      if( *lowest_ub < *bestscore + prunegap )
      {
         /* score for (current \cup \{ newpa \}) and
            ALL ITS SUPERSETS are upper bounded by lowest_ub
            which is below the score of one of its subsets, so
            none of these supersets should be considered
         */
         result = TRUE;
         break;
      }
   }

   SCIPfreeBufferArray(scip, &subset);
   return result;
}

/** create a trie node, corresponding to a parent set, and set its members
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE makenode(
   SCIP* scip,           /**< SCIP instance */
   TRIENODE** nodeptr,   /**< (pointer to) trie node to create */
   VARIABLE elt,         /**< element (BN variable) to associate with node
                            (last element in parent set associated with node) */
   SCIP_Bool keep,       /**< whether parent set should be kept to make a family variable from */
   SCORE score,          /**< score of highest-scoring subset (proper or improper) of this parent set */
   SCORE ub,             /**< upper bound on scores of supersets of this parent set */
   TRIENODE* mother      /**< trie node for parent set produced by removing highest-indexed variable in current parent set */
   )
{

   SCIP_CALL( SCIPallocBlockMemory(scip, nodeptr) );
   (*nodeptr)->elt = elt;
   (*nodeptr)->keep = keep;
   (*nodeptr)->score = score;
   (*nodeptr)->ub = ub;
   (*nodeptr)->mother = mother;
   (*nodeptr)->child = NULL;
   (*nodeptr)->next = NULL;

   return SCIP_OKAY;
}

/** compare two trie nodes by score
 * @return -1,0,1 depending on whether 1st node score is bigger, equal, less than (resp.) 2nd node score
 */
static
SCIP_DECL_SORTPTRCOMP(nodeScoreComp)
{
   TRIENODE* node1;  /**< 1st tree node */
   TRIENODE* node2;  /**< 2nd tree node */

   node1 = (TRIENODE*)elem1;
   node2 = (TRIENODE*)elem2;

   if(node1->score > node2->score)
      return 1;
   else if (node1->score < node2->score)
      return -1;
   else
      return 0;
}

/** set quantities specific to a child */
static
void init_forchild(
   FORCHILD* forchild,      /**< stores child-specific info */
   SCORE_CODE score_type,   /**< score e.g. BDeu, BGe or BIC */
   PRIOR* prior,            /**< prior quantities */
   POSTERIOR* data,         /**< posterior quantities, e.g. the data */
   int child                /**< the child */
   )
{
   if( score_type == SCORE_BDEU )
      forchild->neglogaritychild = -log( (double) prior->arity[child]);
   if( score_type == SCORE_BIC )
      forchild->penalty = data->logndiv2 * ( (int) (prior->arity)[child] - 1);
}

static void build_varynode(VARYNODE *varynode, VARIABLE variable, ROW *theserows, COUNT count, int rmin, const int depth, int *n_nodes, const int adtreedepthlim, const int adtreenodeslim, VARIABLE nvars, VALUE **data, ARITY *arity);

/** Build an AD tree from (a subset of) the data */
static
void build_adtree(
   ADTREE *adtree,            /**< pointer to ADTREE being built */
   const VARIABLE variable,   /**< first variable to specialise further on, if variable=nvars then there is none */
   ROW *theserows,            /**< datapoint indices for this tree */
   COUNT count,               /**< number of datapoints for this tree */
   const int rmin,            /**< if count below this then create a leaflist */
   const int depth,           /**< the depth of this node */
   int *n_nodes,              /**< (pointer to) the number of nodes in the ADTREE */
   const int adtreedepthlim,  /**< limit on the depth of the ADTREE */
   const int adtreenodeslim,  /**< limit on the number of nodes in the ADTREE */
   const VARIABLE nvars,      /**< Number of variables in the data */
   VALUE **data,              /**< data[i][j] is the value of variable i in row j */
   ARITY *arity               /**< arity[i] is the arity of variable i, */
)
{
   COUNT j;
   VARIABLE var;

   assert(variable < nvars + 1);
   assert(count > 0);
   assert(theserows != NULL);
   assert(adtree != NULL);

   adtree->count = count;
   adtree->leaflist = NULL;
   adtree->children = NULL;

   /* if there can be no further splitting just record count */
   if( variable < nvars )
   {
      /* if count small enough then make a leaflist which is a copy of theserows */
      /* if depth too large or number of nodes too large, similarly just dump records in a leaflist */
      if( (int) count < rmin || depth > adtreedepthlim || *n_nodes > adtreenodeslim )
      {
         adtree->leaflist = (ROW *) malloc(count * sizeof(ROW));
         for( j = 0; j < count; ++j )
            adtree->leaflist[j] = theserows[j];
      }
      /* or create vary nodes - one for each further variable - and recurse */
      else
      {
         adtree->children = (VARYNODE *) malloc((nvars - variable) * sizeof(VARYNODE));
         if( adtree->children == NULL )
            printf("Couldn't allocate memory for vary nodes\n");
         for( var = variable; var < nvars; ++var )
            build_varynode((adtree->children) + (var - variable), var, theserows, count, rmin, depth, n_nodes, adtreedepthlim, adtreenodeslim, nvars, data, arity);
      }
   }

   /* can always free since data indices are always copied */
   free(theserows);
   return;
}
/** Build an vary node from (a subset of) the data */
static
void build_varynode(
   VARYNODE *varynode,        /**< varynode being built */
   VARIABLE variable,         /**< which variable is being split */
   ROW *theserows,            /**< datapoint indices to divide between values of variable */
   COUNT count,               /**< number of datapoints for this tree */
   const int rmin,            /**< if count below this then create a leaflist */
   int depth,                 /**< the depth of this node */
   int *n_nodes,              /**< (pointer to) the number of nodes in the ADTREE */
   const int adtreedepthlim,  /**< limit on the depth of the ADTREE */
   const int adtreenodeslim,  /**< limit on the number of nodes in the ADTREE */
   VARIABLE nvars,            /**< Number of variables in the data */
   VALUE **data,              /**< data[i][j] is the value of variable i in row j */
   ARITY *arity               /**< arity[i] is the arity of variable i, */
)
{

   const VALUE *thisdata = data[variable];
   const ARITY thisarity = arity[variable];
   ROW **childdata;
   COUNT *childcount;
   VALUE val;
   COUNT j;
   VALUE mcv = 0;
   COUNT countmcv = 0;
   ROW row;

   assert(variable < nvars);
   assert(varynode != NULL);
   assert(theserows != NULL);
   assert(count > 0);


   /* initialise data structures for splitting data on values of the variable */
   childdata = (ROW **) malloc(thisarity * sizeof(ROW *));
   if( childdata == NULL )
      printf("Couldn't allocate childdata\n");
   childcount = (COUNT *) malloc(thisarity * sizeof(COUNT));
   if( childcount == NULL )
      printf("Couldn't allocate childcount\n");

   for( val = 0; val < thisarity; ++val )
   {
      /* lazily allocate space of size 'count' for each val
         (which is certainly big enough), perhaps should
         do allocate in small blocks, on demand
      */
      childdata[val] = (ROW *) malloc(count * sizeof(ROW));
      if( childdata[val] == NULL )
         printf("Couldn't allocate childdata_val\n");

      childcount[val] = 0;
   }

   /* split the data for this tree on values of the variable */
   for( j = 0; j < count; ++j )
   {
      row = theserows[j];
      val = thisdata[row];
      childdata[val][childcount[val]++] = row;
   }


   /* find most common value */
   for( val = 0; val < thisarity; ++val )
      if( childcount[val] > countmcv )
      {
         countmcv = childcount[val];
         mcv = val;
      }
   assert(countmcv > 0);
   varynode->mcv = mcv;

   /* throw away rows for mcv and any zero counts resize the others */
   /* resize as soon as possible */
   for( val = 0; val < thisarity; ++val )
   {
      if( val == mcv || childcount[val] == 0 )
         free(childdata[val]);
      else
      {
         childdata[val] = (ROW *) realloc(childdata[val], childcount[val] * sizeof(ROW));
         if( childdata[val] == NULL )
            printf("Couldn't re-allocate childdata_val\n");
      }
   }

   varynode->children = (ADTREE **) malloc(thisarity * sizeof(ADTREE *));
   if( varynode->children == NULL )
      printf("Couldn't allocate memory for AD trees\n");

   variable++;          /* can lead to variable=nvars, ie a fake 'extra' variable */
   for( val = 0; val < thisarity; ++val )
   {
      if( val == mcv || childcount[val] == 0 )
         varynode->children[val] = NULL;
      else
      {
         /* childdata[val] freed in build_adtree (unless it becomes a leaflist) */
         varynode->children[val] = (ADTREE *) malloc(sizeof(ADTREE));
         (*n_nodes)++;
         if( varynode->children[val] == NULL )
            printf("Couldn't allocate memory for AD tree\n");

         build_adtree(varynode->children[val], variable, childdata[val], childcount[val],
            rmin, depth + 1, n_nodes, adtreedepthlim, adtreenodeslim, nvars, data, arity);
      }
   }

   free(childdata);
   free(childcount);

   return;

}
/** Delete an AD tree */
static
void delete_adtree(
   ADTREE *adtree,          /**< pointer to ADTREE being deleted */
   const VARIABLE variable, /**< first variable to specialise further on, if variable=nvars then there is none */
   VARIABLE nvars,          /**< Number of variables in the data */
   ARITY *arity             /**< arity[i] is the arity of variable i, */
)
{

   VARIABLE var;
   VARYNODE *vn_ptr = adtree->children;
   VARYNODE vn;
   VALUE val;

   if( adtree->leaflist != NULL )
   {
      assert(vn_ptr == NULL);
      free(adtree->leaflist);
   }
   else if( vn_ptr != NULL )
   {
      assert(adtree->leaflist == NULL);
      for( var = variable; var < nvars; ++var )
      {
         vn = vn_ptr[var - variable];
         for( val = 0; val < arity[var]; ++val )
            if( vn.children[val] != NULL )
               delete_adtree(vn.children[val], var + 1, nvars, arity);
         free(vn.children);
      }
      free(vn_ptr);
   }
   free(adtree);
   return;
}

/** Construct a tree-shaped contingency table from a leaflist
 *  ( contingency table must be for at least one variable )
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE maketreecontableaf(
   SCIP* scip,                 /**< SCIP instance */
   const ROW *leaflist,        /**< datapoints for this query  */
   const COUNT count,          /**< number of datapoints (in leaflist) */
   const VARIABLE *variables,  /**< variables in the contingency table (sorted) */
   int nvariables,             /**< number of variables in the contingency table */
   TREECONTAB* treecontab,     /**< (pointer to ) returned contingency table */
   VALUE **data,               /**< data[i][j] is the value of variable i in row j */
   const ARITY *arity          /**< arity[i] is the arity of variable i, */
)
{

   VARIABLE firstvar;
   VALUE *firstvardata;
   ARITY firstarity;

   VALUE val;
   int valj;

   ROW **bin;
   COUNT *bin_size;

   COUNT j;
   ROW row;

   TREECONTAB* adpt;
   TREECONTAB* adpt2;
   TREECONTAB* adpt3;

   VARIABLE secondvar;
   VALUE *secondvardata;
   ARITY secondarity;
   VALUE val2;

   VARIABLE thirdvar;
   VALUE *thirdvardata;
   ARITY thirdarity;
   VALUE val3;

   SCIP_Bool* poscount;
   VALUE firstvardatarow;

   assert(leaflist != NULL);
   assert(treecontab != NULL);
   assert(variables != NULL);
   assert(data != NULL);
   assert(arity != NULL);
   assert(nvariables > 0);

   firstvar = variables[0];
   firstvardata = data[firstvar];
   firstarity = arity[firstvar];

   if( count == 0 )
   {
      treecontab->children = NULL;
      return SCIP_OKAY;
   }

   treecontab->children = (TREECONTAB *) malloc(firstarity * sizeof(TREECONTAB));
   assert(treecontab->children != NULL);
   adpt = treecontab->children;

   /* specialising for one variable is a big win re performance
    ( since no memory allocation )
   */
   if( nvariables == 1 )
   {
      for( val = 0; val < firstarity; ++val )
         (adpt++)->count = 0;
      adpt = treecontab->children;

      for( j = 0; j < count; ++j )
         ((adpt + firstvardata[leaflist[j]])->count)++;

      return SCIP_OKAY;
   }

   if( nvariables == 2 )
   {
      secondvar = variables[1];
      secondvardata = data[secondvar];
      secondarity = arity[secondvar];

      SCIP_CALL( SCIPallocBufferArray(scip, &poscount, firstarity) );
      /* poscount = (SCIP_Bool *) malloc(firstarity * sizeof(SCIP_Bool)); */
      assert(poscount != NULL);

      for( val = 0; val < firstarity; ++val )
      {
         adpt->children = (TREECONTAB *) malloc(secondarity * sizeof(TREECONTAB));
         assert(adpt->children != NULL);
         adpt2 = adpt->children;
         for( val2 = 0; val2 < secondarity; ++val2 )
            (adpt2++)->count = 0;
         adpt++;
         poscount[val] = FALSE;
      }
      adpt = treecontab->children;

      for( j = 0; j < count; ++j )
      {
         row = leaflist[j];
         firstvardatarow = firstvardata[row];
         poscount[firstvardatarow] = TRUE;
         (((adpt+firstvardatarow)->children)+secondvardata[row])->count++;
      }

      for( val = 0; val < firstarity; ++val )
      {
         if( !poscount[val] )
         {
            free(adpt->children);
            adpt->children = NULL;
         }
         adpt++;
      }
      SCIPfreeBufferArray(scip,&poscount);
      /* free(poscount); */

      return SCIP_OKAY;
   }

   if( nvariables == 3 )
   {
      /* COUNT*** tmpcount; */

      secondvar = variables[1];
      secondvardata = data[secondvar];
      secondarity = arity[secondvar];

      thirdvar = variables[2];
      thirdvardata = data[thirdvar];
      thirdarity = arity[thirdvar];

      SCIP_CALL( SCIPallocBufferArray(scip, &poscount, firstarity) );
      /* SCIP_CALL( SCIPallocBufferArray(scip, &tmpcount, firstarity) ); */
      /* poscount = (SCIP_Bool *) malloc(firstarity * sizeof(SCIP_Bool)); */
      assert(poscount != NULL);

      for( val = 0; val < firstarity; ++val )
      {
         adpt->children = (TREECONTAB *) malloc(secondarity * sizeof(TREECONTAB));
         /* SCIP_CALL( SCIPallocBufferArray(scip, &(tmpcount[val]), secondarity) ); */
         assert(adpt->children != NULL);
         adpt2 = adpt->children;
         for( val2 = 0; val2 < secondarity; ++val2 )
         {
            adpt2->children = (TREECONTAB *) malloc(thirdarity * sizeof(TREECONTAB));
            /* SCIP_CALL( SCIPallocBufferArray(scip, &(tmpcount[val][val2]), thirdarity) ); */
            assert(adpt2->children != NULL);
            adpt3 = adpt2->children;
            for( val3 = 0; val3 < thirdarity; ++val3 )
            {
               /* tmpcount[val][val2][val3] = 0; */
               (adpt3++)->count = 0;
            }
            adpt2++;
         }
         adpt++;
         poscount[val] = FALSE;
      }
      adpt = treecontab->children;

      for( j = 0; j < count; ++j )
      {
         row = leaflist[j];
         firstvardatarow = firstvardata[row];
         poscount[firstvardatarow] = TRUE;
         (((((adpt+firstvardata[row])->children)+secondvardata[row])->children)+thirdvardata[row])->count++;
         /* tmpcount[firstvardata[row]][secondvardata[row]][thirdvardata[row]]++; */
      }

      /* for( val = 0; val < firstarity; ++val ) */
      /* { */
      /*    for( val2 = 0; val2 < secondarity; ++val2 ) */
      /*    { */
      /*       for( val3 = 0; val3 < thirdarity; ++val3 ) */
      /*       { */
      /*          (((((adpt+val)->children)+val2)->children)+val3)->count =  tmpcount[val][val2][val3]; */
      /*       } */
      /*       SCIPfreeBufferArray(scip, &(tmpcount[val][val2])); */
      /*    } */
      /*    SCIPfreeBufferArray(scip, &(tmpcount[val])); */
      /* } */
      /* SCIPfreeBufferArray(scip, &tmpcount); */


      for( val = 0; val < firstarity; ++val )
      {
         if( !poscount[val] )
         {
            adpt2 = adpt->children;
            for( val2 = 0; val2 < secondarity; ++val2 )
            {
               free(adpt2->children);
               adpt2++;
            }

            free(adpt->children);
            adpt->children = NULL;
         }
         adpt++;
      }
      SCIPfreeBufferArray(scip,&poscount);
      /* free(poscount); */


      return SCIP_OKAY;
   }

   /* initialise temporary space
      one bin for each value of the first variable */
   SCIP_CALL( SCIPallocBufferArray(scip, &bin, firstarity) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bin_size, firstarity) );
   /* bin = (ROW **) malloc(firstarity * sizeof(ROW *)); */
   /* bin_size = (COUNT *) malloc(firstarity * sizeof(COUNT)); */
   assert(bin != NULL);
   assert(bin_size != NULL);

   /* use valj to prevent warnings */
   for( valj = firstarity-1; valj >= 0; --valj )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(bin[valj]), count) );
      /* bin[val] = (ROW *) malloc(count * sizeof(ROW)); */
      assert(bin[valj] != NULL);
      bin_size[valj] = 0;
   }

   /* assign each datapoint according to its value for the first variable */
   for( j = 0; j < count; ++j )
   {
      /* row = leaflist[j]; */
      row = *(leaflist++);
      val = firstvardata[row];
      bin[val][bin_size[val]++] = row;
   }

   nvariables--;
   variables++;
   for( val = 0; val < firstarity; ++val )
   {
      SCIP_CALL( maketreecontableaf(scip, bin[val], bin_size[val], variables, nvariables, adpt++, data, arity) );
      SCIPfreeBufferArray(scip, &(bin[val]));
      /* free(bin[val]); */
   }

   SCIPfreeBufferArray(scip, &bin_size);
   SCIPfreeBufferArray(scip, &bin);

   /* free(bin); */
   /* free(bin_size); */

   return SCIP_OKAY;
}

#if 0
/** print a tree-shaped contingency table (only for debugging)
 * pointer version
*/
static void print_treecontabptr(
   const TREECONTAB* treecontab, /**< ( pointer to a ) tree-shaped contingency table */
   const VARIABLE* variables,    /**< variables for the contingency table */
   COUNT nvariables,             /**< number of variables in the contingency table */
   const ARITY* arity            /**< arity[i] is arity of variable i */
   )
{
   VALUE val;

   assert(nvariables == 0 || variables != NULL);
   assert(arity != NULL);

   /* printf("fofof %d\n",nvariables); */

   if( nvariables == 0 )
   {
      printf("count = %d\n", treecontab->count);
   }
   else
   {
      if( treecontab->children == NULL )
      {
         printf("NULL\n");
      }
      else
      {
         for(val = 0; val < arity[variables[0]]; val++ )
         {
            printf("X%d=%d,",variables[0],val);
            fflush(stdout);
            print_treecontabptr(treecontab->children+val,variables+1,nvariables-1,arity);
         }
      }
   }
}
#endif


#if 0
/** In a tree-shaped contingency table, reorder the variables
 * (and thus levels of the tree) by moving a given variable to a given
 * later position in the order
 */
static void pushvar_treecontab_reccall(
   TREECONTAB* treecontab,     /**< (pointer to) tree-shaped contingency table */
   COUNT nvariables,           /**< number of variables in the contingency table */
   ARITY arity0,               /**< arity of the the first (ie 'top') variable */
   const ARITY* arity,         /**< arity[i] is the arity of (i+1)th variable */
   int frompos,                /**< the variable in position (ie layer) frompos is the variable to move */
   int topos                   /**< final position for the variable to move  */
   )
{

   VALUE val0;
   int arity1;
   int val1;
   TREECONTAB* above;

   /* printf("nvars=%d,frompos=%d,topos=%d,%p,%p\n",nvariables,frompos,topos,(void*)treecontab,(void*)treecontab->children); */

   /* just do necessary checks via assert statements */
   assert(treecontab != NULL);
   /* top-level call should never have treecontab->children == NULL
      and neither should any recursive call due to (*GUARD*) below */
   assert(treecontab->children != NULL);
   assert(arity != NULL);

   assert(-1 < frompos);
   assert(frompos < topos);
   assert(topos < (int) nvariables);


   if( frompos != 0 )
   {
      for(val0 = 0; val0 < arity0; val0++ )
      {
         if( treecontab->children[val0].children != NULL ) /* (*GUARD*) */
            pushvar_treecontab_reccall(treecontab->children+val0,nvariables-1,arity[0],arity+1,frompos-1,topos-1);
      }
   }
   else
   {
      /* start by swapping variables in positions 0 and 1 */
      SCIP_Bool counts;

      counts = (nvariables == 2);
      arity1 = arity[0];
      above = (TREECONTAB *) malloc(arity1 * sizeof(TREECONTAB));
      for( val1 = 0; val1 < arity1; val1++ )
      {
         SCIP_Bool all_null_or_zero = TRUE;
         above[val1].children = (TREECONTAB *) malloc(arity0 * sizeof(TREECONTAB));

         for( val0 = 0; val0 < arity0; val0++ )
         {
            if( treecontab->children[val0].children != NULL )
            {
               TREECONTAB tmp;
               tmp = treecontab->children[val0].children[val1];
               /* no need to check "counts" here */
               above[val1].children[val0] = tmp;

               if( all_null_or_zero && ((!counts && tmp.children != NULL) || (counts && tmp.count != 0)))
                  all_null_or_zero = FALSE;
            }
            else
            {
               if( counts )
                  above[val1].children[val0].count = 0;
               else
                  above[val1].children[val0].children = NULL;
            }
         }

         if( all_null_or_zero )
         {
            /* sub-tree at val1 is all zeroes
               so just set to NULL */
            free(above[val1].children);
            above[val1].children = NULL;
         }
      }
      /* get rid of old branches */
      for( val0 = 0; val0 < arity0; val0++ )
      {
         if( treecontab->children[val0].children != NULL )
            free(treecontab->children[val0].children);
      }
      free(treecontab->children);

      /* and replace with new sub-tree */
      treecontab->children = above;

      if( topos > 1 )
      {
         for( val1 = 0; val1 < arity1; val1++ )
         {
            if( treecontab->children[val1].children != NULL ) /* (*GUARD*) */
            {
               /* printf("calling with val = %d <%p>\n", val1, (void*)treecontab->children[val1].children); */
               pushvar_treecontab_reccall(treecontab->children+val1,nvariables-1,arity0,arity+1,0,topos-1);
            }
            /* else */
            /* { */
            /*    printf("val = %d was NULL\n", val1); */
            /* } */
         }
      }
   }
}
#endif

#if 0
/** add up all the counts of the last variable in a contingency table for a given inst of the instvars */
static
void counts_inst_treecontab(
   TREECONTAB treecontab,      /**< tree-shaped contingency table */
   VARIABLE* variables,        /**< variables for the contingency table */
   COUNT nvariables,           /**< number of variables in the contingency table */
   const ARITY* arity,         /**< arity[i] is arity of variable i */
   VARIABLE* instvars,         /**< instantiated variables (must be ordered) */
   VALUE* inst,                /**< the instantiation */
   COUNT* lastvarcounts        /**< lastvarcounts[val] will be count for value val for last variable in contingency table
                                  before initial call to this function lastvarcounts to be initialised to all zeroes */
   )
{
   VALUE val;

   assert(variables != NULL);
   assert(arity != NULL);
   assert(nvariables > 0);

   if( treecontab.children == NULL )
      /* all zero, so nothing to add */
      return;

   if( nvariables == 1 )
   {
      /* add counts for a given full instantiation */
      for(val = 0; val < arity[variables[0]]; val++ )
         lastvarcounts[val] += treecontab.children[val].count;
      return;
   }

   if( variables[0] == instvars[0] )
   {
      /* just add counts for specified value of the instantiated variable */
      counts_inst_treecontab(treecontab.children[inst[0]], variables+1, nvariables-1, arity, instvars+1, inst+1,
         lastvarcounts);
   }
   else
   {
      for(val = 0; val < arity[variables[0]]; val++ )
      {
         /* variables[0] is being summed out so add counts for all its values */
         counts_inst_treecontab(treecontab.children[val], variables+1, nvariables-1, arity, instvars, inst,
            lastvarcounts);
      }
   }
}
#endif

#if 0
/** find the (lexicographically) first instantation of a given subset of variables that has a positive count
 *  in the given contingency table.
 *  returns FALSE if there is no such instantiation
 * @return TRUE if there is an instantiation of given variables with a positive count, else FALSE
 */
static
SCIP_Bool firstposinst_treecontab(
   TREECONTAB treecontab,      /**< tree-shaped contingency table */
   VARIABLE* variables,        /**< variables for the contingency table */
   COUNT nvariables,           /**< number of variables in the contingency table */
   const ARITY* arity,         /**< arity[i] is arity of variable i */
   VARIABLE* instvars,         /**< variables 'to be instantiated' (must be ordered) */
   COUNT ninstvars,            /**< number of variables 'to be instantiated' */
   VALUE* firstinst            /**< returned first instantiation */
   )
{
   VALUE val;
   SCIP_Bool ok;
   SCIP_Bool isinstvar;

   assert(variables != NULL);
   assert(arity != NULL);
   assert(instvars != NULL);
   /* must be at least 2 variables in the contingency table */
   assert(nvariables > 1);
   assert(ninstvars > 0);
   assert(nvariables >= ninstvars);

   isinstvar = ( variables[0] == instvars[0] );
   for(val = 0; val < arity[variables[0]]; val++ )
   {
      if( treecontab.children[val].children != NULL )
      {
         if( isinstvar )
         {
            firstinst[0] = val;
            ok = firstposinst_treecontab(treecontab.children[val], variables+1, nvariables-1, arity, instvars+1, ninstvars-1, firstinst+1);
         }
         else
         {
            ok = firstposinst_treecontab(treecontab.children[val], variables+1, nvariables-1, arity, instvars, ninstvars, firstinst);
         }
         if( ok )
            return TRUE;
      }
   }
   return FALSE;
}
#endif


#if 0
/** compute the bound based on MLEs from a tree-shaped contingency table where the child is at the
 * bottom
 * @return MLE-based bound
 */
static
SCORE mlebound_treecontab(
   TREECONTAB treecontab,      /**< tree-shaped contingency table */
   VARIABLE* variables,        /**< variables for the contingency table */
   COUNT nvariables,           /**< number of variables in the contingency table */
   const ARITY* arity,         /**< arity[i] is arity of variable i */
   SCIP_Bool* parent           /**< parent set in a bitset style representation (all other variables are marginalised over ) */
   )
{

   SCORE mlebound;
   VALUE val;

   assert(arity != NULL);
   assert(nvariables > 0);

   if( treecontab.children == NULL )
      return 0.0;

   if( nvariables > 1)
   {
      mlebound = 0.0;
      if( parent[variables[0]] )
      {
         for(val = 0; val < arity[variables[0]]; val++ )
            mlebound += mlebound_treecontab(treecontab.children[val], variables+1, nvariables-1, arity, parent);
      }
      return 0.0;
   }
   else
   {
      assert(nvariables == 1);

   }


}
#endif


#if 0
/** In a tree-shaped contingency table, reorder the variables
 * (and thus levels of the tree) by moving a given variable to a given
 * later position in the order
 */
static void pushvar_treecontab(
   TREECONTAB* treecontab,     /**< (pointer to) tree-shaped contingency table */
   VARIABLE* variables,        /**< variables for the contingency table */
   COUNT nvariables,           /**< number of variables in the contingency table */
   const ARITY* arity,         /**< arity[i] is arity of variable i */
   int frompos,                /**< the variable in position (ie layer) frompos is the variable to move */
   int topos                   /**< final position for the variable to move  */
   )
{

   ARITY* these_arities;
   int i;
   VARIABLE tmpvar;

   assert(treecontab != NULL);
   assert(arity != NULL);
   assert(-1 < frompos);
   assert(frompos <= topos);
   assert(topos < (int) nvariables);

   if( frompos == topos )
      return;

   these_arities = (ARITY *) malloc( (nvariables-1) * sizeof(ARITY));
   for( i = 1; i < (int) nvariables; i++ )
      these_arities[i-1] = arity[variables[i]];

   pushvar_treecontab_reccall(treecontab, nvariables, arity[variables[0]], these_arities, frompos, topos);

   free(these_arities);

   tmpvar = variables[frompos];
   for( i = frompos; i < topos; i++ )
      variables[i] = variables[i+1];
   variables[topos] = tmpvar;
}
#endif

/* /\** In a tree-shaped contingency table, reorder the variables */
/*  * (and thus levels of the tree) by having a given variable swap places */
/*  * with the one coming after it in the order */
/*  *\/ */
/* static void swapvars_treecontab_reccall( */
/*    TREECONTAB* treecontab,     /\**< (pointer to) tree-shaped contingency table *\/ */
/*    const VARIABLE* variables,  /\**< variables for the contingency table *\/ */
/*    COUNT nvariables,           /\**< number of variables in the contingency table *\/ */
/*    const ARITY* arity,         /\**< arity[i] is arity of variable i *\/ */
/*    VARIABLE swapvar            /\**< variable whose position should be incremented *\/ */
/*    ) */
/* { */

/*    VALUE val0; */
/*    int arity0; */
/*    int arity1; */
/*    int val1; */
/*    TREECONTAB* above; */

/*    /\* just do necessary checks via assert statements *\/ */
/*    assert(treecontab != NULL); */
/*    /\* top-level call should never have treecontab->children == NULL */
/*       and neither should any recursive call due to (*GUARD*) below *\/ */
/*    assert(treecontab->children != NULL); */
/*    assert(variables != NULL); */
/*    assert(arity != NULL); */
/*    assert(nvariables > 1); */
/*    assert(swapvar < nvariables - 1); */

/*    arity0 = arity[variables[0]]; */

/*    if( swapvar != 0 ) */
/*    { */
/*       for(val0 = 0; val0 < arity0; val0++ ) */
/*       { */
/*          if( treecontab->children+val0 != NULL ) /\* (*GUARD*) *\/ */
/*             swapvars_treecontab_reccall(treecontab->children+val0,variables+1,nvariables-1,arity,swapvar-1); */
/*       } */
/*    } */
/*    else */
/*    { */
/*       SCIP_Bool counts; */

/*       counts = (nvariables == 2); */
/*       arity1 = arity[variables[1]]; */
/*       above = (TREECONTAB *) malloc(arity1 * sizeof(TREECONTAB)); */
/*       for( val1 = 0; val1 < arity1; val1++ ) */
/*       { */
/*          SCIP_Bool all_null_or_zero = TRUE; */
/*          above[val1].children = (TREECONTAB *) malloc(arity0 * sizeof(TREECONTAB)); */

/*          for( val0 = 0; val0 < arity0; val0++ ) */
/*          { */
/*             if( treecontab->children[val0].children != NULL ) */
/*             { */
/*                TREECONTAB tmp; */
/*                tmp = treecontab->children[val0].children[val1]; */
/*                /\* no need to check "counts" here *\/ */
/*                above[val1].children[val0] = tmp; */

/*                if( all_null_or_zero && ((!counts && tmp.children != NULL) || (counts && tmp.count != 0))) */
/*                   all_null_or_zero = FALSE; */
/*             } */
/*             else */
/*             { */
/*                if( counts ) */
/*                   above[val1].children[val0].count = 0; */
/*                else */
/*                   above[val1].children[val0].children = NULL; */
/*             } */
/*          } */

/*          if( all_null_or_zero ) */
/*          { */
/*             /\* sub-tree at val1 is all zeroes */
/*                so just set to NULL *\/ */
/*             free(above[val1].children); */
/*             above[val1].children = NULL; */
/*          } */
/*       } */
/*       /\* get rid of old branches *\/ */
/*       for( val0 = 0; val0 < arity0; val0++ ) */
/*       { */
/*          if( treecontab->children[val0].children != NULL ) */
/*             free(treecontab->children[val0].children); */
/*       } */
/*       free(treecontab->children); */

/*       /\* and replace with new sub-tree *\/ */
/*       treecontab->children = above; */
/*    } */
/* } */


/* /\** In a tree-shaped contingency table, reorder the variables */
/*  * (and thus levels of the tree) by having a given variable swap places */
/*  * with the one coming after it in the order */
/*  *\/ */
/* static void swapvars_treecontab( */
/*    TREECONTAB* treecontab,     /\**< (pointer to) tree-shaped contingency table *\/ */
/*    VARIABLE* variables,        /\**< variables for the contingency table *\/ */
/*    COUNT nvariables,           /\**< number of variables in the contingency table *\/ */
/*    const ARITY* arity,         /\**< arity[i] is arity of variable i *\/ */
/*    VARIABLE swapvar            /\**< variable whose position should be incremented *\/ */
/*    ) */
/* { */

/*    VARIABLE tmpvar; */

/*    /\* just do necessary checks via assert statements *\/ */
/*    assert(treecontab != NULL); */
/*    /\* top-level call should never have treecontab->children == NULL */
/*       and neither should any recursive call due to (*GUARD*) below *\/ */
/*    assert(treecontab->children != NULL); */
/*    assert(variables != NULL); */
/*    assert(arity != NULL); */
/*    assert(nvariables > 1); */
/*    assert(swapvar < nvariables - 1); */

/*    /\* re-arrange the tree *\/ */

/*    swapvars_treecontab_reccall(treecontab, variables, nvariables, arity, swapvar); */

/*    /\* just update the variables array to reflect new order *\/ */

/*    tmpvar = variables[swapvar]; */
/*    variables[swapvar] = variables[swapvar+1]; */
/*    variables[swapvar+1] = tmpvar; */
/* } */

#ifdef SCIP_DEBUG
/** print a tree-shaped contingency table (only for debugging) */
static void print_treecontab(
   const TREECONTAB treecontab,  /**< tree-shaped contingency table */
   const VARIABLE *variables,    /**< variables for the contingency table */
   COUNT nvariables,             /**< number of variables in the contingency table */
   const ARITY* arity            /**< arity[i] is arity of variable i */
   )
{
   VALUE val;

   assert(nvariables == 0 || variables != NULL);
   assert(arity != NULL);

   if( nvariables == 0 )
   {
      printf("count = %d\n", treecontab.count);
   }
   else
   {
      if( treecontab.children == NULL )
      {
         printf("NULL\n");
      }
      else
      {
         for(val = 0; val < arity[variables[0]]; val++ )
         {
            printf("X%d=%d,",variables[0],val);
            print_treecontab(treecontab.children[val],variables+1,nvariables-1,arity);
         }
      }
   }
}
#endif





/** produce a flat contingency table from a tree-shaped one */
static
void flatten_treecontab(
   const TREECONTAB treecontab,  /**< tree-shaped contingency table (N.B. not a pointer to one ) */
   const VARIABLE *variables,    /**< variables in the contingency table */
   COUNT nvariables,             /**< number of variables in the contingency table */
   const ARITY* arity,           /**< arity[i] is arity of variable i */
   COUNT** counts                /**< pointer to returned flat contingency table */
   )
{
   VALUE val;
   VALUE val2;
   COUNT j;
   ARITY i_arity;

   assert(nvariables == 0 || variables != NULL);
   assert(arity != NULL);
   assert(counts != NULL);
   assert(*counts != NULL);

   if(  nvariables == 0  )
   {
      **counts = treecontab.count;
      (*counts)++;
   }
   else
   {
      if( treecontab.children == NULL )
      {
         i_arity = 1;
         for( j = 0; j < nvariables; ++j )
            i_arity *= arity[variables[j]];

         for( val2 = 0; val2 < i_arity; ++val2 )
         {
            /* just fill up with zeroes */
            **counts = 0;
            (*counts)++;
         }
      }
      else
      {
         for(  val = 0; val < arity[variables[0]]; ++val  )
            flatten_treecontab((treecontab.children)[val],variables+1,nvariables-1,arity,counts);
      }
   }
}


/** Subtract one contingency table from another. The two contingency tables must be over the same variables
 */
static
void subtract_contab(
   TREECONTAB *treecontab1,        /**< Contingency table from which treecontab2 will be subtracted */
   const TREECONTAB *treecontab2,  /**< Contingency table which will be subtracted from treecontab1 */
   const VARIABLE *variables,      /**< Variables common to both treecontab1 and treecontab2 */
   int nvariables,                 /**< Number of variables */
   const ARITY *arity              /**< arity[i] is the arity of variable i, */
)
{

   VALUE val;

   TREECONTAB *pt1;
   TREECONTAB *pt2;
   ARITY firstarity;

   assert(treecontab1 != NULL);
   assert(treecontab2 != NULL);

   if( nvariables == 0 )
   {
      /* *treecontab1 and *treecontab2 must both be counts */
      assert(variables == NULL);
      treecontab1->count -= treecontab2->count;
      return;
   }

   /* *treecontab1 and *treecontab2 must both be TREECONTAB*
    * arrays of the same length */

   assert(variables != NULL);
   assert(nvariables > 0);
   firstarity = arity[variables[0]];
   nvariables--;

   pt1 = treecontab1->children;
   pt2 = treecontab2->children;
   if( nvariables == 0 )
      for( val = 0; val < firstarity; ++val )
         (pt1++)->count -= (pt2++)->count;
   else
   {
      variables++;
      for( val = 0; val < firstarity; ++val )
      {
         if( pt1->children != NULL && pt2->children != NULL )
            subtract_contab(pt1, pt2, variables, nvariables, arity);
         pt1++;
         pt2++;
      }
   }
   return;
}

/** Compute the (marginal) log-likelihood for a subset of the variables (modulo a constant)
 *
 * For a flat contingency table
 *
 * To get the true log-likelihood one would need to add lgamma(alpha) - log(alpha+N)
 * to the result of this function, where N is the number of samples. Since local scores
 * are always computed as this_function(family)-this_function(parents), there is no need to compute
 * the constant term
 *
 * @return the (marginal) log-likelihood for the given subset
 */
static SCORE log_likelihood_flat(
   const FLATCONTAB flatcontab,  /**< contingency table for variables */
   int flatcontabsize,           /**< contingency table size */
   const SCOREARG aijk,          /**< effective samples size (ESS) divided by the number of
                                    joint instantiations of the subset of variables */
   const SCORE lgaijk,           /**< log(Gamma(aijk)) */
   COUNT *npos_cells_ptr,        /**< pointer to number of cells in flatcontab with a positive count */
   SCIP_Bool fast                /**< If true then C's lgamma function is used for computing scores,
                                    otherwise a slower, more accurate sum of logs */
)
{

   int i;
   COUNT j;
   COUNT count;
   SCORE skore = 0;

   if( fast )
   {
      for(i = 0; i < flatcontabsize; i++)
      {
         count = flatcontab[i];
         if( count > 0 )
         {
            (*npos_cells_ptr)++;
            skore += lgamma(count + aijk) - lgaijk;
         }
      }
   }
   else
   {
      for(i = 0; i < flatcontabsize; i++)
      {
         count = flatcontab[i];
         if( count > 0 )
         {
            (*npos_cells_ptr)++;
            for( j = 0; j < count; ++j )
                  skore += log(j + aijk);
         }
      }
   }
      return skore;
}

/** Compute the (marginal) log-likelihood for a subset of the variables (modulo a constant)
 *
 * For a tree-shaped contingency table
 *
 * To get the true log-likelihood one would need to add lgamma(alpha) - log(alpha+N)
 * to the result of this function, where N is the number of samples. Since local scores
 * are always computed as this_function(family)-this_function(parents), there is no need to compute
 * the constant term
 * @return the (marginal) log-likelihood for the given subset
 */
static SCORE log_likelihood(
   const TREECONTAB *treecontab,  /**< contingency table for variables */
   const VARIABLE *variables,     /**< the subset of the variables */
   int nvariables,                /**< the size of the subset of the variables */
   const SCOREARG aijk,           /**< effective samples size (ESS) divided by the number of
                                     joint instantiations of the subset of variables */
   const SCORE lgaijk,            /**< log(Gamma(aijk)) */
   COUNT *npos_cells_ptr,         /**< pointer to number of cells in treecontab with a positive count */
   const ARITY *arity,            /**< arity[i] is the arity of variable i, */
   SCIP_Bool fast                 /**< If true then C's lgamma function is used for computing scores,
                                     otherwise a slower, more accurate sum of logs */
)
{

   VARIABLE firstvar;
   VALUE val;
   SCORE skore;
   TREECONTAB *adval;
   COUNT i;

   if( nvariables == 0 )
   {
      assert(variables == NULL);
      if( treecontab->count > 0 )
      {
         (*npos_cells_ptr)++;
         if( fast )
            return lgamma(treecontab->count + aijk) - lgaijk;
         else
         {
            skore = 0;
            for( i = 0; i < treecontab->count; ++i )
               skore += log(i + aijk);
            return skore;
         }
      }
      else
         return 0;
   }

   assert(variables != NULL);
   assert(nvariables > 0);
   firstvar = variables[0];
   nvariables--;
   skore = 0;
   adval = treecontab->children;
   if( nvariables == 0 )
      /* avoiding a recursive call here provides a surprisingly big speed-up */
      for( val = 0; val < arity[firstvar]; ++val )
      {
         if( adval->count > 0 )
         {
            (*npos_cells_ptr)++;
            if( fast )
               skore += (lgamma(adval->count + aijk) - lgaijk);
            else
               for( i = 0; i < adval->count; ++i )
                  skore += log(i + aijk);
         }
         adval++;
         /*skore += score_contab2(adval++,NULL,0,aijk,lgaijk);*/
      }
   else
   {
      variables++;
      for( val = 0; val < arity[firstvar]; ++val )
      {
         if( adval->children != NULL )
            skore += log_likelihood(adval, variables, nvariables, aijk, lgaijk, npos_cells_ptr, arity, fast);
         adval++;
      }
   }
   return skore;
}




/** Delete a tree-shaped contingency table.
 *  The contingency table must have at least one variable (i.e.\ it cannot be a count), but this only checked for in DEBUG mode.
 */
static
void delete_contab(
   TREECONTAB treecontab,       /**< Contingency table to delete */
   const VARIABLE *variables,   /**< Variables in the contingency table */
   int nvariables,              /**< Number of variables in the contingency table (must be positive) */
   const ARITY *arity           /**< arity[i] is the arity of variable i, */
)
{
   VALUE val;
   VALUE val2;
   const VARIABLE firstvar = variables[0];
   VARIABLE secondvar;

   assert(variables != NULL);
   assert(nvariables > 0);

   if( nvariables > 1 )
   {
      if( nvariables == 2 )
      {
         for( val = 0; val < arity[firstvar]; ++val )
            if( (treecontab.children[val]).children != NULL )
               free((treecontab.children)[val].children);
      }
      else if( nvariables == 3 )
      {
         secondvar = variables[1];
         for( val = 0; val < arity[firstvar]; ++val )
            if( (treecontab.children[val]).children != NULL )
            {
               for( val2 = 0; val2 < arity[secondvar]; ++val2 )
                  if( ((treecontab.children[val]).children)[val2].children != NULL )
                     free(((treecontab.children)[val].children)[val2].children);
               free((treecontab.children)[val].children);
            }
      }
      else
      {
         variables++;
         nvariables--;
         assert(variables != NULL);
         assert(nvariables > 0);
         for( val = 0; val < arity[firstvar]; ++val )
            if( (treecontab.children[val]).children != NULL )
               delete_contab((treecontab.children)[val], variables, nvariables, arity);
      }
   }
   free(treecontab.children);
   return;
}

/** Construct a flat contingency table from a leaflist
 *  ( contingency table must be for at least one variable )
 */
static
void makeflatcontableaf(
   const ROW *leaflist,        /**< datapoints for this query  */
   COUNT count,                /**< number of datapoints (in leaflist) */
   const VARIABLE *variables,  /**< variables in the contingency table (sorted) */
   const int *strides,         /**< stride sizes for each variable */
   int nvariables,             /**< number of variables in the contingency table */
   FLATCONTAB flatcontab,      /**< Contingency table initialised to zero */
   VALUE **data                /**< data[i][j] is the value of variable i in row j */
)
{
   COUNT j;
   int i;
   int k;
   ROW row;

   int stride0;
   const VALUE* data0;
   int stride1;
   const VALUE* data1;
   int stride2;
   const VALUE* data2;
   int stride3;
   const VALUE* data3;
   int stride4;
   const VALUE* data4;

   switch( nvariables )
   {
   case 1:
   {
       stride0 = strides[0];
       data0 = data[variables[0]];

      for( j = 0; j < count; ++j )
         flatcontab[stride0*data0[leaflist[j]]]++;

      break;
   }
   case 2:
   {
        stride0 = strides[0];
        data0 = data[variables[0]];
        stride1 = strides[1];
        data1 = data[variables[1]];


        for( j = 0; j < count; ++j )
        {
           row = leaflist[j];
           flatcontab[stride0*data0[row] + stride1*data1[row]]++;
        }

        break;
   }
   case 3:
   {
        stride0 = strides[0];
        data0 = data[variables[0]];
        stride1 = strides[1];
        data1 = data[variables[1]];
        stride2 = strides[2];
        data2 = data[variables[2]];


        for( j = 0; j < count; ++j )
        {
           row = leaflist[j];
           flatcontab[stride0*data0[row] + stride1*data1[row] + stride2*data2[row]]++;
        }

        break;
   }
   case 4:
   {
        stride0 = strides[0];
        data0 = data[variables[0]];
        stride1 = strides[1];
        data1 = data[variables[1]];
        stride2 = strides[2];
        data2 = data[variables[2]];
        stride3 = strides[3];
        data3 = data[variables[3]];


        for( j = 0; j < count; ++j )
        {
           row = leaflist[j];
           flatcontab[stride0*data0[row] + stride1*data1[row] + stride2*data2[row] + stride3*data3[row]]++;
        }

        break;
   }
   case 5:
   {
        stride0 = strides[0];
        data0 = data[variables[0]];
        stride1 = strides[1];
        data1 = data[variables[1]];
        stride2 = strides[2];
        data2 = data[variables[2]];
        stride3 = strides[3];
        data3 = data[variables[3]];
        stride4 = strides[4];
        data4 = data[variables[4]];


        for( j = 0; j < count; ++j )
        {
           row = leaflist[j];
           flatcontab[stride0*data0[row] + stride1*data1[row] + stride2*data2[row] + stride3*data3[row] + stride4*data4[row]]++;
        }

        break;
   }
   default:
   {
      for( j = 0; j < count; ++j )
      {
         row = leaflist[j];
         i = 0;
         for(k = 0; k < nvariables; k++)
            i += strides[k]*data[variables[k]][row];
         flatcontab[i]++;
      }

      break;
   }
   }
}

#if 0
/** Construct an indexed contingency table for the entire dataset
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE makeindexedcontab(
   SCIP* scip,                    /**< SCIP instance */
   int nvars,                     /**< number of variables in the data */
   const POSTERIOR* data,         /**< the data */
   INDEXEDCONTAB* indexcontab     /**< the indexed contingency table */
   )
{
   int row;
   COUNT* rows;
   COUNT* thisrow;
   COUNT* last;
   SCIP_Bool same;
   int nrows;
   const int rowlength = nvars + 1; /* last place is for the count */
   int v;

   /* make big enough to allow for no duplicates */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rows, (data->nrows)*rowlength) );
   nrows = 0;

   last = rows;
   for( row = 0; row < data->nrows; row++)
   {
      for( thisrow = rows; thisrow < last; thisrow += rowlength)
      {
         /* does data->data[row] have same index as thisrow? */
         same = TRUE;
         for( v = 0; v < nvars; v++ )
         {
            if( data->data[v][row] != thisrow[v] )
            {
               same = FALSE;
               break;
            }
         }

         if( same )
         {
            thisrow[nvars]++;
            break;
         }
      }

      if( !same )
      {
         for( v = 0; v < nvars; v++ )
         {
            thisrow[v] = data->data[v][row];
         }
         thisrow[nvars] = 1;
         last += rowlength;
         nrows++;
      }
   }
   /* shrink, hopefully */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &rows, (data->nrows)*rowlength, nrows*rowlength) );
   indexcontab->nrows = nrows;
   indexcontab->rows = rows;
   indexcontab->nvariables = nvars;

   return SCIP_OKAY;
}
#endif

#if 0
/** Construct a flat contingency table from an index contingency table
 */
static
void makeflatcontabfromidxcontab(
   const INDEXEDCONTAB *indexcontab,         /**< (Pointer to) the indexed source contingency table */
   const VARIABLE *variables,                /**< Variables in the sought contingency table (sorted) */
   const int *strides,                       /**< stride sizes for each variable */
   int nvariables,                           /**< Number of variables in the contingency table */
   FLATCONTAB flatcontab                     /**< Contingency table initialised with zeroes */
)
{
   const int countidx = indexcontab->nvariables;
   const int rowlength = indexcontab->nvariables + 1;
   int idx;
   COUNT* thisrow;
   COUNT* last;
   int i;

   last = indexcontab->rows + rowlength*indexcontab->nrows;

   /* just move the pointer along */
   for( thisrow = indexcontab->rows; thisrow < last; thisrow += rowlength)
   {
      idx = 0;
      for( i = 0; i < nvariables; i++ )
         idx += thisrow[variables[i]]*strides[i];
      flatcontab[idx] += thisrow[countidx];
   }
}
#endif

/** Construct a flat contingency table from an adtree
 */
static
void makeflatcontab(
   const ADTREE *adtree,       /**< (Pointer to) the ADTREE */
   VARIABLE offset,            /**< Offset for first variable (to identify correct vary nodes) */
   const VARIABLE *variables,  /**< Variables in the sought contingency table (sorted) */
   const int *strides,         /**< stride sizes for each variable */
   int nvariables,             /**< Number of variables in the contingency table */
   FLATCONTAB flatcontab,      /**< Contingency table initialised with zeroes */
   VALUE **data,               /**< data[i][j] is the value of variable i in row j */
   const ARITY *arity          /**< arity[i] is the arity of variable i, */
)
{

   VARYNODE vn;
   FLATCONTAB flatcontabmcv;
   VALUE val;
   FLATCONTAB flatcontabval;
   int i;

   assert(adtree != NULL);
   assert(flatcontab != NULL);
   assert(data != NULL);
   assert(arity != NULL);
   assert(nvariables == 0 || variables != NULL);
   assert(nvariables == 0 || strides != NULL);
   assert(adtree->children != NULL || adtree->leaflist != NULL );
   assert(adtree->children == NULL || adtree->leaflist == NULL );

   if( nvariables == 0 )
   {
      flatcontab[0] = adtree->count;
      return;
   }

   if( adtree->leaflist != NULL )
   {
      /* construct contingency table directly from data in leaf list */
      makeflatcontableaf(adtree->leaflist, adtree->count, variables, strides, nvariables, flatcontab, data);
      return;
   }

   /* find varynode for firstvar */
   vn = adtree->children[variables[0] - offset];

   flatcontabmcv = flatcontab+strides[0]*vn.mcv;

   /* make contingency table where variables[0] is marginalised away (and store at flatcontab + stride*vn.mcv) */
   makeflatcontab(adtree, offset, variables+1, strides+1, nvariables-1, flatcontabmcv, data, arity);

   for( val = 0; val < arity[variables[0]]; ++val )
   {
      /* if vn.children[val] == NULL then either val=vn.mcv or the contingency table for val is all zeroes
         so do nothing */
      if( vn.children[val] != NULL)
      {
         flatcontabval = flatcontab + strides[0]*val;
         makeflatcontab(vn.children[val], variables[0]+1, variables+1, strides+1, nvariables-1, flatcontabval, data, arity);
         /* subtract contingency table for val from that for mcv */
         for(i = 0; i < strides[0]; ++i)
         {
            flatcontabmcv[i] -= flatcontabval[i];
         }
      }
   }
}

/** Construct a tree-shaped contingency table from an adtree
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE maketreecontab(
   SCIP* scip,                 /**< SCIP instance */
   const ADTREE *adtree,       /**< (Pointer to) the ADTREE or NULL */
   VARIABLE offset,            /**< Offset for first variable (to identify correct vary nodes) */
   const VARIABLE *variables,  /**< Variables in the sought contingency table (sorted) */
   int nvariables,             /**< Number of variables in the contingency table */
   TREECONTAB *treecontab,         /**< Returned contingency table */
   VALUE **data,               /**< data[i][j] is the value of variable i in row j */
   const ARITY *arity          /**< arity[i] is the arity of variable i, */
)
{

   VARIABLE firstvar;
   ARITY firstarity;

   VALUE val;
   VALUE mcv;

   VARYNODE vn;
   ADTREE **vnchildren;

   TREECONTAB *ptmcv;
   TREECONTAB *ptval;


   assert(adtree != NULL);
   assert(treecontab != NULL);
   assert(adtree->leaflist == NULL || adtree->children == NULL);
   assert(nvariables > 0 || variables == NULL);

   if( nvariables == 0 )
   {
      treecontab->count = adtree->count;
      return SCIP_OKAY;
   }

   assert(nvariables > 0);
   assert(variables != NULL);

   if( adtree->leaflist != NULL )
   {
      assert(adtree->children == NULL);
      /* construct contingency table directly from data in leaf list */
      SCIP_CALL( maketreecontableaf(scip, adtree->leaflist, adtree->count, variables, nvariables, treecontab, data, arity) );
      return SCIP_OKAY;
   }

   assert(adtree->children != NULL);
   firstvar = variables[0];
   firstarity = arity[firstvar];

   /* need to create a contab for each value of firstvar */

   /* find varynode for firstvar */
   vn = adtree->children[firstvar - offset];
   mcv = vn.mcv;

   assert(vn.children != NULL);
   assert(mcv < firstarity);

   nvariables--;
   if( nvariables == 0 )
      variables = NULL;
   else
      variables++;

   treecontab->children = (TREECONTAB *) malloc(firstarity * sizeof(TREECONTAB));

   maketreecontab(scip, adtree, offset, variables, nvariables, (treecontab->children) + mcv, data, arity);

   ptmcv = (treecontab->children) + mcv;
   ptval = treecontab->children;
   offset = firstvar + 1;
   vnchildren = vn.children;
   for( val = 0; val < firstarity; ++val )
   {
      if( val != mcv )
      {
         if( *vnchildren == NULL )
            (*ptval).children = NULL;
         else
         {
            maketreecontab(scip, *vnchildren, offset, variables, nvariables, ptval, data, arity);
            subtract_contab(ptmcv, ptval, variables, nvariables, arity);
         }
      }
      ptval++;
      vnchildren++;
   }
   return SCIP_OKAY;
}


/** Compute the (marginal) log-likelihood for a subset of the variables and store in cache, or retrieve from cache if already computed
 * @return the (marginal) log-likelihood for the given subset
 */
static SCORE log_likelihood_cache(
   SCIP* scip,                  /**< SCIP instance */
   const ADTREE* adtree,        /**< the data */
   const VARIABLE* variables,   /**< the subset of the variables */
   int nvariables,              /**< the size of the subset of the variables */
   COUNT* npos_cells_ptr,       /**< pointer to number of cells in treecontab with a positive count */
   SCORECACHE* scorecache,      /**< score cache */
   const double alpha,          /**< The 'effective sample size' governing the BDeu score */
   const ARITY* arity,          /**< arity[i] is the arity of variable i, */
   SCIP_Bool fast,              /**< If true then C's lgamma function is used for computing scores,
                                   otherwise a slower, more accurate sum of logs */
   VALUE** data,                /**< data[i][j] is the value of variable i in row j */
   int maxflatcontabsize        /**< maximum allowed size for a flat contingency table (if exceeded a tree-shaped contingency
                                   table is used) */
)
{

   SCORE llh;
   TREECONTAB treecontab;
   FLATCONTAB flatcontab;
   double aijk;
   double lgaijk;
   int i;
   int v;
   int flatcontabsize;
   int* strides;

   int rank;

   llh = get_score_count_rank(scorecache, variables, nvariables, &rank, npos_cells_ptr);

   /* if log-likelihood in cache, just return it ... */
   if( *npos_cells_ptr != 0 )
      return llh;

   /* ... else have to compute log-likelihood */
   aijk = alpha;
   for( v = 0; v < nvariables; ++v )
      aijk /= arity[variables[v]];
   lgaijk = lgamma(aijk);
   (*npos_cells_ptr) = 0;

   /* compute size for a flat contingency table (and associated 'strides') */
   SCIP_CALL( SCIPallocBufferArray(scip, &strides, nvariables) );
   flatcontabsize = 1;
   for( i = nvariables - 1; i >= 0; i--)
   {
      strides[i] = flatcontabsize;
      flatcontabsize *= arity[variables[i]];
   }

   if( flatcontabsize > maxflatcontabsize )
   {
      maketreecontab(scip, adtree, 0, variables, nvariables, &treecontab, data, arity);
      llh = log_likelihood(&treecontab, variables, nvariables, aijk, lgaijk, npos_cells_ptr, arity, fast);
      if( nvariables > 0 )
         delete_contab(treecontab, variables, nvariables, arity);
   }
   else
   {
      SCIP_CALL( SCIPallocClearBufferArray(scip, &flatcontab, flatcontabsize) );
      makeflatcontab(adtree, 0, variables, strides, nvariables, flatcontab, data, arity);
      llh = log_likelihood_flat(flatcontab, flatcontabsize, aijk, lgaijk, npos_cells_ptr, fast);
      SCIPfreeBufferArray(scip, &flatcontab);
   }
   SCIPfreeBufferArray(scip, &strides);

   assert(*npos_cells_ptr > 0);

   /* if possible store log-likelihood and count in cache, using previously computed rank
      if rank == -1, then no room in the cache, so not possible
   */
   if( rank != -1)
      SCIP_CALL( set_score_count_from_rank(scip, scorecache, rank, llh, *npos_cells_ptr) );

   return llh;
}

/** Add scoring parameters
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
SCIP_RETCODE SC_addScoringParameters(
   SCIP* scip  /**< SCIP data structure */
)
{
   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/scoring/prune",
                             "whether to prune during scoring",
                             TRUE
                            ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/rmin",
                            "minimum number of datapoints to create a branch in the AD tree",
                            2000, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/adtreedepthlim",
                            "limit on the depth of the AD tree",
                            1000, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/adtreenodeslim",
                            "limit on the number of nodes in  the AD tree",
                            10000, 0, INT_MAX
                           ));


   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/cachesizelimit",
                            "limit on number of scores that are cached",
                            10000000, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/cacheblocksize",
                            "how much to increase the cache when needed and allowed",
                            10000, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/nvarscachelimit",
                            "subsets must have size below this to be cached",
                            50, 0, INT_MAX
                           ));

   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/maxflatcontabsize",
                            "maximum size for a flat contingency table (as opposed to ADTREE tree contingency table) ",
                            10000, 0, INT_MAX
                           ));


   SCIP_CALL(UT_addIntParam(scip,
                            "gobnilp/scoring/palim",
                            "maximum number of parents for each node (-1 for no limit)",
                            3, -1, INT_MAX
                           ));

   SCIP_CALL(UT_addRealParam(scip,
                             "gobnilp/scoring/alpha",
                             "alpha value for use in BDeu scoring",
                             1, 0, SCIP_REAL_MAX
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/scoring/names",
                             "whether variable names are given in the data file",
                             TRUE
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/scoring/arities",
                             "whether variable arities are given in the data file",
                             TRUE
                            ));

   SCIP_CALL(UT_addRealParam(scip,
                             "gobnilp/scoring/prunegap",
                             "If Pa1 and Pa2 are parent sets for some variable, Pa1 is a subset of Pa2 and local_score(Pa1) >= local_score(Pa2) - prunegap then Pa2 will be pruned",
                             0, SCIP_REAL_MIN, SCIP_REAL_MAX
                            ));

   SCIP_CALL(UT_addBoolParam(scip,
                             "gobnilp/scoring/fast",
                             "whether to fast scoring (which is inaccurate for high values of ESS)",
                             TRUE
                            ));

   SCIP_CALL(UT_addStringParam(scip,
         "gobnilp/scoring/score_type",
         " Set to either \"BDeu\" (default) or \"BGe\" or \"BIC\" ",
         "BDeu"
         ));


   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/scoring/continuous",
         "whether the data is continuous (as opposed to entirely discrete)",
         FALSE
         ));

   SCIP_CALL(UT_addIntParam(scip,
         "gobnilp/scoring/alpha_omega_minus_nvars",
         "Added to the number of variables to get the degrees of freedom (alpha_omega) for Wishart prior for BGe scoring",
         2, 2, INT_MAX
         ));


   SCIP_CALL(UT_addRealParam(scip,
         "gobnilp/scoring/alpha_mu",
         "alpha_mu scaling parameter in normal-Wishart prior distribution for BGe scoring",
         1.0, 0.000001, SCIP_REAL_MAX
         ));


   SCIP_CALL(UT_addStringParam(scip,
         "gobnilp/scoring/nu",
         "If non-empty a comma separated sequence of values defining the prior mean  for BGe scoring. If empty (the default) the sample mean is used ",
         ""
         ));

   SCIP_CALL(UT_addStringParam(scip,
         "gobnilp/scoring/T",
         "Currently this value is ignored",
         /* "If non-empty a comma and newline separated sequence of values defining the prior parametric matrix for BGe scoring. If empty (the default) then this matrix is tI as suggested in Section B.2 of Kuipers et al ", */
         ""
         ));


   SCIP_CALL(UT_addBoolParam(scip,
         "gobnilp/scoring/pricing",
         "whether pricing is enable",
         FALSE
         ));

   SCIP_CALL(UT_addIntParam(scip,
         "gobnilp/scoring/initpalim",
         "maximum number of parents for each node in initial scoring (before any pricing ) (-1 for no limit)",
         3, -1, INT_MAX
         ));


   return SCIP_OKAY;
}
/** Record local scores for parent sets of a particular child
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE addParentSets(
   SCIP* scip,                      /**< SCIP data structure */
   VARIABLE child,                  /**< Child variable whose scored parent sets are being added */
   TRIENODE** parent_sets,          /**< Trie nodes representing parent sets */
   int* parent_sets_sizes,          /**< Parent set sizes */
   int num_parent_sets,             /**< The number of scored parent sets for child */
   ParentSetData* psd,              /**< Parent set data to populate */
   SCIP_Real*** scores              /**< (*scores)[child][j] will be the score for the jth parent set of child */
)
{
   int i, j;

   psd->nParentSets[child] = num_parent_sets;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->nParents[child]),   num_parent_sets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->ParentSets[child]), num_parent_sets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*scores)[child]),       num_parent_sets) );

   for( i = 0; i < num_parent_sets; i++ )
   {
      TRIENODE* node = parent_sets[i];
      int npa = parent_sets_sizes[i];

      assert(npa > -1);

      assert(node != NULL);
      (*scores)[child][i] = node->score;
      psd->nParents[child][i] = npa;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->ParentSets[child][i]), npa) );
      for( j = npa-1; j >= 0; j-- )
      {
         assert(node != NULL);
         psd->ParentSets[child][i][j] = node->elt;
         node = node->mother;
      }
   }
   return SCIP_OKAY;
}

/** Find the integer encoding for (the string representation of) a value of a variable.
 *  If the input string @c name does not occurs in @c names, then it will be added to it
 *  (as long as doing so would not exceed the arity), and  <tt>(*used)</tt> will be incremented.
 *  @return the integer encoding the value, or -1 not possible due to exceeding the given variable arity @c length.
*/
static
int lookup(
   char** names,        /**< @c names[r] is (the string representation of) the rth value of the variable */
   ARITY length,        /**< arity of the variable (number of values the variable has) */
   ARITY* used,         /**< the number of strings stored in  @c names */
   char* name           /**< (string representation of) a value of a variable */
   )
{
   ARITY i;
   for( i = 0; i < (*used); i++ )
      if( strcmp(names[i], name) == 0 )
         /* Found it */
         return i;
   if( i < length )
   {
      /* Not found but space to create it */
      sprintf(names[i], "%s", name);
      (*used)++;
      return i;
   }
   else
   {
      /* Not found and no space for it */
      return -1;
   }
}

#ifdef BLAS
/** Read continuous data from a file.
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE readContinuousProblemFromFile(
      SCIP* scip,                                   /**< SCIP data structure */
      const char* filename,                         /**< File containing the data */
      int num_delims,                               /**< The number of field delimiters to use */
      char* delims,                                 /**< The field delimiters to use */
      SCIP_Bool merge_delims,                       /**< Whether multiple field delimiters should be merged in to one */
      char*** nodeNames,                            /**< (Pointer to) node (i.e.\ variable) names */
      Bge_Matrix** data,                            /**  (Pointer to) an array storing the data of the variables */
      int* num_vars,                                /**< (Pointer to) the number of variables */
      int* num_rows                                 /**< (Pointer to) the number of rows (i.e.\ datapoints) in the data */
)
{
  FILE* file;

  char*** lines;
  int num_lines;
  int* line_lengths;

  SCIP_Bool hasNames;

  int i, j, k;

  /* Read the data from the file */
  if( strcmp(filename, "-") == 0 )
    file = stdin;
  else
    file = fopen(filename, "r");
  if( file == NULL )
  {
    SCIPerrorMessage("Could not open file %s.\n", filename);
    return SCIP_NOFILE;
  }
  SCIP_CALL( UT_readFileAndSplit(scip, file, delims, num_delims, merge_delims,
             &lines, &num_lines, &line_lengths) );
  fclose(file);

  /* Check the number of lines is ok */
  SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/names", &hasNames) );

  (*num_rows) = 0;

  for( i = 0; i < num_lines; i++ )
  {
    if( lines[i][0][0] != '#' )
      (*num_rows)++;
  }

  if( hasNames )
    (*num_rows)--;

  if( (ROW)(*num_rows) > UINT_MAX )
  {
    SCIPerrorMessage("Warning: Too many rows to store them as unsigned ints.\n");
    return SCIP_READERROR;
  }

  /* Check the line lengths are all ok */
  (*num_vars) = -1;
  for( i = 0; i < num_lines; i++ )
  {
    if( lines[i][0][0] == '#' )
    {
      /* Comment line, so do nothing */
    }
    else if( (*num_vars) < 0 )
    {
      /* Found the first line which wasn't a comment */
      (*num_vars) = line_lengths[i];
    }
    else if( line_lengths[i] != (*num_vars) )
    {
      /* Line is wrong length */
      SCIPerrorMessage("Wrong number of data items on line %d.  Found %d when %d  were expected.\n",
                        i + 1, line_lengths[i], (*num_vars));
      return SCIP_READERROR;
    }
    else
    {
      /* Line is correct length.  Do nothing */
    }
  }

  /* Get names */
  SCIP_CALL( SCIPallocMemoryArray(scip, nodeNames, (*num_vars)) );

  for( i = 0; i < (*num_vars); i++ )
    SCIP_CALL( SCIPallocMemoryArray(scip, &((*nodeNames)[i]), SCIP_MAXSTRLEN) );

  if( hasNames )
  {
    /* First non-comment line will be names */
    int name_line = 0;

    while( lines[name_line][0][0] == '#' )
      name_line++;

    for( j = 0; j < (*num_vars); j++ )
      sprintf((*nodeNames)[j], "%s", lines[name_line][j]);
  }
  else
  {
    /* No names - so use column numbers */
    for( i = 0; i < (*num_vars); i++ )
      sprintf((*nodeNames)[i], "%d", i);
  }

  /* get values */
  *data = BgeMatrixCreate(*num_vars, *num_rows);

  k = 0;

  for(i = 0; i < *num_rows; i++)
  {
    if(hasNames)
    {
      /* if the variables have name skip the first line */
      k = i + 1;
    }

    if(lines[k][0][0] != '#') /* if the line isn't a comment */
    {
      for(j = 0; j < *num_vars; j++)
      {
        /* check if the string is a float format */
        if(lines[k][j][0] != '0' && atof(lines[k][j]) == 0.0F)
        {
          printf("Incorrect format for column entry: %i and row entry: %i\n", j + 1, k + 1);
          return SCIP_READERROR;
        }

        /* read data into the data matrix */
        (*data)->items[i * (*num_vars) + j] = atof(lines[k][j]);
      }
    }
  }

   /* Tidy up */
   for( i = 0; i < num_lines; i++ )
   {
      for( j = 0; j < line_lengths[i]; j++ )
         SCIPfreeMemoryArray(scip, &(lines[i][j]));
      SCIPfreeMemoryArray(scip, &(lines[i]));
   }
   SCIPfreeMemoryArray(scip, &lines);
   SCIPfreeMemoryArray(scip, &line_lengths);

  return SCIP_OKAY;

}
#endif

/** Read discrete data from a file.
 *  The last 7 arguments are outputs for which we use pointers.
 *  This function allocates space for these outputs.
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static SCIP_RETCODE readProblemFromFile(
   SCIP* scip,                           /**< SCIP data structure */
   const char* filename,                 /**< File containing the data */
   int num_delims,                       /**< The number of field delimiters to use */
   char* delims,                         /**< The field delimiters to use */
   SCIP_Bool merge_delims,               /**< Whether multiple field delimiters should be merged in to one */
   char*** nodeNames,                    /**< (Pointer to) node (i.e.\ variable) names */
   ARITY** arities,                      /**< (Pointer to) node variable arities */
   VALUE*** items,                       /**< (Pointer to) the data. &c (*items)[i][j] is the
                                            integer-encoded value of the ith variable in the jth datapoint */
   int* num_vars,                        /**< (Pointer to) the number of variables */
   int* num_rows,                        /**< (Pointer to) the number of rows (i.e.\ datapoints) in the data */
   char**** labels,                      /**< (Pointer to) a mapping from the integer encoding of a value to the corresponding string
                                              @c (*labels)[i][r] is the string corresponding to the rth value of the ith variable */
   int** num_observed                    /**< (Pointer to) the number of distinct values observed in the data for each variable */
)
{
   FILE* file;

   char*** lines;
   int num_lines;
   int* line_lengths;

   SCIP_Bool hasNames;
   SCIP_Bool hasArities;

   char*** label_map;
   ARITY* label_map_count;

   int i, j, k;

   ARITY i_arity;

   int tmp_arity;

   SCIP_Bool floatwarninggiven = FALSE;

   /* Read the data from the file */
   if( strcmp(filename, "-") == 0 )
      file = stdin;
   else
      file = fopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open file %s.\n", filename);
      return SCIP_NOFILE;
   }
   SCIP_CALL( UT_readFileAndSplit(scip, file, delims, num_delims, merge_delims, &lines, &num_lines, &line_lengths) );
   fclose(file);

   /* Check the number of lines is ok */
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/names", &hasNames) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/arities", &hasArities) );
   (*num_rows) = 0;
   for( i = 0; i < num_lines; i++ )
      if( lines[i][0][0] != '#' )
         (*num_rows)++;
   if( hasNames )
      (*num_rows)--;
   if( hasArities )
      (*num_rows)--;
   if( (ROW)(*num_rows) > UINT_MAX )
   {
      SCIPerrorMessage("Warning: Too many rows to store them as unsigned  ints.\n");
      return SCIP_READERROR;
   }
   /* Check the line lengths are all ok */
   (*num_vars) = -1;
   for( i = 0; i < num_lines; i++ )
   {
      if( lines[i][0][0] == '#' )
      {
         /* Comment line, so do nothing */
      }
      else if( (*num_vars) < 0 )
      {
         /* Found the first line which wasn't a comment */
         (*num_vars) = line_lengths[i];
      }
      else if( line_lengths[i] != (*num_vars) )
      {
         /* Line is wrong length */
         SCIPerrorMessage("Wrong number of data items on line %d.  Found %d when %d were expected.\n",
            i + 1, line_lengths[i], (*num_vars));
         return SCIP_READERROR;
      }
      else
      {
         /* Line is correct length.  Do nothing */
      }
   }

   /* Get names */
   SCIP_CALL( SCIPallocMemoryArray(scip, nodeNames, (*num_vars)) );
   for( i = 0; i < (*num_vars); i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*nodeNames)[i]), SCIP_MAXSTRLEN) );
   if( hasNames )
   {
      /* First non-comment line will be names */
      int name_line = 0;
      while( lines[name_line][0][0] == '#' )
         name_line++;
      for( j = 0; j < (*num_vars); j++ )
         sprintf((*nodeNames)[j], "%s", lines[name_line][j]);
   }
   else
   {
      /* No names - so use column numbers */
      for( i = 0; i < (*num_vars); i++ )
         sprintf((*nodeNames)[i], "%d", i);
   }

   /* Get arities */
   SCIP_CALL( SCIPallocMemoryArray(scip, arities, (*num_vars)) );
   if( hasArities )
   {
      /* First non-comment line will be arities, or second if names are being used */
      int arity_line = 0;
      if( hasNames )
      {
         while( lines[arity_line][0][0] == '#' )
            arity_line++;
         arity_line++;
      }
      while( lines[arity_line][0][0] == '#' )
         arity_line++;
      for( j = 0; j < (*num_vars); j++ )
      {

         if( !isarity(lines[arity_line][j]) )
         {
            SCIPerrorMessage("Line %d of %s: %s does not represent a valid arity for a discrete variable.\n",
               arity_line+1, filename, lines[arity_line][j]);
            return SCIP_READERROR;
         }

         /* use %d to skip whitespace */
         sscanf(lines[arity_line][j], "%d", &tmp_arity);
         (*arities)[j] = (ARITY) tmp_arity;
      }
   }
   else
   {
      /* No arities - so use biggest possible values for now */
      for( i = 0; i < (*num_vars); i++ )
         (*arities)[i] = MAXARITY;
   }

   /* Get values */
   SCIP_CALL( SCIPallocMemoryArray(scip, items, (*num_vars)) );
   for( i = 0; i < (*num_vars); i++ )
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*items)[i]), (*num_rows)) );

   SCIP_CALL( SCIPallocMemoryArray(scip, &label_map, (*num_vars)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &label_map_count, (*num_vars)) );
   for( i = 0; i < (*num_vars); i++ )
   {
      label_map_count[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &(label_map[i]), (*arities)[i]) );
      for( i_arity = 0; i_arity < (*arities)[i]; i_arity++ )
         SCIP_CALL( SCIPallocMemoryArray(scip, &(label_map[i][i_arity]), SCIP_MAXSTRLEN) );
   }

   i = 0;
   if( hasNames )
   {
      while( lines[i][0][0] == '#' )
         i++;
      i++;
   }
   if( hasArities )
   {
      while( lines[i][0][0] == '#' )
         i++;
      i++;
   }
   for( j = 0; j < (*num_rows); j++ )
   {
      while( lines[i][0][0] == '#' )
         i++;
      for( k = 0; k < (*num_vars); k++ )
      {
         int value;

         if( !floatwarninggiven && lookslikefloat(lines[i][k]) )
         {
            SCIPerrorMessage("Warning: %s on row %d looks like a float (but interpreted as discrete value).\n",
               lines[i][k], i);
            floatwarninggiven = TRUE;
         }

         value = lookup(label_map[k], (*arities)[k], &(label_map_count[k]), lines[i][k]);
         if( value == -1 )
         {
            SCIPerrorMessage("Invalid value %s for discrete variable %s on row %d (%d values have already been seen)\n",
               lines[i][k], (*nodeNames)[k], i, (*arities)[i]);
            return SCIP_READERROR;
         }
         (*items)[k][j] = value;
      }
      i++;
   }

   /* Work out arities if they weren't given */
   if( !hasArities )
   {
      for( i = 0; i < (*num_vars); i++ )
         (*arities)[i] = label_map_count[i];
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, labels, (*num_vars)) );
   for( i = 0; i < (*num_vars); i++ )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &((*labels)[i]), (*arities)[i]) );
      for( i_arity = 0; i_arity < (*arities)[i]; i_arity++ )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &((*labels)[i][i_arity]), SCIP_MAXSTRLEN) );
         if( i_arity < label_map_count[i] )
            sprintf((*labels)[i][i_arity], "%s", label_map[i][i_arity]);
         else
            sprintf((*labels)[i][i_arity], "UNKNOWN");
      }
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, num_observed, (*num_vars)) );
   for( i = 0; i < (*num_vars); i++ )
   {
      (*num_observed)[i] = label_map_count[i];
   }

   /* Tidy up */
   for( i = 0; i < num_lines; i++ )
   {
      for( j = 0; j < line_lengths[i]; j++ )
         SCIPfreeMemoryArray(scip, &(lines[i][j]));
      SCIPfreeMemoryArray(scip, &(lines[i]));
   }
   SCIPfreeMemoryArray(scip, &lines);
   SCIPfreeMemoryArray(scip, &line_lengths);
   for( i = 0; i < (*num_vars); i++ )
   {
      for( i_arity = 0; i_arity < (*arities)[i]; i_arity++ )
         SCIPfreeMemoryArray(scip, &(label_map[i][i_arity]));
      SCIPfreeMemoryArray(scip, &(label_map[i]));
   }
   SCIPfreeMemoryArray(scip, &label_map);
   SCIPfreeMemoryArray(scip, &label_map_count);

   return SCIP_OKAY;
}

/** Compute the (marginal) log-likelihood for a subset of the variables and store in cache, or retrieve from cache if already computed
 * @return the (marginal) log-likelihood for the given subset
 */
static
SCORE log_likelihood_cache_compact(
   const VARIABLE* family, /**< family */
   int nfamily,            /**< family size */
   const PRIOR* prior,     /**< prior information for scoring, e.g. hyperparameters */
   const POSTERIOR* data,  /**< the data */
   SCORING_EXTRAS* extras  /**< any necessary extra arguments required for scoring */
   )
{
   return log_likelihood_cache(extras->scip,data->adtree,family,nfamily,&(extras->npos_cells),extras->scorecache,
      prior->alpha,prior->arity,extras->fast,data->data,extras->maxflatcontabsize);
}

/**< make a 'family' by inserting a child into a parent set
 * (correct amount of space) for family must already have been allocateed
 */
static
void makefamily(
   VARIABLE* family,        /**< family to populate */
   const VARIABLE* parents, /**< parent set */
   int nparents,            /**< size of parent set */
   VARIABLE child           /**< child variable */
   )
{

   int i;

   /* make family by inserting child into parents */
   i = 0;
   while( i < nparents && parents[i] < child )
   {
      assert(parents[i] != child);
      family[i] = parents[i];
      i++;
   }
   family[i] = child;          /* insert child */
   while( i < nparents )
   {
      assert(parents[i] != child);
      family[i+1] = parents[i];
      i++;
   }
}


/** compute a (fitted) log-likelihood score from a flat contingency table for a family with given child
 * @return (fitted) log-likelihood score
*/
static
SCORE ll_score_flat(
   FLATCONTAB flatcontab,    /**< contingency table with relevant counts for the family */
   int flatcontabsize,       /**< size of contingency table */
   const VARIABLE* family,   /**< family */
   int nfamily,              /**< size of family */
   VARIABLE child,           /**< child of family */
   const ARITY *arity        /**< arity[i] is the arity of variable i, */
   )
{
   SCORE skore;
   int block;
   int miniblock;
   int blocksize;
   int nonzeroes;
   COUNT total;
   COUNT tmp;
   int i;
   int stepsize;
   double floattotal;

   /* stepsize is how far apart child values are for a given instantiation of the parents */
   stepsize = 1;
   for(i = nfamily-1; i >= 0; i--)
   {
      if( family[i] == child )
         break;
      else
         stepsize *= arity[family[i]];
   }

   /* blocksize is how far apart are insts of the subset of those parents with lower index
      than the child */
   blocksize = stepsize * arity[child];
   skore = 0;

   /* treat each block separately */
   for(block = 0; block < flatcontabsize; block += blocksize)
   {
      /* within each block iterate through insts of those parents with index higher than
         the child */
      for( miniblock = block; miniblock < block + stepsize; miniblock++)
      {
         total = 0;
         nonzeroes = 0;

         /* sum over child values for a particular inst of parents */
         for( i = miniblock; i < miniblock + blocksize; i+= stepsize )
         {
            tmp = flatcontab[i];
            if( tmp > 0 )
            {
               total += tmp;
               nonzeroes++;
            }
         }

         /* need a least two child values to be non-zero for Nijk*log(Nijk/Nij) to be non-zero */
         if( nonzeroes > 1 )
         {
            floattotal = (double) total;
            for( i = miniblock; i < miniblock + blocksize; i+= stepsize )
            {
               tmp = flatcontab[i];
               if( tmp > 0 )
                  skore += tmp * log(tmp / floattotal);
            }
         }
      }
   }
   return skore;
}

/** compute a (fitted) log-likelihood score from a tree-shaped contingency table for a family with given
 * @return (fitted) log-likelihood score
 */
static
SCORE ll_score(
   SCIP* scip,             /**< SCIP instance */
   TREECONTAB treecontab,  /**< contingency table with relevant counts */
   const VARIABLE* family, /**< family */
   int nfamily,            /**< size of family */
   VARIABLE child,         /**< child of family */
   const ARITY *arity      /**< arity[i] is the arity of variable i, */
   )
{
   /* not using VALUE val; to stop warnings */
   int val;
   SCORE skore;
   COUNT** counts;
   int nparentinsts;
   int i;
   COUNT nij;
   COUNT* tmp;

   assert(scip != NULL);
   assert(family != NULL);
   assert(nfamily > 0);
   assert(arity != NULL);

   /* since treecontab is the contingency table for family which is non-empty
      (nfamily > 0), treecontab is not a count, so use treecontab.children */


   /* all counts in the contingency table are zero, so is the log-likelihood score */
   if( treecontab.children == NULL )
      return 0;

   skore = 0;

   if( family[0] != child )
   {
      for( val = 0; val < arity[family[0]]; ++val )
         skore += ll_score(scip,(treecontab.children)[val], family+1, nfamily-1, child, arity);
   }
   else
   {
      assert(family[0] == child);

      /* find the number of instantiations of the parents for each value of child */
      nparentinsts = 1;
      for(i = 1; i < nfamily; i++)
         nparentinsts *= arity[family[i]];

      SCIP_CALL( SCIPallocBufferArray(scip, &counts, arity[child]) );
      for( val = 0; val < arity[child]; ++val )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(counts[val]), nparentinsts) );
         tmp = counts[val];
         flatten_treecontab((treecontab.children)[val],family+1,nfamily-1,arity,&tmp);
      }

      /* for each of the 'nparentinsts' instantiations of the parents ... */
      for( i = 0; i < nparentinsts; i++ )
      {
         nij = 0;
         for( val = 0; val < arity[child]; val++ )
            nij += counts[val][i];

         if( nij == 0 )
            continue;

         /* counts[val][i] = "n_ijk" */
         for( val = 0; val < arity[child]; val++ )
            if( counts[val][i] != 0 )
               skore += counts[val][i] * log( (float) counts[val][i] / nij );
      }

      for( val = arity[child]-1;  val > -1; --val )
         SCIPfreeBufferArray(scip, &(counts[val]));
      SCIPfreeBufferArray(scip, &counts);

   }
   return skore;
}

/** Compute the BIC score for given child and parent set
 * @return BIC score
 */
static
SCORE bic_score(
   VARIABLE child,          /**< the child of the family */
   const VARIABLE* parents, /**< the parents of the family */
   int nparents,            /**< the number of parents in the family */
   const PRIOR* prior,      /**< prior information for scoring, e.g. hyperparameters */
   const POSTERIOR* data,   /**< the data */
   SCORING_EXTRAS* extras,  /**< any necessary extra arguments required for scoring */
   const FORCHILD* forchild /**< useful quantities specific to the child */
   )
{
   VARIABLE* family;
   int nfamily = nparents + 1;
   TREECONTAB treecontab;
   int nparentinsts;
   int i;
   SCORE ll;

   int flatcontabsize;
   FLATCONTAB flatcontab;
   int* strides;

   /* make family by inserting child into parents */
   SCIP_CALL( SCIPallocBufferArray(extras->scip, &family, nfamily) );
   makefamily(family,parents,nparents,child);

   /* compute number of joint instantiations of the parents */
   nparentinsts = 1;
   for( i = 0; i < nparents; i++)
      nparentinsts *= (prior->arity)[parents[i]];

   /* compute size for a flat contingency table (and associated 'strides') */
   SCIP_CALL( SCIPallocBufferArray(extras->scip, &strides, nfamily) );
   flatcontabsize = 1;
   for( i = nfamily - 1; i >= 0; i--)
   {
      strides[i] = flatcontabsize;
      flatcontabsize *= prior->arity[family[i]];
   }


   if( flatcontabsize > extras->maxflatcontabsize )
   {
      /* create TREECONTAB tree contingency table for family (parents + child) */
      maketreecontab(extras->scip, data->adtree, 0, family, nfamily, &treecontab, data->data, prior->arity);

      /* need to check this is the correct treecontab */
#ifdef SCIP_DEBUG
      SCIPdebugMessage("This child: %d (%s)\n",child,prior->nodeNames[child]);
      SCIPdebugMessage("This family\n");
      for(i = 0; i < nfamily; i++)
         SCIPdebugMessage("\t %d (%s)\n",family[i],prior->nodeNames[family[i]]);
      SCIPdebugMessage("\nThis contingency table\n");
      print_treecontab(treecontab, family, nfamily,prior->arity);
#endif
      ll = ll_score(extras->scip, treecontab, family, nfamily, child, prior->arity);
      delete_contab(treecontab,family,nfamily,prior->arity);
   }
   else
   {
      /* using SCIPallocClearBlockMemoryArray initialises to zeroes */
      SCIP_CALL( SCIPallocClearBlockMemoryArray(extras->scip, &flatcontab, flatcontabsize) );
      makeflatcontab(data->adtree, 0, family, strides, nfamily, flatcontab, data->data, prior->arity);
      ll = ll_score_flat(flatcontab, flatcontabsize, family, nfamily, child, prior->arity);
      SCIPfreeBlockMemoryArray(extras->scip, &flatcontab, flatcontabsize);
   }
   SCIPfreeBufferArray(extras->scip, &strides);
   SCIPfreeBufferArray(extras->scip, &family);
   return ll - forchild->penalty * nparentinsts;
}

#if 0
/** The input to this function is a set of counts n1, n2, ... nr and a positive real alpha.
 * alpha defines a Dirichlet distribution with r parameters, all equal to alpha. There is a marginal probability
 * of observing the counts n1, n2, ... nr given this Dirichlet, a value that can be viewed as a function of alpha: f(alpha).
 * @return This function returns TRUE if it establishes that f'(alpha) >= 0 and FALSE otherwise.
 * A return value of FALSE does not guarantee that f'(alpha) < 0. FALSE can be returned if determining whether f'(alpha) >= 0
 * would take too long.
 */
static
SCIP_Bool positive_gradient(
   COUNT* counts,    /**< the set of counts */
   int ncounts,      /**< the number of counts */
   COUNT total,      /**< the sum of the counts */
   SCIP_Real alpha   /**< Common value for all `ncounts` Dirichlet parameters */
   )
{

   int i;
   SCIP_Real chisq;
   SCIP_Real mean;
   SCIP_Real diff;

   assert( counts != NULL );
   assert( ncounts > 1 );
   assert( alpha > 0.0 );

   /* shouldn't have all zero counts */
   assert( total > 0 );

   chisq = 0.0;
   mean = ((SCIP_Real) total) / ncounts;
   for( i = 0; i < ncounts; i++)
   {
      diff = counts[i] - mean;
      chisq += (diff*diff);
   }
   chisq /= mean;

   /* From Levin & Reeds
    * NB. This does not depend on alpha
    */
   if( chisq < ncounts - 1 )
      return TRUE;

   /* sort counts into non-increasing order */
   SCIPsortDownInt((int*)counts,ncounts);

   /* always negative gradient when there is only one positive count */
   if( counts[1] == 0 )
      return FALSE;

   /* common special case */
   if( ncounts == 2 && alpha <= 0.5 )
   {
      /* Validity of following inference established by some
         offline computation. This could be generalised with some more work!
      */
      if( counts[1] - counts[0] < 20 )
      {
         /* not too close to a 'deterministic' distribution */
         return TRUE;
      }
   }

   /* other cases where gradient can be shown to be positive TODO */

   /* failed to establish a positive gradient ... */
   return FALSE;

}
#endif

#if 0
/** Used to compute an upper bound for a particular instantiation of some parent set
 * for some child, from a contingency table whose variables are that child and all potential additional
 * parents with higher index than any of the current parents. Counts in the contingency table
 * are only for datapoints satisfying the particular instantation (which is not explicitly
 * represented). The child variable is always the last variable in the contingency table.
 * @return The upper bound
 */
static SCORE local_bdeu_ub(
   const TREECONTAB treecontab,  /**< tree-shaped contingency table */
   const VARIABLE *variables,    /**< variables for the contingency table */
   COUNT nvariables,             /**< number of variables in the contingency table */
   const ARITY* arity,           /**< arity[i] is arity of variable i */
   SCIP_Real alpha,              /**< ESS divided by number of parent set instantiations */
   SCORE* bestfirstdiff          /**< best improvement by replacing a general-purpose upper bound
                                    for child counts by that obtained by having those child counts
                                    first in the chain rule */
   )
{
   SCORE ub = 0.0;
   VALUE val;
   const int nvals = arity[variables[0]];
   assert(nvariables > 0);
   assert(variables != NULL);
   assert(arity != NULL);

   /* NULL value indicates a contingency table where all counts are zero */
   if( treecontab.children == NULL )
      return 0.0;

   /* print_treecontabptr(&treecontab, variables, nvariables, arity); */

   if( nvariables == 1 )
   {

      /* since the child variable is always last variable in the (tree-shaped)
         contingency table the only variable in this contingency table is the child variable
      */
      COUNT* counts;
      COUNT nij;
      COUNT nijk;


      nij = 0;


      for(val = 0; val < nvals; val++ )
         nij += treecontab.children[val].count;

      /* special case of no counts even though treecontab.children not NULL */
      if( nij == 0 )
         return 0.0;

      /* TODO avoid this malloc */
      counts = (COUNT*) malloc(nvals * sizeof(COUNT));
      for(val = 0; val < nvals; val++ )
      {
         nijk = treecontab.children[val].count;
         counts[val] = nijk;
         if( nijk > 0 )
         {
            /* simple MLE-based upper bound */
            ub += nijk * log(nijk / (SCIP_Real) nij);
            /* printf("nijk=%d,nij=%d,ubinc=%g\n",nijk,nij,nijk * log(nijk / (SCIP_Real) nij)); */
         }
      }

      /* if( positive_gradient(counts,nvals,nij,alpha) ) */
      /* { */
      /*    /\* TODO, compute local local BDeu score for supersets and see whether better than 'ub' */
      /*       computed above  */
      /*    *\/ */
      /*    printf("TODO\n"); */
      /* } */
      free(counts);
   }
   else
   {
      for(val = 0; val < nvals; val++ )
         ub += local_bdeu_ub(treecontab.children[val], variables+1, nvariables-1, arity, alpha, bestfirstdiff);
   }
   /* printf("some ub = %g\n",ub); */
   return ub;
}
#endif

#if 0
/** Compute an upper bound on the BDeu score of (a subset of ) the supersets of a parent set from a tree-shaped
 * contingency table. The parents come first in the tree, and the child
 * comes last. In between are all non-child variables with a higher index than any parent.
 * The upper bound is on all superset parent sets formed by adding parents with a higher
 * index than any of the current parents
 * @return The upper bound
 */
static SCORE add_bdeu_ubs(
   const TREECONTAB treecontab,  /**< tree-shaped contingency table */
   const VARIABLE *variables,    /**< variables for the contingency table */
   COUNT nvariables,             /**< number of variables in the contingency table */
   const ARITY* arity,           /**< arity[i] is arity of variable i */
   int nparents,                 /**< the first nparents variables in the tree are the current parents */
   SCIP_Real alpha               /**< ESS divided by number of parent set instantiations */
   )
{

   assert(nvariables > 0);
   assert(treecontab.children != NULL);
   assert(variables != NULL);
   assert(arity != NULL);
   assert(nparents < (int) nvariables);

   if( nparents == 0 )
   {
      /* printf("START HERE\n"); */
      SCORE bestfirstdiff = 0.0;
      return local_bdeu_ub(treecontab, variables, nvariables, arity, alpha, &bestfirstdiff);
   }
   else
   {
      SCORE ub;
      VALUE val;

      ub = 0.0;
      for(val = 0; val < arity[variables[0]]; val++ )
         if( treecontab.children[val].children != NULL )
            ub += add_bdeu_ubs(treecontab.children[val], variables+1, nvariables-1, arity, nparents-1, alpha);
      return ub;
   }
}
#endif

#if 0
/** Compute an upper bound on the BDeu scores of supersets of given parent set for given child
 * where only parent sets consisting of new parents of higher index than any in given parent set are considered
 * @return Upper bound on supserset BDeu scores
 */
static
SCORE bde_score_upperbound(
   VARIABLE child,          /**< the child of the family */
   const VARIABLE* parents, /**< the parents of the family */
   int nparents,            /**< the number of parents in the family */
   const PRIOR* prior,      /**< prior information for scoring, e.g. hyperparameters */
   const POSTERIOR* data,   /**< the data */
   SCORING_EXTRAS* extras   /**< any necessary extra arguments required for scoring */
   )
{
   int i;
   int j;
   VARIABLE last_parent;
   VARIABLE jvar;
   VARIABLE* full_family;
   int full_family_size;
   SCIP_Bool early_child_to_insert;
   int childpos;
   SCORE ub;
   int flatcontabsize;
   int* strides;
   FLATCONTAB flatcontab;

   last_parent = parents[nparents-1];

   /* compute size of 'full family' when child comes after (has higher index)
      than the last parent
    */
   full_family_size = nparents + (prior->nvars - (int) last_parent - 1);

   /* if( full_family_size > 0 ) */
   /*    return 0.0; */


   /* if child comes before the last parent then it will not be included
      in the set of later variables, so need to correct value for size of full family */
   if( child < last_parent )
   {
      early_child_to_insert = TRUE;
      full_family_size++;
   }
   else
      early_child_to_insert = FALSE;

   SCIP_CALL( SCIPallocBufferArray(extras->scip, &full_family, full_family_size) );

   /* construct full family so that variables are in order
      ensure that the child variable is inserted correctly
      and its position recorded
   */
   i = 0;
   j = 0;
   while( j < nparents )
   {
      if( early_child_to_insert && child < parents[j] )
      {
         childpos = i;
         full_family[i++] = child;
         early_child_to_insert = FALSE;
      }
      else
         full_family[i++] = parents[j++];
   }
   jvar = last_parent + 1;
   while( i < full_family_size )
   {
      if( jvar == child)
         childpos = i;
      full_family[i++] = jvar++;
   }

   /* compute size for a flat contingency table (and associated 'strides') */
   SCIP_CALL( SCIPallocBufferArray(extras->scip, &strides, full_family_size) );
   flatcontabsize = 1;
   for( i = full_family_size - 1; flatcontabsize <= extras->maxflatcontabsize && i >= 0; i--)
   {
      strides[i] = flatcontabsize;
      flatcontabsize *= (prior->arity)[full_family[i]];
   }

   ub = 0.0;

   /* only compute a useful upper bound if needed contingency table
    * would not be too big
    */
   if( flatcontabsize <= extras->maxflatcontabsize )
   {

      COUNT nij;
      COUNT nijk;
      int k;

      int abovechildstride;
      int childstride;

      SCORE thisub;

      childstride = strides[childpos];
      abovechildstride = childstride * (prior->arity)[child];
      /* entries in flatcontab childstride apart only differ in the value of the child variable */

      SCIP_CALL( SCIPallocClearBufferArray(extras->scip, &flatcontab, flatcontabsize) );
      makeflatcontab(data->adtree, 0, full_family, strides, full_family_size, flatcontab, data->data, prior->arity);

      /* iterate through insts which differ in values of variables 'above' the child */
      for( i = 0; i < flatcontabsize; i += abovechildstride )
      {
         /* each value of j should correspond to the first value of the child for a particular inst of
            all the other variables */
         for( j = i; j < i + childstride; j++ )
         {
            /* find total of counts */
            nij = 0.0;
            for( k = j; k < j + abovechildstride; k += childstride )
            {
               nij += flatcontab[k];
            }
            if( nij > 0 )
            {
               thisub = 0.0;
               for( k = j; k < j + abovechildstride; k += childstride )
               {
                  nijk = flatcontab[k];
                  if( nijk > 0 )
                     /* simple MLE-based upper bound */
                     thisub += nijk * log(nijk / (SCIP_Real) nij);
               }
               ub += thisub;
            }
         }
      }
      SCIPfreeBufferArray(extras->scip, &flatcontab);
   }
   SCIPfreeBufferArray(extras->scip, &strides);
   SCIPfreeBufferArray(extras->scip, &full_family);

   /* printf("ub is %g\n", ub); */

   return ub;
}
#endif


/**< Compute the BDeu score for given child and parent set
 * @return BDeu score
 */
static
SCORE bde_score(
   VARIABLE child,          /**< the child of the family */
   const VARIABLE* parents, /**< the parents of the family */
   int nparents,            /**< the number of parents in the family */
   const PRIOR* prior,      /**< prior information for scoring, e.g. hyperparameters */
   const POSTERIOR* data,   /**< the data */
   SCORING_EXTRAS* extras   /**< any necessary extra arguments required for scoring */
   )
{

   SCIP* scip = extras->scip;
   VARIABLE* family;
   int nfamily = nparents + 1;
   SCORE skore;

   assert(nparents > -1);
   assert(parents != NULL || nparents == 0);
   assert(prior != NULL);
   assert(data != NULL);
   assert(extras != NULL);

   /* make family by inserting child into parents */
   SCIP_CALL( SCIPallocBufferArray(scip, &family, nfamily) );
   makefamily(family,parents,nparents,child);

   /* because log_likelihood for family is computed AFTER that for parents,
      when this function returns extras->npos_cells_ptr will contain
      a pointer to the number of positive counts in the contingency table for
      the family, which is what we want
   */
   /* improvement thanks to Cassio de Campos */

   skore =  (
      -log_likelihood_cache_compact(parents,nparents,prior,data,extras)
      +log_likelihood_cache_compact(family,nfamily,prior,data,extras)) ;

   SCIPfreeBufferArray(scip, &family);

   return skore;
}

/** Compute a local score for given child and parents
 * and an upper bound on local scores of all supersets of this parent set
 * @return local score
*/
static
SCORE local_score(
   SCIP* scip,               /**< SCIP instance */
   SCORE_CODE score_type,    /**< which local score? (either SCORE_BDEU, SCORE_BGE or SCORE_BIC )*/
   VARIABLE child,           /**< the child of the family */
   const VARIABLE* parents,  /**< the parents of the family */
   int nparents,             /**< the number of parents in the family */
   const PRIOR* prior,       /**< prior information for scoring, e.g. hyperparameters */
   const POSTERIOR* data,    /**< the data */
   SCORING_EXTRAS* extras,   /**< any necessary extra arguments required for scoring */
   const FORCHILD* forchild, /**< useful quantities for this child */
   SCORE* ub                 /**< *ub is an upper bound on local scores of supersets of this parent set */
   )
{

   double skore;
   int i;

#ifdef BLAS
   VARIABLE* family;
   VARIABLE* ordered_family;
#endif

   assert(scip != NULL);
   assert(nparents == 0 || parents != NULL);
   assert(prior != NULL);
   assert(data != NULL);
   assert(extras != NULL);

   if( score_type == SCORE_BDEU )
   {
      skore = bde_score(child, parents, nparents, prior, data, extras);
      *ub = forchild->neglogaritychild * extras->npos_cells;
      /* printf("naive=%g ",*ub); */
      /* if( nparents > 0 )  */
      /*    *ub = MIN(*ub,bde_score_upperbound(child, parents, nparents, prior, data, extras)); */
      /* printf("after=%g\n",*ub); */
      return skore;
   }
   else if( score_type == SCORE_BIC )
   {
      ARITY nparentinsts = 2; /* (*) */
      for(i = 0; i < nparents; i++)
         nparentinsts *= prior->arity[i];
      /* any superset parent set will have at least nparentinsts parent instantiations
         since 2 (*) is smallest arity possible, ll_score is upperbounded by 0,
         so forchild->penalty * nparentinsts upper bounds best possible score for superset
         parent sets */
      *ub = -forchild->penalty * nparentinsts;
      return bic_score(child, parents, nparents, prior, data, extras, forchild);
   }
   else if(score_type == SCORE_BGE )
   {
#ifdef BLAS

      assert(data->log_gamma_ratio_table != NULL);
      assert(prior->prior_matrix != NULL);
      assert(data->posterior_matrix != NULL);
      assert(data->continuous_data != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &family,  nparents+1) );
      family[0] = child;
      for(i = 0; i < nparents; i++)
         family[i+1] = parents[i];

      SCIP_CALL( SCIPallocBufferArray(scip, &ordered_family,  nparents+1) );
      makefamily(ordered_family,parents,nparents,child);

      skore = LogBgeScore(child, family, ordered_family, nparents, prior->alpha_mu, prior->alpha_omega, data->log_prefactor, data->log_gamma_ratio_table,
         prior->prior_matrix, data->posterior_matrix, data->continuous_data, extras->scorecache, scip);
      SCIPfreeBufferArray(scip, &family);
      /* UPPER BOUND NOT CURRENTLY COMPUTED FOR BGE!! */
      *ub = SCIPinfinity(scip);
      return (SCORE) skore;
#else
      SCIPerrorMessage("You can't use BGe scoring without BLAS support\n");
      return 0;
#endif
   }
   else
   {
      SCIPerrorMessage("Unrecognised score function\n");
      return 0;
   }
}

/** Free the data associated with posterior quantities, e.g. ADtrees, the data itself, etc
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE freedata(
   SCIP* scip,                    /**< SCIP instance */
   SCORINGPARAMS* scoringparams,  /**< scoring parameters */
   PRIOR* prior,                  /**< prior quantities */
   POSTERIOR* data                /**< posterior quantities, e.g. the data */
   )
{
   assert(scoringparams->score_type != SCORE_BGE || scoringparams->continuous);
   assert(scoringparams->score_type != SCORE_BDEU || !scoringparams->continuous);
   assert(scoringparams->score_type != SCORE_BIC || !scoringparams->continuous);

   if( scoringparams->continuous )
   {
      if( scoringparams->score_type == SCORE_BGE )
      {
#ifdef BLAS
         BgeMatrixDelete(&(data->continuous_data));
         BgeMatrixDelete(&(data->posterior_matrix));
         free(data->log_gamma_ratio_table);

#else
         SCIPerrorMessage("You can't use BGe scoring without BLAS support\n");
         return SCIP_ERROR;
#endif
      }
   }
   else
   {
      int i;

      assert(prior->arity != NULL);

      delete_adtree(data->adtree, 0, prior->nvars, prior->arity);
      for(i = 0; i < prior->nvars; ++i)
         SCIPfreeMemoryArray(scip, &((data->data)[i]));
      SCIPfreeMemoryArray(scip, &(data->data));
      SCIPfreeMemoryArray(scip, &(data->num_observed));
   }
   return SCIP_OKAY;
}

/** Set the data associated with posterior quantities, e.g. ADtrees, the data itself, etc
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE setdata(
   SCIP* scip,                   /**< SCIP instance */
   SCORINGPARAMS* scoringparams, /**< scoring parameters */
   PRIOR* prior,                 /**< prior quantities */
   POSTERIOR* data               /**< posterior quantities, e.g. the data */
   )
{

   assert(scoringparams->score_type != SCORE_BGE || scoringparams->continuous);
   assert(scoringparams->score_type != SCORE_BDEU || !(scoringparams->continuous));
   assert(scoringparams->score_type != SCORE_BIC || !(scoringparams->continuous));

   data->logndiv2 = log(data->nrows) / 2;

   if( scoringparams->continuous )
   {
      if( scoringparams->score_type == SCORE_BGE )
      {
#ifdef BLAS
         data->posterior_matrix = BgeMatrixCreate(prior->nvars, prior->nvars);
         data->log_prefactor = 0.5 * (log(prior->alpha_mu) - log(prior->alpha_mu + data->nrows));
         data->log_gamma_ratio_table = create_log_gamma_ratio_table(data->nrows, prior->alpha_omega, prior->nvars);
         SetPosteriorParametricMatrix(data->continuous_data, prior->prior_matrix, data->posterior_matrix,
            prior->alpha_mu, prior->alpha_omega, prior->nu);
#else
         SCIPerrorMessage("You can't use BGe scoring without BLAS support\n");
         return SCIP_ERROR;
#endif
      }
      else
      {
         SCIPerrorMessage("\nFor continuous data only BGe scoring is currently available. You need to explicitly specify BGe scoring.\nAdd the following line:\ngobnilp/scoring/score_type = \"BGe\"\nto your settings file\n");
         return SCIP_ERROR;
      }
   }
   else
   {
      /* Build the AD tree, not specific to a particular (discrete) score type */

      int i;
      ADTREE* adtree = malloc(sizeof(ADTREE));
      ROW* allrows;
      int rmin;
      int adtreedepthlim;
      int adtreenodeslim;
      int n_nodes;

      allrows = (ROW *) malloc(data->nrows * sizeof(ROW));
      for( i = 0; i < data->nrows; ++i )
         allrows[i] = i;

      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/rmin", &rmin) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/adtreedepthlim", &adtreedepthlim) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/adtreenodeslim", &adtreenodeslim) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Building the AD tree ...");
      /* build_adtree will free allrows */
      build_adtree(adtree, 0, allrows, data->nrows, rmin, 0, &n_nodes,
         adtreedepthlim, adtreenodeslim, prior->nvars, data->data, prior->arity);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, " Done\n");
      data->adtree = adtree;
   }
   return SCIP_OKAY;
}

/**< free memory for extra scoring information */
static
void freeextras(
   SCIP* scip,                    /**< SCIP instance */
   SCORINGPARAMS* scoringparams,  /**< scoring parameters */
   SCORING_EXTRAS* extras         /**< extra information for scoring */
   )
{


   assert(scip != NULL );
   assert(scoringparams != NULL );
   assert(extras != NULL );

   /* no extras to free if using BIC */
   if( scoringparams->score_type == SCORE_BIC )
      return;

   delete_scorecache(scip, extras->scorecache);
   SCIPfreeBlockMemory(scip, &(extras->scorecache));

}

/** set extra information for scoring
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE setextras(
   SCIP* scip,                    /**< SCIP instance */
   SCORINGPARAMS* scoringparams,  /**< scoring parameters */
   PRIOR* prior,                  /**< prior quantities */
   POSTERIOR* data,               /**< posterior quantities, e.g. the data */
   SCORING_EXTRAS* extras         /**< extra scoring information */
   )
{

   int cachesizelimit;
   int cacheblocksize;
   int nvarscachelimit;

   assert(scip != NULL);
   assert(scoringparams != NULL);
   assert(prior != NULL);
   assert(data != NULL);
   assert(extras != NULL);

   assert(scoringparams->score_type != SCORE_BGE || scoringparams->continuous);
   assert(scoringparams->score_type != SCORE_BDEU || !(scoringparams->continuous));
   assert(scoringparams->score_type != SCORE_BIC || !(scoringparams->continuous));

   extras->scip = scip;
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/maxflatcontabsize", &(extras->maxflatcontabsize)) );

   if( scoringparams->score_type == SCORE_BDEU )
      SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/fast", &(extras->fast)) );

   if( scoringparams->score_type == SCORE_BDEU || scoringparams->score_type == SCORE_BGE )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &(extras->scorecache)) );

      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/cachesizelimit", &cachesizelimit) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/cacheblocksize", &cacheblocksize) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/nvarscachelimit", &nvarscachelimit) );

      SCIP_CALL( make_scorecache(scip, extras->scorecache, prior->nvars,
            nvarscachelimit, cachesizelimit, cacheblocksize) );
   }
   return SCIP_OKAY;
}

/** free memory for prior information
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE freeprior(
   SCIP* scip,                   /**< SCIP instance */
   SCORINGPARAMS* scoringparams, /**< scoring parameters */
   PRIOR* prior                  /**< prior quantities */
   )
{
   int i;
   int j;

   for( i = 0; i < prior->nvars; ++i )
   {
      SCIPfreeMemoryArray(scip, &(prior->is_parent[i]) );
      SCIPfreeMemoryArray(scip, &(prior->is_parents[i]) );
   }
   SCIPfreeMemoryArray(scip, &(prior->n_is_parents) );
   SCIPfreeMemoryArray(scip, &(prior->is_parents) );
   SCIPfreeMemoryArray(scip, &(prior->is_parent) );

   if( !scoringparams->continuous )
   {
      for(i = 0; i < prior->nvars; ++i)
      {
         for(j = 0; j < (prior->arity)[i]; ++j)
            SCIPfreeMemoryArray(scip, &((prior->labels)[i][j]));
         SCIPfreeMemoryArray(scip, &((prior->labels)[i]));
      }
      SCIPfreeMemoryArray(scip, &(prior->labels));
      SCIPfreeMemoryArray(scip, &(prior->arity));
   }
   else if (scoringparams->score_type == SCORE_BGE)
#ifdef BLAS
   {
      BgeMatrixDelete(&(prior->prior_matrix));
      if( prior->nu != NULL )
         BgeVectorDelete(&(prior->nu));
   }
#else
   {
      SCIPerrorMessage("You can't use BGe scoring without BLAS support\n");
      return SCIP_ERROR;
   }
#endif
   return SCIP_OKAY;
}

/** set prior quantities
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE setprior(
   SCIP* scip,                   /**< SCIP instance */
   SCORINGPARAMS* scoringparams, /**< scoring parameters */
   PRIOR* prior                  /**< prior quantities */
   )
{
   int i;
   int j;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(prior->is_parent), prior->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prior->is_parents), prior->nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(prior->n_is_parents), prior->nvars) );
   for( i = 0; i < prior->nvars; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(prior->is_parent[i]), prior->nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(prior->is_parents[i]), prior->nvars) );
      (prior->n_is_parents)[i] = 0;
      for( j = 0; j < prior->nvars; ++j )
         (prior->is_parent)[i][j] = 0;
   }

   /* find any parents/nonparents */
   SCIP_CALL( addGeneralDAGConstraints(scip, prior->nvars, prior->nodeNames, prior->is_parent) );

   /* set is_parents */
   for( i = 0; i < prior->nvars; ++i )
   {
      for( j = 0; j < prior->nvars; ++j )
         if( (prior->is_parent)[i][j] == 1)
         {
            (prior->is_parents)[i][((prior->n_is_parents)[i])++] = j;
         }
      SCIP_CALL( SCIPreallocMemoryArray(scip, &((prior->is_parents)[i]), (prior->n_is_parents)[i]) );
   }

#ifdef SCIP_DEBUG
   for( i = 0; i < prior->nvars; ++i )
   {
      for( j = 0; j < prior->nvars; ++j )
      {
         if( (prior->is_parent)[i][j] == -1)
            SCIPdebugMessage("%d is a banned parent for %d\n", j, i);
         if( (prior->is_parent)[i][j] == 1)
            SCIPdebugMessage("%d is a required parent for %d\n", j, i);
      }
      for( j = 0; j < (prior->n_is_parents)[i]; ++j )
         SCIPdebugMessage("%d is a required parent for %d\n", (prior->is_parents)[i][j], i);
   }
#endif


   if( scoringparams->score_type == SCORE_BDEU)
   {
      SCIP_CALL(SCIPgetRealParam(scip, "gobnilp/scoring/alpha", &(prior->alpha)) );
   }
   else if (scoringparams->score_type == SCORE_BGE)
   {
#ifdef BLAS
      int alpha_omega_minus_nvars;
      char* nu_str;
      char* t_str;

      SCIP_CALL( SCIPgetRealParam(scip, "gobnilp/scoring/alpha_mu", &(prior->alpha_mu)) );
      SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/alpha_omega_minus_nvars", &alpha_omega_minus_nvars) );
      prior->alpha_omega = prior->nvars + alpha_omega_minus_nvars;

      SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/scoring/nu", &nu_str) );
      if( strcmp(nu_str,"") == 0 )
      {
         /* sample mean will be used */
         prior->nu = NULL;
      }
      else
      {
         const char delim[1] = ",";
         char* token;

         prior->nu = BgeVectorCreate(prior->nvars);

         /* get first value */
         token = strtok(nu_str,delim);
         if( token == NULL )
         {
            SCIPerrorMessage("Could not find component 0 of nu prior vector\n");
            return SCIP_ERROR;
         }
         prior->nu->items[0] = atof(token);

         /* get other values */
         for(i = 1; i < prior->nvars; i++)
         {
            token = strtok(NULL,delim);
            if( token == NULL )
            {
               SCIPerrorMessage("Could not find component %d of nu prior vector\n",i);
               return SCIP_ERROR;
            }
            prior->nu->items[i] = atof(token);
         }
      }

      prior->prior_matrix = BgeMatrixCreate(prior->nvars,prior->nvars);
      SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/scoring/T", &t_str) );
      /* currently ignore value of T and always use default */
      if( TRUE || strcmp(t_str,"") == 0 )
      {
         /* Use 'default' prior matrix */
         SetPriorParametricMatrix(prior->nvars, prior->alpha_mu, prior->alpha_omega, prior->prior_matrix);
      }


#else
      SCIPerrorMessage("You can't use BGe scoring without BLAS support\n");
      return SCIP_ERROR;
#endif
   }
   return SCIP_OKAY;
}

/** checks whether all required parents are in a parent set
 *
 * @return TRUE if all required parents present in parent set, else FALSE
 */
static
SCIP_Bool requiredparentspresent(
   VARIABLE child,                 /**< child */
   VARIABLE* ps,                   /**< parent set */
   int npa,                        /**< size of parent set */
   PRIOR* prior                    /**< prior information e.g. required parent sets */
   )
{
   int i;
   int j;
   SCIP_Bool ok;

   /* check that all required parents are there */
   for( i = 0; i < prior->n_is_parents[child]; ++i )
   {
      ok = FALSE;
      for( j = 0; j < npa; ++j )
      {
         if( ps[j] == prior->is_parents[child][i] )
         {
            ok = TRUE;
            break;
         }
      }
      if( !ok )
      {
         /* ith required parent is missing */
         return FALSE;
      }
   }
   /* all there */
   return TRUE;
}

/** compute local scores for a particular child by extending an existing parent set
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE makelocalscores_child_leaf(
   SCIP* scip,                     /**< SCIP instance */
   PRIOR* prior,                   /**< prior quantities */
   POSTERIOR* data,                /**< posterior quantities, e.g. the data */
   SCORING_EXTRAS* extras,         /**< extra information for scoring */
   SCORINGPARAMS* scoringparams,   /**< scoring parameters */
   VARIABLE child,                 /**< child for local scores */
   FORCHILD* forchild,             /**< stores prior information specific to a child */
   int npa,                        /**< size of new parent set */
   int* nnodes,                    /**< (*nnodess) is the number of trie nodes for parent sets of size npa for this child */
   int* nodessize,                 /**< (*nodessize) is the length of the array layer */
   TRIENODE*** layer,              /**< (*layer)[i] is the ith trie node of size npa for this child */
   int* nkeptall,                  /**< total number of kept nodes */
   TRIENODE* leaf                  /**< leaf to extend ( corresponds to a parent set of size npa-1 ) */
   )
{

   VARIABLE lb;                 /* lower bound (=<) for creating new parent sets */
   const int noldpa = npa - 1;
   VARIABLE newpa;              /* a new parent added to an old parent set to make a new one */
   SCIP_Bool last_layer;        /* indicates whether current layer is the final one */
   SCORE bestsubsetscore;       /* score of best scoring subset of a parent set */
   SCORE lowest_ub;             /* score of lowest upper bound of supersets of a parent set */
   TRIENODE* nodebefore;
   SCIP_Bool first;

   VARIABLE* ps;                /* parent sets of size npa */
   SCORE skore;                 /* local score */
   SCORE ub;                    /* upper bound on superset local scores */
   SCIP_Bool keep;              /* whether to keep (i.e. not prune) a scored parent set */
   TRIENODE* node;              /* a new node */

   assert(leaf->child == NULL);
   assert(npa > 0);

   /* used both during initial scoring and also pricing */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM || SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   SCIP_CALL( SCIPallocBufferArray(scip, &ps,  npa) );
   getset(leaf,ps,noldpa);

   last_layer = ( npa == scoringparams->palim );

   /* extend old parent set to create new parent sets
      by adding new parents strictly greater
      than any in the old parent set
   */

   first = TRUE;
   lb = (leaf->mother == NULL) ? 0 : leaf->elt + 1;
   for( newpa = lb; (int) newpa < prior->nvars; ++newpa )
   {
      /* just skip over parents forbidden by, say,
         prior constraints */
      if( newpa == child || prior->is_parent[child][newpa] == -1)
         continue;

      /* if pruning then look for all subsets with cardinality
         just one less than current parent set
         if any missing this must be due to exponential
         pruning so don't score, otherwise return the best score of
         all subsets, and also lowest upper bound (currently available) on this parent set's supersets
      */
      if( scoringparams->pruning && exppruned(scip,leaf,noldpa,newpa,scoringparams->prunegap,&bestsubsetscore,&lowest_ub) )
         continue;

      /* create new parent set by adding new parent to old parent set */
      ps[noldpa] = newpa;

      /* score the parent set
         upon return 'extras' may contain additional useful information
         for example, when BDeu scoring, extras->npos_cells is the number of positive cells in
         the contingency table for the family */
      skore = local_score(scip,scoringparams->score_type,child,ps,npa,prior,data,extras,forchild,&ub);

      /* printf("child=%d,parents=", child); */
      /* for( i = 0; i < npa; i++ ) */
      /*    printf("%d,", ps[i]); */
      /* printf(":skore=%g,ub=%g\n", skore, ub); */


      /* just in case we have a better (i.e. lower upper bound) from subsets
         (NB. a good upper bound finding routine should ensure that ub <= lowest_ub, anyway ) */
      /* Only if we are pruning do we even have a proper value for lowest_ub */
      if( scoringparams->pruning )
         ub = MIN(ub,lowest_ub);


      /* check that upper bound used for 'exponential' pruning is correct */
      assert(skore <= ub - scoringparams->prunegap );

      if( scoringparams->pruning && skore <= bestsubsetscore + scoringparams->prunegap )
      {
         /* score not good enough to keep */

         if( last_layer || bestsubsetscore + scoringparams->prunegap > ub  )
         {
            /* this is the last layer or parent set was
               exponentially pruned, don't create node, continue to next new parent */
            continue;
         }
         else
         {
            /* create node since supersets might be worth keeping
               but won't keep this parent set to make a family variable
               overwrite score with score of best scoring subset parent set */
            skore = bestsubsetscore;
            keep = FALSE;
         }
      }
      else
      {
         /* only keep parent set if all required parents are present in this
            parent set */
         keep = requiredparentspresent(child, ps, npa, prior);
         if( keep )
            (*nkeptall)++;
      }

      /* create node and store it in current layer
         ( allocating more space if required ) */
      SCIP_CALL( makenode(scip,&node,newpa,keep,skore,ub,leaf) );
      if( *nodessize <= *nnodes )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, layer, *nodessize, *nodessize + BLOCKSIZE) );
         *nodessize += BLOCKSIZE;
      }
      (*layer)[(*nnodes)++] = node;


      /* if this node is the first extension of the leaf
         make it the child of the leaf, else ensure the previous
         extension of leaf points to this node as 'next' (sibling)
      */
      if( first )
      {
         leaf->child = node;
         first = FALSE;
      }
      else
      {
         nodebefore->next = node;
      }
      nodebefore = node;

   }
   SCIPfreeBufferArray(scip,&ps);
   return SCIP_OKAY;
}

/** compute local scores for a particular child and particular parent set size
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE makelocalscores_child_npa(
   SCIP* scip,                     /**< SCIP instance */
   PRIOR* prior,                   /**< prior quantities */
   POSTERIOR* data,                /**< posterior quantities, e.g. the data */
   SCORING_EXTRAS* extras,         /**< extra information for scoring */
   SCORINGPARAMS* scoringparams,   /**< scoring parameters */
   VARIABLE child,                 /**< child for local scores */
   FORCHILD* forchild,             /**< stores prior information specific to a child */
   int npa,                        /**< size of parent set */
   int* nnodess,                   /**< nnodess[noldpa] is the number of trie nodes for parent sets of size noldpa for this child */
   TRIENODE*** layers,             /**< layers[noldpa][i] is the ith trie node of size noldpa for this child */
   int* nkeptall                   /**< total number of kept nodes */
   )
{

   const int noldpa = npa - 1;      /* size of an old parent set */
   int nodessize;
   int nnodes;                      /* number of nodes created in this layer */
   int i;
   SCORE ub;

   /* only to be used during initial scoring */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM);

   /* deal with empty parent set as a special case */
   if( npa == 0 )
   {

      SCORE skore = local_score(scip,scoringparams->score_type,child,NULL,0,prior,data,extras,forchild,&ub);
      /* empty parent set always kept unless prior disallows it */
      SCIP_Bool keep = requiredparentspresent(child, NULL, 0, prior );
      TRIENODE* node;

      if( keep )
         (*nkeptall)++;

      /* root node corresponds to empty parent set
       * it has no last element, so put in prior->nvars as a dummy value */
      SCIP_CALL( makenode(scip,&node,prior->nvars,keep,skore,ub,NULL) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(layers[0]), 1) );
      layers[0][0] = node;
      nnodess[0] = 1;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(layers[npa]), BLOCKSIZE) );
   nnodes = 0;
   nodessize = BLOCKSIZE;

   /* extend each parent set ('leaf') from the previous layer */
   for(i = 0; i < nnodess[noldpa]; i++)
   {
      SCIP_CALL( makelocalscores_child_leaf(scip, prior, data, extras, scoringparams,
            child, forchild, npa, &nnodes, &nodessize, &(layers[npa]), nkeptall, layers[noldpa][i]) );
   }
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(layers[npa]), nodessize, nnodes) );
   nnodess[npa] = nnodes;


   return SCIP_OKAY;
}

#if 0
/** create new nodes by adding parents to nodes on an existing open list */
static
SCIP_RETCODE makelocalscores_child_openlist(
   SCIP* scip,                     /**< SCIP instance */
   PRIOR* prior,                   /**< prior quantities */
   POSTERIOR* data,                /**< posterior quantities, e.g. the data */
   SCORING_EXTRAS* extras,         /**< extra information for scoring */
   SCORINGPARAMS* scoringparams,   /**< scoring parameters */
   VARIABLE child,                 /**< child for local scores */
   FORCHILD* forchild,             /**< stores prior information specific to a child */
   TRIENODE** openlist,            /**< the open list */
   int nopenlist,                  /**< number of nodes in the open list */
   TRIENODE*** newnodes,           /**< (*newnodes) is the list of new nodes */
   int* nnewnodes                  /**< number of new nodes */
   )
{

   int newnodessize = BLOCKSIZE;
   int i;
   /* not using nkeptall at present */
   int nkeptall = 0;

   /* only for use during price-and-cut loop */
   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, newnodes, BLOCKSIZE) );
   *nnewnodes = 0;

   /* extend each parent set ('leaf') from the open list */
   for(i = 0; i < nopenlist; i++)
   {
      TRIENODE* node = openlist[i];

      SCIP_CALL( makelocalscores_child_leaf(scip, prior, data, extras, scoringparams,
            child, forchild, getnpa(node), nnewnodes, &newnodessize, newnodes, &nkeptall, node) );
   }
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, newnodes, newnodessize, *nnewnodes) );

   return SCIP_OKAY;
}
#endif

/** compute local scores for a particular child
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE makelocalscores_child(
   SCIP* scip,                     /**< SCIP instance */
   PRIOR* prior,                   /**< prior quantities */
   POSTERIOR* data,                /**< posterior quantities, e.g. the data */
   SCORING_EXTRAS* extras,         /**< extra information for scoring */
   SCORINGPARAMS* scoringparams,   /**< scoring parameters */
   ParentSetData* psd,             /**< to be populated, does not include scores */
   SCIP_Real*** scores,            /**< (*scores)[child][j] will be the score for the jth parent set of child */
   int* nscores,                   /**< *nscores will contain the total number of local scores created */
   VARIABLE child                  /**< child for local scores */
   )
{

   TRIENODE*** layers;
   int* nnodess;
   int nnodes;

   FORCHILD forchild;         /* stores prior information specific to a child */
   TRIENODE** allkeptnodes;
   int* allkeptnodessizes;
   int npa;
   int i;
   int nkeptall;              /* number of kept nodes in all layers */

   /* only to be used during initial scoring */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM);

   SCIPdebugMessage("Computing local scores for child %d (%s)\n",child,prior->nodeNames[child]);

   /* set values for any useful quantities specific to this child */
   init_forchild(&forchild,scoringparams->score_type,prior,data,child);

   /* arrays long enough for any possible parent set size */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nnodess, (scoringparams->initpalim)+1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &layers, (scoringparams->initpalim)+1) );

   /* now compute local scores in layers
      parent sets in each layer are one bigger than those
      in the preceding layer */

   nkeptall = 0;
   for( npa = 0; npa <= scoringparams->initpalim; ++npa )
   {
      SCIP_CALL( makelocalscores_child_npa(scip, prior, data, extras, scoringparams,
            child, &forchild, npa, nnodess, layers, &nkeptall) );

   }
   if( nkeptall == 0 )
   {
      /* just crash out, don't bother tidying up allocated memory */
      SCIPerrorMessage("Infeasibility detected during scoring\n");
      SCIPerrorMessage("No possible parent set for %s \n", prior->nodeNames[child]);
      return SCIP_ERROR;
   }


   /* collect together all *KEPT* (to make family variable) nodes from all layers */
   SCIP_CALL( SCIPallocBufferArray(scip, &allkeptnodes, nkeptall) );
   SCIP_CALL( SCIPallocBufferArray(scip, &allkeptnodessizes, nkeptall) );
   nnodes = 0;
   /* why is this loop backwards ? */
   for( npa = scoringparams->initpalim; npa >= 0; --npa )
   {
      for( i = 0; i < nnodess[npa]; i++ )
         if( (layers[npa][i])->keep )
         {
            allkeptnodessizes[nnodes] = npa;
            allkeptnodes[nnodes++] = layers[npa][i];
         }
   }
   assert(nkeptall == nnodes);

   /* sort nodes in non-increasing order by score
      (also sort corresponding array of parent set sizes) */
   SCIPsortDownPtrInt((void**)allkeptnodes, allkeptnodessizes, nodeScoreComp, nnodes);

   /* store parent sets for this child */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d candidate parent sets for %s.\n", nkeptall, prior->nodeNames[child]);
   SCIPdebugMessage("%d candidate parent sets for child %d.\n", nkeptall, child);
   SCIP_CALL( addParentSets(scip, child, allkeptnodes, allkeptnodessizes, nkeptall, psd, scores) );
   *nscores += nkeptall;

   SCIPfreeBufferArray(scip, &allkeptnodessizes);
   SCIPfreeBufferArray(scip, &allkeptnodes);

   if( scoringparams->pricing )
   {
      SCIP_PROBDATA* probdata;

      /* get problem data */
      probdata = SCIPgetProbData(scip);
      assert( probdata != NULL );

      /* put nodes in last layer in an open list for later expansion */
      npa = scoringparams->initpalim;
      probdata->data_etc->nopenlist[child] = nnodess[npa];
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((probdata->data_etc->openlist)[child]), nnodess[npa]) );
      for( i = 0; i < nnodess[npa]; i++ )
         (probdata->data_etc->openlist)[child][i] = layers[npa][i];
   }
   else
   {
      /* clear memory in reverse order to allocation  */
      for(npa = scoringparams->initpalim; npa > -1; npa--)
      {
         for(i = nnodess[npa]-1; i > -1; i--)
            SCIPfreeBlockMemory(scip,&(layers[npa][i]));
      }
      for(npa = scoringparams->initpalim; npa > -1; npa--)
         SCIPfreeBlockMemoryArray(scip,&(layers[npa]),nnodess[npa]);
      SCIPfreeBlockMemoryArray(scip,&layers,(scoringparams->initpalim)+1 );
      SCIPfreeBlockMemoryArray(scip,&nnodess,(scoringparams->initpalim)+1 );
   }


   return SCIP_OKAY;
}

/** main function for computing local scores for each child and setting up associated parent set data structure
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE makelocalscores(
   SCIP* scip,                     /**< SCIP instance */
   PRIOR* prior,                   /**< prior quantities */
   POSTERIOR* data,                /**< posterior quantities, e.g. the data */
   SCORING_EXTRAS* extras,         /**< extra information for scoring */
   SCORINGPARAMS* scoringparams,   /**< scoring parameters */
   ParentSetData* psd,             /**< to be populated, does not include scores */
   SCIP_Real*** scores,            /**< (*scores)[child][j] will be the score for the jth parent set of child */
   int* nscores                    /**< *nscores will contain the total number of local scores created */
   )
{

   VARIABLE child;

   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM);

   /* set parent set data as much as is possible before scoring */
   psd->n = prior->nvars;
   psd->nodeNames = prior->nodeNames;
   SCIP_CALL( SCIPallocMemoryArray(scip, scores, psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->nParentSets), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->nParents), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(psd->ParentSets), psd->n) );

   *nscores = 0;

   /* compute local scores, one child at a time */
   for( child = 0; (int) child < prior->nvars; ++child )
   {
      SCIP_CALL( makelocalscores_child(scip, prior, data, extras, scoringparams, psd, scores, nscores, child) );
   }

   return SCIP_OKAY;

}

/** set global properties which record information about current problem
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
static
SCIP_RETCODE setglobalproperties(
   SCIP* scip,                   /**< SCIP instance */
   PropertyData* prop,           /**< Property data structure  */
   PRIOR* prior,                 /**< prior quantities */
   POSTERIOR* data,              /**< posterior quantities, e.g. the data */
   SCORINGPARAMS* scoringparams, /**< scoring parameters */
   const char* filename          /**< File containing the data */
   )
{
   /* Set global properties */
   SCIP_CALL( PR_setGlobalProperty(scip, prop,           "scorer_name",    "GOBNILP") );
   SCIP_CALL( PR_setGlobalProperty(scip, prop,           "scorer_version", GOBNILP_VERSION) );
   SCIP_CALL( PR_setGlobalProperty(scip, prop,           "scorer_url",     "www.cs.york.ac.uk/aig/sw/gobnilp/") );

   if( scoringparams->score_type == SCORE_BDEU )
   {
      SCIP_CALL( PR_setGlobalProperty(scip, prop,           "score_type",     "BDeu") );
      SCIP_CALL( PR_setGlobalPropertyFromReal(scip, prop,   "ess",            prior->alpha) );
   }
   else if ( scoringparams->score_type == SCORE_BIC )
   {
      SCIP_CALL( PR_setGlobalProperty(scip, prop,           "score_type",     "BIC") );
   }
   else if ( scoringparams->score_type == SCORE_BGE )
   {
#ifdef BLAS
      SCIP_CALL( PR_setGlobalProperty(scip, prop,           "score_type",     "BGe") );
      SCIP_CALL( PR_setGlobalPropertyFromReal(scip, prop,   "alpha_mu",            prior->alpha_mu) );
      SCIP_CALL( PR_setGlobalPropertyFromReal(scip, prop,   "alpha_omega",            prior->alpha_omega) );
#else
      SCIPerrorMessage("You can't use BGe scoring without BLAS support\n");
      return SCIP_ERROR;
#endif
   }
   else
      SCIP_CALL( PR_setGlobalProperty(scip, prop,           "score_type",     "**ERROR**") );

   SCIP_CALL( PR_setGlobalPropertyFromInt(scip, prop,    "parent_limit",   scoringparams->palim) );
   SCIP_CALL( PR_setGlobalPropertyFromBool(scip, prop,   "pruning",        scoringparams->pruning) );
   SCIP_CALL( PR_setGlobalPropertyFromReal(scip, prop,   "prune_gap",      scoringparams->prunegap) );
   SCIP_CALL( PR_setGlobalPropertyFromInt(scip, prop,    "num_records",    data->nrows) );
   SCIP_CALL( PR_setGlobalProperty(scip, prop,           "input_file",     filename) );

   if( !(scoringparams->continuous) )
   {
      int i;
      /* Set the individual properties */
      prop->n = prior->nvars;
      for( i = 0; i < prop->n; i++ )
      {
         SCIP_CALL( PR_setPropertyFromInt(scip, prop, i, "arity", (prior->arity)[i]) );
         SCIP_CALL( PR_setPropertyFromArray(scip, prop, i, "values", (const char**)(prior->labels)[i], (data->num_observed)[i]) );
         if( (data->num_observed)[i] < (int) (prior->arity)[i] )
            SCIP_CALL( PR_setPropertyFromInt(scip, prop, i, "num_unobserved_values", (prior->arity)[i] - (data->num_observed)[i]) );
      }
   }
   return SCIP_OKAY;
}

/** Used to ensure that the limit on parent set size is strictly less than the number of variables in the data
 * @return Lowest limit on parent set cardinality consistent with palim value supplied by user
 */
static
int fixpalim(
   int palim,  /**< parent set size limit set by user (can be -1 to indicate no limit) */
   int nvars   /**< number of variables in the data */
   )
{
   if( palim == -1 || palim > nvars - 1 )
      return nvars - 1;
   else
      return palim;
}

/** Read data in from a file, generate and store local scores.
 * In addition, information about e.g.\ scoring parameters as well as about the data is stored in a PropertyData data structure
 * @return SCIP_OKAY if all is well, otherwise some error code
 */
SCIP_RETCODE SC_readProblemInDataFormat(
   SCIP* scip,                  /**< SCIP data structure */
   const char* filename,        /**< File containing the data */
   int num_delims,              /**< The number of field delimiters to use */
   char* delims,                /**< The field delimiters to use */
   SCIP_Bool merge_delims,      /**< Whether multiple field delimiters should be merged in to one */
   ParentSetData* psd,          /**< Parent set data to populate */
   SCIP_Real*** scores,         /**< (*scores)[child][j] will be the score for the jth parent set of child */
   PropertyData* prop           /**< Property data structure  */
)
{

   char* score_type_str;
   int palim;
   int initpalim;

   PRIOR* prior;
   POSTERIOR* data;
   SCORING_EXTRAS* extras;
   SCORINGPARAMS* scoringparams;

   SCIP_CLOCK* clock;
   int nscores;

   /* create and start timing */
   SCIP_CALL( SCIPcreateClock(scip, &clock) );
   SCIP_CALL( SCIPstartClock(scip, clock) );

   SCIP_CALL( SCIPallocBlockMemory(scip, &prior) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &data) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &extras) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &scoringparams) );

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/continuous", &(scoringparams->continuous)) );
   if( scoringparams->continuous)
   {
#ifdef BLAS
      SCIP_CALL( readContinuousProblemFromFile(scip, filename, num_delims, delims,
            merge_delims, &(prior->nodeNames), &(data->continuous_data), &(prior->nvars), &(data->nrows)) );
      assert(data->continuous_data != NULL);
#else
      SCIPerrorMessage("You can't read in continuous data without BLAS support\n");
      return SCIP_ERROR;
#endif
   }
   else
   {
      SCIP_CALL( readProblemFromFile(scip, filename, num_delims, delims, merge_delims,
            &(prior->nodeNames), &(prior->arity), &(data->data), &(prior->nvars), &(data->nrows), &(prior->labels), &(data->num_observed)) );

      /* check that memory allocated as expected */
      assert(prior->arity != NULL);
      assert(prior->labels != NULL);
      assert(data->num_observed != NULL);
      assert(data->data != NULL);
   }

   assert(prior->nodeNames != NULL);

   SCIP_CALL( SCIPgetStringParam(scip, "gobnilp/scoring/score_type", &score_type_str) );
   if( strcmp(score_type_str,"BDeu") == 0 )
      scoringparams->score_type = SCORE_BDEU;
   else if( strcmp(score_type_str,"BGe") == 0 )
      scoringparams->score_type = SCORE_BGE;
   else if( strcmp(score_type_str,"BIC") == 0 )
      scoringparams->score_type = SCORE_BIC;
   else
   {
      SCIPerrorMessage("Unrecognised scoring function \"%s\". Use either \"BDeu\", \"BGe\" or \"BIC\"\n", score_type_str);
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/palim", &palim) );
   SCIP_CALL( SCIPgetIntParam(scip, "gobnilp/scoring/initpalim", &initpalim) );
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/pricing", &(scoringparams->pricing)) );

   scoringparams->palim = fixpalim(palim,prior->nvars);
   scoringparams->initpalim = fixpalim(initpalim,prior->nvars);

   if( scoringparams->pricing )
   {
      if( scoringparams->initpalim >= scoringparams->palim )
      {
         SCIPerrorMessage("\"gobnilp/scoring/initpalim\" (set to %d) must be strictly less than  \"gobnilp/scoring/palim\" (set to %d) when pricing is enabled", scoringparams->initpalim, scoringparams->palim);
         return SCIP_ERROR;
      }
   }
   else
   {
      /* suppress warning since there is no pricing on this branch! */
      /* if( scoringparams->initpalim != scoringparams->palim ) */
      /* { */
      /*    SCIPwarningMessage(scip, "Value of \"gobnilp/scoring/initpalim\" (set to %d) which differs from that \"gobnilp/scoring/palim\" (set to %d) being ignored since pricing is not enabled", scoringparams->initpalim, scoringparams->palim); */
      /* } */
      scoringparams->initpalim = scoringparams->palim;
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/scoring/prune", &(scoringparams->pruning)) );
   SCIP_CALL( SCIPgetRealParam(scip, "gobnilp/scoring/prunegap", &(scoringparams->prunegap)) );

   SCIP_CALL( setprior(scip, scoringparams, prior) );
   SCIP_CALL( setdata(scip, scoringparams, prior, data) );
   SCIP_CALL( setextras(scip, scoringparams, prior, data, extras) );

   if( scoringparams->pricing )
   {
      SCIP_PROBDATA* probdata;

      /* get problem data */
      probdata = SCIPgetProbData(scip);
      assert( probdata != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->data_etc->openlist), psd->n) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->data_etc->nopenlist), psd->n) );
   }

   SCIP_CALL( makelocalscores(scip, prior, data, extras, scoringparams, psd, scores, &nscores) );

   SCIP_CALL( setglobalproperties(scip, prop, prior, data, scoringparams, filename) );

   if( scoringparams->pricing )
   {
      SCIP_PROBDATA* probdata;

      /* get problem data */
      probdata = SCIPgetProbData(scip);
      assert( probdata != NULL );
      probdata->data_etc->prior = prior;
      probdata->data_etc->posterior = data;
      probdata->data_etc->scoring_extras = extras;
      probdata->data_etc->scoring_params = scoringparams;
   }
   else
   {
      freeextras(scip, scoringparams, extras);
      SCIP_CALL( freedata(scip, scoringparams, prior, data) );
      SCIP_CALL( freeprior(scip, scoringparams, prior) );

      SCIPfreeBlockMemory(scip, &scoringparams);
      SCIPfreeBlockMemory(scip, &extras);
      SCIPfreeBlockMemory(scip, &data);
      SCIPfreeBlockMemory(scip, &prior);
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%s: %d candidate parent sets created in %.1f seconds\n", filename, nscores, SCIPgetClockTime(scip, clock));

   SCIP_CALL( SCIPfreeClock(scip, &clock) );

   return SCIP_OKAY;

}
