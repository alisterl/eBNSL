/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2016 James Cussens, Mark Bartlett        *
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
  * Describes the main functions for scoring a node within a Gaussian Network
  * using the bge scoring metric.
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "bge_score.h"
#include "bge_matrix.h"

/** Computes The log of the bge score for a node either with or without parents
  * @param node The integer representation of the node in relation to the row/column in the dataset
  * @param family The family of the node being score where the value family[0] = node
  * @param no_parents The number of parents of the node
  * @param alpha_mu A hyper parameter used in computing the posterior_matrix
  * @param alpha_omega A hyper parameter used in computing the posterior_matrix
  * @param log_prefactor The ratio of logarithms for the score
  * @param log_gamma_ratio_table A table that contains all the possible gamma_ratios for each different size parent set
  * @param prior_matrix The prior matrix T
  * @param posterior_matrix The posterior matrix R
  * @param data - The continuous data that the network is being learned from
 */
double LogBgeScore(
   int node,
   unsigned int * family,
   int no_parents,
   int alpha_mu,
   int alpha_omega,
   double log_prefactor,
   double* log_gamma_ratio_table,
   Bge_Matrix* prior_matrix,
   Bge_Matrix* posterior_matrix,
   Bge_Matrix* data
)
{
   double score;

   if( no_parents == 0 )
   {
     score = LogBgeScoreWithoutParents(node, alpha_mu, alpha_omega, log_prefactor,
             log_gamma_ratio_table,prior_matrix, posterior_matrix, data);
   }
   else
   {
    score = LogBgeScoreWithParents(node, family, no_parents, alpha_mu, alpha_omega,
            log_prefactor, log_gamma_ratio_table,  prior_matrix, posterior_matrix, data);
   }

  return score;
}

/**
 *
 */
double log_prefactor(
   int no_samples,
   int alpha_mu
)
{
  /* Only needs to computed once for the entire search procedure */
  return 0.5 * (log(alpha_mu) - log(alpha_mu + no_samples));
}

double * create_log_gamma_ratio_table(
   int samples,
   int alpha_omega,
   int no_vars
)
{

   double * table;
   double log_gamma_ratio;

   int i;

   table = malloc(sizeof(double) * no_vars);

   for(i = 0; i < no_vars; i++)
   {
     log_gamma_ratio  = lgamma(0.5 * (samples + alpha_omega + - no_vars + i + 1));
     log_gamma_ratio -= lgamma(0.5 * (alpha_omega - no_vars + i + 1));
     log_gamma_ratio -= ((samples * 0.5) * log(M_PI));

     table[i] = log_gamma_ratio;
   }

   return table;
}


/** Computes The log of the bge score for a node without parents
  * @param node The integer representation of the node in relation to the row/column in the dataset
  * @param alpha_mu A hyper parameter used to in computing the posterior_matrix
  * @param alpha_omega A hyper parameter used to in computing the posterior_matrix
  * @param log_prefactor The ratio of logarithms for the score
  * @param log_gamma_ratio_table A table that contains all the possible gamma_ratios for each different size parent set
  * @param prior_matrix The prior matrix T
  * @param posterior_matrix The posterior matrix R
  * @param data - The continuous data that the network is being learned from
 */
double LogBgeScoreWithoutParents(
   int node,
   int alpha_mu,
   int alpha_omega,
   double log_prefactor,
   double* log_gamma_ratio_table,
   Bge_Matrix* prior_matrix,
   Bge_Matrix* posterior_matrix,
   Bge_Matrix* data
)
{

   double score;
   double sub_prior;
   double sub_posterior;

   int no_samples;
   int no_variables;

   no_samples   = data->rows; /* N number of samples */
   no_variables = data->cols; /* n number of variables */

  /* get sub prior_matrix and sub posterior_matrix (only a single element as there are no parents) */
   sub_prior     = prior_matrix->items[( node * no_variables ) + node];
   sub_posterior = posterior_matrix->items[( node * no_variables ) + node];

  /* compute the score */
   score  = log_prefactor;  /* log_prefactor is the same for every score */

   score += log_gamma_ratio_table[0]; /* size there are no parents get the first element in the table */

   score += (0.5 * (alpha_omega - no_variables + 1)) * log(sub_prior); /* determinant is equal to the single element
                                                                          in sub_prior/sub_posterior */
   score -= (0.5 * (no_samples + alpha_omega - no_variables + 1)) * log(sub_posterior);

   return score;
}

/** Computes The log of the bge score for a node with parents
  * @param node The integer representation of the node in relation to the row/column in the dataset
  * @param family The family of the node being score where the value family[0] = node
  * @param no_parents The number of parents of the node
  * @param alpha_mu A hyper parameter used to in computing the posterior_matrix
  * @param alpha_omega A hyper parameter used to in computing the posterior_matrix
  * @param log_prefactor The ratio of logarithms for the score
  * @param log_gamma_ratio_table A table that contains all the possible gamma_ratios for each different size parent set
  * @param prior_matrix The prior matrix T
  * @param posterior_matrix The posterior matrix R
  * @param data - The continuous data that the network is being learned from
 */
double LogBgeScoreWithParents(
   int node,
   unsigned int * family,
   int no_parents,
   int alpha_mu,
   int alpha_omega,
   double log_prefactor,
   double * log_gamma_ratio_table,
   Bge_Matrix* prior_matrix,
   Bge_Matrix* posterior_matrix,
   Bge_Matrix* data
)
{
  double score;

   int no_samples;
   int no_variables;

   Bge_Matrix* parent_posterior;
   Bge_Matrix *sub_posterior;

   no_samples   = data->rows; /* N number of total samples */
   no_variables = data->cols; /* n number of total variables */

   parent_posterior = BgeMatrixCreate(no_parents,no_parents);
   sub_posterior    = BgeMatrixCreate(no_parents + 1, no_parents + 1);

   /* set the sub_posterior matrix and parent_posterior matrix */
   BgeMatrixSetPosteriorSubMatrix(family, no_parents + 1, posterior_matrix, sub_posterior);
   BgeMatrixSetPosteriorParentMatrix(sub_posterior, parent_posterior);

   score  = log_prefactor; /* log_prefactor same for every score */

   score += log_gamma_ratio_table[no_parents]; /* get the gamma ration corresponding to the size of the parent set */

   score += (0.5 * (alpha_omega - no_variables + (2 * no_parents) + 1)) *
            log((alpha_mu * (alpha_omega - no_variables - 1.0)) / (alpha_mu + 1.0));

   score += ((0.5 * (no_samples + alpha_omega - no_variables + no_parents)))
            * BgeMatrixLogDeterminant(parent_posterior);

   score -= ((0.5 * (no_samples + alpha_omega - no_variables + no_parents + 1)))
            * BgeMatrixLogDeterminant(sub_posterior);

   BgeMatrixDelete(&parent_posterior);
   BgeMatrixDelete(&sub_posterior);

   return score;
}
