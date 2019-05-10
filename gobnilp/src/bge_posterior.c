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
  * Describes a set of functions used for computing the posterior
  * matrix of a continuous data set
  */

 #include <stdio.h>
 #include <string.h>
 #include <assert.h>
 #include "bge_posterior.h"

 /** Finds the sample mean of the given data
 * @param data the sample data used to calculate the mean
 * @param mean_vec the output mean vector (capacity = data->cols)
*/
void SetSampleMean(
   Bge_Matrix* data,
   Bge_Vector* mean_vec
)
 {
   int row_index, col_index;
   double mean_entry = 0.0;

   // Compute the mean for each variable given the samples in the data set, then
   // insert them into the mean_vec
   for( col_index = 0; col_index < data->cols; col_index++ )
   {
     for( row_index = 0; row_index < data->rows; row_index++ )
     {
       mean_entry += data->items[(row_index * data->cols) + col_index];
     }

     mean_entry = mean_entry / data->rows;
     mean_vec->items[col_index] = mean_entry;
     mean_entry = 0.0;
   }
 }



 /** Finds the Sample Variance (multiplied by N-1) of a set of data
 * @param data the data set for the problem where rows = samples and cols = vars
 * @param mean_vec the vector containing the means for each variable
 * @param variance_matrix the output variance matrix  (dimensions = data->cols * data->cols)
*/
void SetSampleVariance(
   Bge_Matrix* data,
   Bge_Vector* mean_vec,
   Bge_Matrix* variance_matrix
)
{
   int sample_index;
   double* save;
   
   Bge_Vector* sample_vec   = BgeVectorCreate( data->cols ); // no of cols = no of variables
   Bge_Vector* diff_vec     = BgeVectorCreate( data->cols );
   Bge_Matrix* outer_matrix = BgeMatrixCreate( data->cols, data->cols );

   /* save the pointer to the space allocated to sample_vec->items */
   save = sample_vec->items;

   /* each sample is subtracted by the mean vector, then perform an
    * outer product with the resulting vector. This output matrix from the outer product
    * is then added to the variance_matrix */
   for( sample_index = 0; sample_index < data->rows; sample_index++ )  // for each sample in the data
   {
     sample_vec->items = &data->items[sample_index * data->cols]; // represents the current data sample (row vector)

     BgeVectorSubtraction(sample_vec, mean_vec, diff_vec); // subtract mean vector from data sample
     BgeVectorOuterProduct(diff_vec, diff_vec, outer_matrix);

     if( sample_index != 0 )
       BgeMatrixAddition(variance_matrix, outer_matrix, variance_matrix);
     else
       memcpy( variance_matrix->items, outer_matrix->items, sizeof(double) * data->cols * data->cols );
   }

   BgeMatrixDelete(&outer_matrix);
   BgeVectorDelete(&diff_vec);
   sample_vec->items = save;
   BgeVectorDelete(&sample_vec);
}

 /** Defines the Prior matrix for the Wishart disturbution over precision matrix W
 * @param n the number of variables in the gaussian network
 * @param alpha_mu hyper parameter for the wishart disturbution (normally set to 1)
 * @param alpha_omega hyper parameter for the wishart disturbution (normally set to n + 2)
 * @param prior_matrix prior matrix T for the wishart disturbtion over precision matrix W (dimensions = n * n)
*/
void SetPriorParametricMatrix(
   int n,
   int alpha_mu,
   int alpha_omega,
   Bge_Matrix* prior_matrix
   )
 {
   int i;
   int j;

   double t = (alpha_mu * (alpha_omega - n - 1.0)) / (alpha_mu + 1.0);

   for( i = 0; i < n; i++ )
   {
     for( j = 0; j < n; j++ )
     {
       if(i == j)
         prior_matrix->items[(n * i) + j] = t;
       else
        prior_matrix->items[(n * i) + j] = 0;
     }

   }
 }

/** Sets the posterior matrix of the normal wishart disturbution over precision matrix W
* and mean v
* @param data the problem data set
* @param the prior matrix for the wishart disturbution over the precision matrix W
* @param posterior_matrix the posterior matrix that is set for the normal wishart joint
* disturbution over the precision matrix and mean
*/
void SetPosteriorParametricMatrix(
   Bge_Matrix* data,
   Bge_Matrix* prior_matrix,
   Bge_Matrix* posterior_matrix,
   int alpha_mu,
   int alpha_omega
)
{

  int n         = data->cols; /* number of varibles in the data set */
  int obvs      = data->rows; /* the number of observations in the data set*/
  double scalar = (double) (obvs * alpha_mu) / (double) (obvs + alpha_mu);

  /* allocate memory to all the neccessary matrices/vectors for computing the posterior matrix */
  Bge_Vector* mean_vec        = BgeVectorCreate(n);
  Bge_Matrix* variance_matrix = BgeMatrixCreate(n,n);
  Bge_Matrix* tmp_matrix      = BgeMatrixCreate(n,n); /* the output of outerproduct(mean_vec, mean_vec) */

  assert(prior_matrix->cols == posterior_matrix->cols);
  assert(prior_matrix->rows == posterior_matrix->rows);

  
  /* initialize the mean vector, the variance matrix */
  SetSampleMean(data, mean_vec);
  SetSampleVariance(data, mean_vec, variance_matrix);
  SetPriorParametricMatrix(n, alpha_mu, alpha_omega, prior_matrix);

  /* computations to set the posterior_matrix */
  BgeVectorOuterProduct(mean_vec, mean_vec, tmp_matrix); /* assumes that v vector = [0,0...0] */
  //print_matrix_rowmajor( "temp matrix",tmp_matrix->rows, tmp_matrix->rows, tmp_matrix->items, tmp_matrix->rows );
  BgeMatrixScalarMultipliciation(scalar, tmp_matrix);

  //print_matrix_rowmajor( "temp matrix",tmp_matrix->rows, tmp_matrix->rows, tmp_matrix->items, tmp_matrix->rows );
  BgeMatrixAddition(variance_matrix, prior_matrix, posterior_matrix);
  BgeMatrixAddition(tmp_matrix, posterior_matrix, posterior_matrix);

  /* Free all memory allocated to matrices/vectors no longer needed */
  BgeMatrixDelete(&variance_matrix);
  BgeMatrixDelete(&tmp_matrix);
  BgeVectorDelete(&mean_vec);
}
