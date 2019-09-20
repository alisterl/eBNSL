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
 * Implements a Matrix Structure for floating point data and a set of matrix routines
 * for scoring a Gaussian Network
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> /* required for memcpy */
#include <assert.h>
#include <lapacke.h>

#include "bge_matrix.h"
#include "bge_vector.h"

/** Creates a new Matrix to be used by lapack routines
* @param col_dimension The number of columns of the matrix
* @param row_dimension The number of rows of the matrix
* @return A pointer to a matrix with dimensions row_dimension * col_dimension
*/
Bge_Matrix* BgeMatrixCreate(
   int cols,
   int rows
)
{
  Bge_Matrix* matrix = malloc(sizeof(Bge_Matrix));

  matrix->items = malloc(rows * cols * sizeof(double));
  matrix->rows = rows;
  matrix->cols = cols;
  matrix->lda = cols; /* sets the matrix layout to type LAPACK_ROW_MAJOR */

  if (matrix->items==NULL)
  {
    printf("error of memory allocation\n"); /* if memory allocation fails exit */
    exit(0);
  }
  return matrix;
}

/** Deletes a particular matrix and frees the memory allocated
* @param matrix The matrix to be deleted
*/
void BgeMatrixDelete(
  Bge_Matrix** matrix
  )
{
  free((*matrix)->items);
  free((*matrix));
}

/** Computes the log determinant of a matrix (Only possible for square matrices)
* @param Bge_Matrix** The matrix for which the log determinant is to be calculated
* @return The resulting value of the log determinant
*/
double BgeMatrixLogDeterminant(
   Bge_Matrix* matrix
)
{
  /* log to trace of the upper matrix after the PLU Factorisation */

   int index; /* Index used when iterating over factorized matrix */

   double log_det = 0; /* The output log determinant */

   int sign = 0; /* the sign of the determinant */

   double value = 0; /* current value being added in the determinant loop */

   lapack_int * ipiv; /* the pivot matrix that denotes the row interchanges made
                       * during the LAPACKE_dgetrf routine */

   lapack_int info; /* Variable that contains information after the execution of
                     * LAPACKE_dgetrf. If info > 0 then an error has occured when
                     * executing the LAPACKE_dgetrf routine */

   double * items_copy; /* since the LU factorisation is an in place operation then to
                         * maintain the existing matrix a copy of the items has to be
                         * made */

   if(matrix->cols == 1 && matrix->rows == 1)
     return log(matrix->items[0]);

   ipiv = malloc(matrix->rows * sizeof(lapack_int));
   items_copy = malloc(matrix->rows * matrix->cols * sizeof(double));

   /* copy the items to new array as LAPACKE_dgetrf is an in place operation */
   memcpy( items_copy, matrix->items, matrix->rows * matrix->cols * sizeof(double));

   /* Computes a LU factorisation making the determinant the sum of the resulting diagonal */
   info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, matrix->rows, matrix->cols, items_copy,
                         matrix->lda, ipiv);

  if(info > 0)
  {
     /* some sort of error has occurred! */
     return -INFINITY;
  }

   
   for( index = 0; index < matrix->rows; index++ )
   {
     value = items_copy[(matrix->lda + 1) * index];

      /* if value in ipiv matrix is not equal to index then change sign of log_det */
      if( ipiv[index] != index + 1 )
        sign += 1;

       if(value < 0)
       {
         sign += 1;
         value = value * -1;
       }

       log_det += log(value);
 }

   free(ipiv);
   free(items_copy);
   /* if determinant was negative return negative infinity */
  if(sign % 2 != 0)
  {
    return -INFINITY;
  }
  /* otherwise just return the log of the determinant */
  else
  {
     return log_det;
  }
}

/** Sets the posterior submatrix for a parent set of a node
* @param Bge_Matrix* the posterior matrix restricted to the node and its parents
* @param Bge_Matrix* the posterior matrix restricted to the parents (output)
*/
void BgeMatrixSetPosteriorParentMatrix(
   Bge_Matrix* sub_posterior,
   Bge_Matrix* parent_posterior
)
{
   int i; /* sub_posterior->cols * index to iterate over rows */
   int j; /* index to iterate over columns */
   int k; /* index to keep track of parent_posterior index */

   k = 0;

   /* originally was this: */
   /* for(i = 1; i < sub_posterior->rows; i++) */
   /* { */
   /*   for(j = 1; j < sub_posterior->cols; j++) */
   /*   { */
   /*     parent_posterior->items[k] = sub_posterior->items[(sub_posterior->cols * i) + j]; */
   /*     k++; */
   /*   } */
   /* } */

   
   for(i = sub_posterior->cols; i < sub_posterior->cols * sub_posterior->rows; i+= sub_posterior->cols)
     for(j = 1; j < sub_posterior->cols; j++)
       parent_posterior->items[k++] = sub_posterior->items[i + j];
}

/** Returns the Posterior matrix restricted to a family of a node
* @param family A set of nodes that contains a node and its parents
* @param no_nodes The number of nodes in the family
* @param Bge_Matrix* Posterior_matrix the posterior matrix
* @param Bge_Matrix* Submatrix the posterior matrix restricted to the variables in the family set
*/
void BgeMatrixSetPosteriorSubMatrix(
   unsigned int * family,
   int no_nodes,
   Bge_Matrix* posterior_matrix,
   Bge_Matrix* submatrix
)
{
   int i;
   int j;

   int matrix_entry; /* matrix entry index for the posterior matrix */
   int submatrix_entry; /* matrix entry for the sub posterior matrix */

   submatrix_entry = 0;

   for( i = 0; i < no_nodes; i++ )
   {
     for( j= 0; j < no_nodes; j++ )
     {
       matrix_entry = ( family[i] * posterior_matrix->cols ) + family[j];

      /* look through matrix and add the to relevant entries to the submatrix */
       submatrix->items[submatrix_entry] = posterior_matrix->items[matrix_entry];
       submatrix_entry++;
     }
   }
}

/** Multiples a given matrix be a scalar (Real number)
* @param scalar Each item in the matrix is mulitlplied by this number
* @param Bge_Matrix* The matrix that is to be multiplied by the scalar
*/
void BgeMatrixScalarMultipliciation(
   double scalar,
   Bge_Matrix* matrix
)
{
   int index; /* index to iterate over input matrix */

   for(index = 0; index < matrix->cols * matrix->rows; index++)
   {
     /* multiplies each value within the matrix by scalar */
     matrix->items[index] = scalar * matrix->items[index];
   }
}

/** Adds two matrices together
* @param matrix_1 matrix to be added
* @param matrix_2 matrix to be added
* @return output matrix for the addition of matrix_1 and matrix_2
*/
void BgeMatrixAddition(
   Bge_Matrix* matrix_1,
   Bge_Matrix* matrix_2,
   Bge_Matrix* output_matrix
  )
{
   int index;

   for( index = 0; index < matrix_1->cols * matrix_1->rows; index++ )
   {
     output_matrix->items[index] = matrix_1->items[index] + matrix_2->items[index];
   }
}

/** Takes two vectors and computes their outer product
 * @param row_vector First vector
 * @param col_vector Second vector
 * @return The outer product of vec_1 and vec_2 stored in a Bge_Matrix
 */
void BgeVectorOuterProduct(
   Bge_Vector* row_vector,
   Bge_Vector* col_vector,
   Bge_Matrix* output_matrix
)
{
   int col_index;
   int row_index;
   int output_matrix_index;

   for( col_index = 0; col_index < col_vector->capacity; col_index++ )
   {
     for( row_index = 0; row_index < row_vector->capacity; row_index++ )
     {
       output_matrix_index = (row_vector->capacity * col_index) + row_index;

       output_matrix->items[output_matrix_index] = col_vector->items[col_index] *
                                                   row_vector->items[row_index];
     }
   }
}
