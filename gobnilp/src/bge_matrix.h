/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2015 James Cussens, Mark Bartlett        *
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
 * Function declarations for bge_matrix.c
 */

#ifndef __BGE_MATRIX_H__
#define __BGE_MATRIX_H__

#include <lapacke.h>
#include "bge_vector.h"

/** A list of double values. */
typedef struct
{
   /* array of double containing the items within the matrix */
  double * items;

  /* number of rows within the matrix */
  lapack_int rows;
   /* number of cols within the matrix */
  lapack_int cols;
  /* defines the space in memory between rows within the matrix (Leading dimension)
   * this shoud be set to same value as col_dimension */
  lapack_int lda;

  /* defines the number of items in  matrix */
  int size;
} Bge_Matrix;

extern Bge_Matrix* BgeMatrixCreate(int cols, int rows);

extern void BgeMatrixDelete(Bge_Matrix** matrix);

extern double BgeMatrixLogDeterminant(Bge_Matrix* matrix);

extern void BgeMatrixSetPosteriorSubMatrix(unsigned int * family, int no_nodes, Bge_Matrix* posterior_matrix, Bge_Matrix* submatrix);

extern void BgeMatrixSetPosteriorParentMatrix(Bge_Matrix* sub_posterior, Bge_Matrix* parent_posterior);

extern void BgeMatrixScalarMultipliciation(double scalar, Bge_Matrix* matrix);

extern void BgeMatrixAddition(Bge_Matrix* matrix_1, Bge_Matrix* matrix_2, Bge_Matrix* output);

extern void BgeVectorOuterProduct(Bge_Vector* row_vector, Bge_Vector* col_vector, Bge_Matrix* output);

#endif
