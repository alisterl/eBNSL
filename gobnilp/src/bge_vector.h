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
 * Function declarations for bge_vector.c
 */

#ifndef __BGE_VECTOR_H__
#define __BGE_VECTOR_H__

/** A list of double values. */
typedef struct
{
   /** The maximum capacity of the vector. */
   int  capacity;
   /** The items currently in the vector. */
   double * items;

} Bge_Vector;

extern Bge_Vector* BgeVectorCreate(int capacity);

extern void BgeVectorDelete(Bge_Vector** vec);

extern double BgeVectorDotProduct(Bge_Vector* vec_1, Bge_Vector* vec_2);

extern void BgeVectorScalarMultiplication(double scalar, Bge_Vector* vec);

extern void BgeVectorSubtraction(Bge_Vector* vec_1, Bge_Vector* vec_2, Bge_Vector* output_vec);

#endif
