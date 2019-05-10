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
 *  Defines a data structure that represents a vector with values of type double,
 *  as well as defining some vector operations.
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include "bge_vector.h"


/** Creates a vector with a given capacity.
 *  @param capacity The maximum capacity of the vector.
 *  @return An empty vector.
 */
Bge_Vector* BgeVectorCreate(
   int capacity
)
{
   Bge_Vector* vector = malloc(sizeof(Bge_Vector));
   vector->capacity = capacity;
   vector->items = malloc(capacity * sizeof(double));

   if (vector->items==NULL)
   {
     printf("error of memory allocation\n");
     exit(0);
   }

   return vector;
}


/** Deallocates the memory associated with a vector.
 *  @param vec The vector to free.
 */
void BgeVectorDelete(
   Bge_Vector** vec
)
{
  free((*vec)->items);
  free((*vec));
}

/** Computes the dot product between two vectors
 * @param vec_1 first input vector
 * @param vec_2 second input vector
 * @return output result of computing the dot product between vec_1 and vec_2
 */
double BgeVectorDotProduct(
   Bge_Vector* vec_1,
   Bge_Vector* vec_2
  )
{
  int index; /* Current index of each vector */
  double output = 0; /* the output of the inner product operation */

  for( index = 0; index < vec_1->capacity; index++ )
  {
      output += vec_1->items[index] * vec_2->items[index];
  }

  return output;
}

/** Performs The scalar multiplication operation for a vector
* @param scalar That the vector is multiplied by
* @param vec The input/output vector that is multiplied by the scalar
*/
void BgeVectorScalarMultiplication(
   double scalar,
   Bge_Vector* vec
  )
{
  int index;

  for(index = 0; index < vec->capacity; index++)
  {
    vec->items[index] = scalar * vec->items[index];
  }
}

/** Performs the vector subtraction operation
* @param vec_1 Input vector that has is being subtracted from
* @param vec_2 Input vector that is subtracted
* @param output_vec The result of the subtraction
*/
void BgeVectorSubtraction(
   Bge_Vector* vec_1,
   Bge_Vector* vec_2,
   Bge_Vector* output_vec
)
{
  int index;

  for(index = 0; index < vec_1->capacity; index ++)
  {
    output_vec->items[index] = vec_1->items[index] - vec_2->items[index];
  }

}
