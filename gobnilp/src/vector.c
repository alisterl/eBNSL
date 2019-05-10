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

/** @file
 *  Implements a data type representing a resizable array of integers.
 */

#include "vector.h"

/** Creates a vector with a given capacity.
 *  @param capacity The maximum capacity of the vector.
 *  @return An empty vector.
 */
Vector* VectorCreate(
   int capacity
   )
{
   Vector* vec = malloc(sizeof(Vector));
   vec->capacity = capacity;
   vec->size = 0;
   vec->items = malloc(sizeof(int) * capacity);
   return vec;
}

/** Appends an item to the end of a vector.
 *  @param vec The vector to extend.
 *  @param item The item to add.
 */
void VectorAppend(
   Vector* vec, 
   int item
   )
{
   vec->items[vec->size] = item;
   vec->size++;
}

/** Determine whether a vector contains a value.
 *  @param vec The vector to search.
 *  @param item The item to search for.
 *  @return Whether the item was in the vector.
 */
SCIP_Bool VectorContains(
   Vector* vec, 
   int item
   )
{
   int i;
   for( i = 0; i < vec->size; i++ )
      if( vec->items[i] == item )
         return TRUE;
   return FALSE;
}

/** Clears all the values from a vector.
 *  @param vec The vector to clear.
 */
void VectorClear(
   Vector* vec
   )
{
   vec->size = 0;
}

/** Deallocates the memory associated with a vector.
 *  @param vec The vector to free.
 */
void VectorDelete(
   Vector** vec
   )
{
   free((*vec)->items);
   free((*vec));
}


