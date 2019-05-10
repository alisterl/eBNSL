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
 *  Implements a data type representing a resizable array of resizable arrays of integers.
 */

#include "vectorlist.h"

/** Creates a vectorlist with a given capacity.
 *  @param capacity The maximum capacity of the vectorlist.
 *  @return An empty vectorlist.
 */
VectorList* VectorListCreate(
   int capacity
   )
{
   VectorList* vl = malloc(sizeof(VectorList));
   vl->capacity = 0;
   vl->size = 0;
   vl->items = malloc(sizeof(Vector**)*capacity);
   return vl;
}

/** Appends an item to the end of a vectorlist.
 *  @param vl The vectorlist to extend.
 *  @param item The item to add.
 */
void VectorListAppend(
   VectorList* vl, 
   Vector* item
   )
{
   vl->items[vl->size] = item;
   vl->size++;
}

/** Clears all the values from a vectorlist.
 *  @param vl The vectorlist to clear.
 */
void VectorListClear(
   VectorList* vl
   )
{
   int i;
   for( i = 0; i < vl->size; i++ )
      VectorDelete(&(vl->items[i]));
   vl->size = 0;
}

/** Deallocates the memory associated with a vectorlist.
 *  @param vl The vectorlist to free.
 */
void VectorListDelete(
   VectorList** vl
   )
{
   int i;
   for( i = 0; i < (*vl)->size; i++ )
      VectorDelete(&((*vl)->items[i]));
   free(((*vl)->items));
   free((*vl));
}

