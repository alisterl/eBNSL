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
 * Implements a data type representing a stack of integers.
 */

#include "stack.h"

/** Creates a new stack.
 *  @param capacity The maximum capacity of the stack.
 *  @return A new stack with with given capacity.
 */
Stack* StackCreate(
   int capacity
   )
{
   Stack* s = malloc(sizeof(Stack));
   s->max_capacity = capacity;
   s->num_items = 0;
   s->items = malloc(sizeof(int) * capacity);
   return s;
}

/** Gets the current number of items in a stack.
 *  @param s The stack to query.
 *  @return The number of items in the stack.
 */
int StackSize(
   Stack* s
   )
{
   return s->num_items;
}

/** Pops the top item from a stack.
 *  @param s The stack to use.
 *  @return The item that was on the top of the stack.
 */
int StackPop(
   Stack* s
   )
{
   int top = s->items[s->num_items - 1];
   s->num_items -= 1;
   return top;
}

/** Pushes an item onto the top of a stack.
 *  @param s The stack to add the item to.
 *  @param item The item to add to the stack.
 */
void StackPush(
   Stack* s,
   int item
   )
{
   s->items[s->num_items] = item;
   s->num_items += 1;
}

/** Determines whether a stack contains an item.
 *  @param s The stack to query.
 *  @param item The item to search for.
 *  @return TRUE if the item is found in the stack.  FALSE otherwise.
 */
SCIP_Bool StackContains(
   Stack* s,
   int item
   )
{
   int i;
   for( i = 0; i < s->num_items; i++ )
      if( s->items[i] == item )
         return TRUE;
   return FALSE;
}

/** Deallocates the memory used by a stack.
 *  @param s The stack to delete.
 */
void StackDelete(
   Stack** s
   )
{
   free((*s)->items);
   free((*s));
}

