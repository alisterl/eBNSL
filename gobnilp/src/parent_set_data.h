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
 *  Function and type declarations for parent_set_data.c
 */

#ifndef __PARENT_SET_DATA_H__
#define __PARENT_SET_DATA_H__

#include <scip/scip.h>

/** The basic data needed to record a collection of parent sets associated with a problem. */
typedef struct
{
   int n;                   /**< Number of elements in the collection. */
   int* nParentSets;        /**< `nParentSets[i]` is the number of parent sets for element `i`. */
   int** nParents;          /**< `nParents[i][k]` is the number of parents in the `k`th parent set for element `i`. */
   int*** ParentSets;       /**< `ParentSets[i][k][l]` is the `l`th parent in the `k`th parent set of the `i`th element. */
   SCIP_VAR*** PaVars;      /**< `PaVars[i][k]` is the variable linked to the `k`th parent set of element `i`. */
   char** nodeNames;        /**< `nodeNames[i]` is the name of the `i`th node */
   SCIP_HASHMAP* arrow;     /**<  `arrow(key=n*i+j)` (if it exists) is an indicator variable for an arrow from j to i */
   SCIP_HASHMAP* edge;      /**<  `edge(key=n*i+j)`  (if it exists) is an indicator variable for an edge (in either direction) between i and j. Only non-NULL for j > i */
} ParentSetData;

extern SCIP_RETCODE PS_deallocateParentSetData(SCIP* scip, ParentSetData** psd, SCIP_Bool releasevars);
extern SCIP_RETCODE PS_copyParentSetData(SCIP* scip, ParentSetData* original, ParentSetData** duplicate);

extern SCIP_RETCODE PS_transformParentSetData(SCIP* scip, ParentSetData* original, ParentSetData** transformed);
extern SCIP_RETCODE PS_deepCopyParentSetData(SCIP* scip, SCIP* sourcescip, ParentSetData* original, ParentSetData** duplicate, SCIP_HASHMAP* consmap, SCIP_HASHMAP* varmap, SCIP_Bool global, SCIP_Bool* valid);


extern SCIP_RETCODE PS_splitToComponents(SCIP* scip, ParentSetData* original, int* num_components, ParentSetData*** components);
extern SCIP_RETCODE PS_specialiseFor(SCIP* scip, ParentSetData* original, int* nodes, int num_nodes, ParentSetData** specialisation);

extern SCIP_RETCODE PS_writeToFile(SCIP* scip, FILE* file, ParentSetData* psd);
extern SCIP_RETCODE PS_parse(SCIP* scip, char* str, ParentSetData** psd);

#endif
