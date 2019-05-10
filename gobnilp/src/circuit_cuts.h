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
 *  Function declarations for circuit_cuts.c
 */

#ifndef __SCIP_CIRCUIT_CUTS_H__
#define __SCIP_CIRCUIT_CUTS_H__

#include "scip/scip.h"
#include "stack.h"
#include "vector.h"
#include "vectorlist.h"
#include "parent_set_data.h"

/** A collection of storage items that can be allocated just once per constraint to save time. */
typedef struct
{
   SCIP_Real** adj_matrix;
   Stack* s;
   SCIP_Bool* blocked;
   VectorList* components;
   VectorList* cycles;
   VectorList* b;
   Vector* index_array;
   Vector* lowlink;
   SCIP_VAR** included_cluster;
   SCIP_VAR** excluded_cluster;
   SCIP_VAR** included_cycle;
   SCIP_VAR** excluded_cycle;
   int max_cycles;                        /**< The maximum number of cycles to find. */
   SCIP_Bool add_fractional_cuts;         /**< Whether fractional circuits should be looked for. */
   SCIP_Bool add_cluster_cuts;            /**< Whether cluster cuts should be added. */
   SCIP_Bool add_cycle_cuts;              /**< Whether cycle cuts should be added. */
   SCIP_Bool add_cluster_cuts_to_pool;    /**< Whether cluster cuts should be added to the cut pool. */
   SCIP_Bool add_cycle_cuts_to_pool;      /**< Whether cycle cuts should be added to the cut pool. */
   int max_cycle_length;                  /**< The maximum length of cycle to look for. */

} CircuitCutsStorage;

extern SCIP_RETCODE CC_addParams(SCIP* scip);
extern SCIP_RETCODE CC_initialise(SCIP* scip, CircuitCutsStorage* ccs, ParentSetData* data);
extern SCIP_RETCODE CC_finalise(SCIP* scip, CircuitCutsStorage* ccs, int n);
extern SCIP_RETCODE CC_findCuts(SCIP* scip, SCIP_CONSHDLR* conshdlr, ParentSetData* psd, SCIP_SOL* sol, CircuitCutsStorage* ccs, int* nGen, SCIP_Bool forcecuts, SCIP_Bool* found_efficacious_ptr, SCIP_Bool* cutoff);

#endif
