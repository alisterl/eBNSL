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
 *  Function declarations for fractional_circuit_cuts.c
 */

#ifndef __SCIP_FRACTIONAL_CIRCUIT_CUTS_H__
#define __SCIP_FRACTIONAL_CIRCUIT_CUTS_H__

#include "scip/scip.h"
#include "circuit_cuts.h"

extern SCIP_RETCODE FC_addParams(SCIP* scip);
/* extern SCIP_RETCODE FC_initialise(SCIP* scip, CircuitCutsStorage* ccs, ParentSetData* data); */
/* extern SCIP_RETCODE FC_finalise(SCIP* scip, CircuitCutsStorage* ccs, int n); */
extern SCIP_RETCODE FC_findCuts(SCIP* scip, SCIP_CONSHDLR* conshdlr, ParentSetData* psd, SCIP_SOL* sol, CircuitCutsStorage* ccs, int* nGen, SCIP_Bool forcecuts, SCIP_Bool* found_efficacious_ptr, SCIP_Bool* cutoff);

#endif

