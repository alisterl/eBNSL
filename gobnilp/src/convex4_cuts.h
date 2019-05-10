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
 *  Function declarations for convex4_cuts.c
 */

#ifndef __SCIP_CONVEX4_CUTS_H__
#define __SCIP_CONVEX4_CUTS_H__

#include "scip/scip.h"
#include "parent_set_data.h"
#include "solution_info.h"

extern SCIP_RETCODE C4_addParams(SCIP* scip);
extern SCIP_RETCODE C4_add_initconvexhull_constraints(SCIP* scip, SCIP_CONSHDLR*  conshdlr, SolutionInfo* solinfo, ParentSetData* psd, SCIP_Bool*** store);
extern SCIP_RETCODE C4_add_convexhull_constraints(SCIP* scip, SolutionInfo* solinfo, ParentSetData* psd, SCIP_Bool*** store, SCIP_CONSHDLR* conshdlr, SCIP_SOL* sol, int* nGen, SCIP_Bool forcecuts, SCIP_Bool* found_efficacious, SCIP_Bool* cutoff);

#endif
