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
 *  Function declarations for solution_info.c
 */

#ifndef __SCIP_SOLUTION_INFO_H__
#define __SCIP_SOLUTION_INFO_H__

#include "scip/scip.h"
#include "parent_set_data.h"

/** Stores information about a solution (typically an LP solution) */
typedef struct
{
   int**        pa;                   /**< pa[i][j] is the jth positive parent of i in the solution */
   int*         npa;                  /**< npa[i] is the number of positive parents of i in the solution */
   int**        ch;                   /**< ch[i][j] is the jth positive child of i in the solution */
   int*         nch;                  /**< nch[i] is the number of positive children of i in the solution */
   SCIP_Bool**  ispa;                 /**< ispa[i][j] = TRUE if j is in a positive parent set for i in the solution */
   int**        posvars;              /**< posvars[i][ki] is the (ki)th positive parent set for i in the solution */
   int*         nposvars;             /**< nposvars[i] is the number of positive parent sets for i in the soltion */
   SCIP_Real**  lpsolvals;            /**< stores values of variable in 'current' solution */
} SolutionInfo;

extern void SI_freesolinfo(SolutionInfo* solinfo, int n);

extern SCIP_RETCODE SI_setsolinfo(SCIP* scip, SolutionInfo* solinfo, ParentSetData* psd, SCIP_SOL* sol, SCIP_Bool augmented, SCIP_Bool dummysol);

#endif
