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
 *  Function and type declarations for pedigree_data.c
 */

#ifndef __PEDIGREE_DATA_H__
#define __PEDIGREE_DATA_H__

#include <scip/scip.h>

/** Problem data that relates solely to pedigree-based constraints.  */
typedef struct
{
   /** The number of individuals with pedigree data. */
   int n;
   /** Variable recording whether each individual is female. */
   SCIP_VAR** SexVars;
   /** The ages of each individual. */
   int* ages;
   /** The initially supplied sex of each individual. (M, F or U) */
   char* sexes;
} PedigreeData;

extern SCIP_RETCODE PE_deallocatePedigreeData(SCIP* scip, PedigreeData** pd);
extern SCIP_RETCODE PE_copyPedigreeData(SCIP* scip, PedigreeData* original, PedigreeData** duplicate);

extern SCIP_RETCODE PE_writeToFile(SCIP* scip, FILE* file, PedigreeData* pd);
extern SCIP_RETCODE PE_parse(SCIP* scip, char* str, PedigreeData** pd);

#endif
