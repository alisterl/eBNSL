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
 *  Function declarations for output.c
 */

#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <scip/scip.h>
#include "parent_set_data.h"
#include "property_data.h"
#include "model_averaging.h"


extern SCIP_RETCODE IO_addOutputParameters(SCIP* scip);
extern SCIP_RETCODE IO_printProblem(SCIP* scip, int run);
extern SCIP_RETCODE IO_doIterativePrint(SCIP* scip, MA_info* ma_info, ParentSetData* psd, int run);
extern SCIP_RETCODE IO_printParameters(SCIP* scip);
extern SCIP_RETCODE IO_printHeader(SCIP* scip);
extern SCIP_RETCODE IO_printScoresInJKLFormat(SCIP* scip, ParentSetData* psd, SCIP_Real** scores);
extern SCIP_RETCODE IO_printScoresInPSSFormat(SCIP* scip, ParentSetData* psd, PropertyData* prop);
extern SCIP_RETCODE IO_printcountsols(SCIP* scip, ParentSetData* psd, char* filename);

#endif
