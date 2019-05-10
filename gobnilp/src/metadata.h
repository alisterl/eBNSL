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
 *  Function declarations for metadata.c
 */

#ifndef __SCIP_METADATA_H__
#define __SCIP_METADATA_H__

#include "scip/scip.h"
#include "parent_set_data.h"
#include "pedigree_data.h"
#include "property_data.h"

extern SCIP_RETCODE MD_initialiseMetadata(SCIP* scip);

extern SCIP_RETCODE MD_setParentSetData(SCIP* scip, ParentSetData* psd);
extern ParentSetData* MD_getParentSetData(SCIP* scip);

extern SCIP_RETCODE MD_setPedigreeData(SCIP* scip, PedigreeData* pd);
extern PedigreeData* MD_getPedigreeData(SCIP* scip);

extern SCIP_RETCODE MD_setPropertyData(SCIP* scip, PropertyData* pd);
extern PropertyData* MD_getPropertyData(SCIP* scip);

#endif
