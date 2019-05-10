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
 *  Function and type declarations for property_data.c
 */

#ifndef __PROPERTY_DATA_H__
#define __PROPERTY_DATA_H__

#include <scip/scip.h>

/** Problem data that is not directly related to the parent sets. */
typedef struct
{

   /* Global Properties */
   /** The number of global properties. */
   int num_global;
   /** The names of the global properties. */
   char** global_property_names;
   /** The values of the global properties. */
   char** global_property_values;

   /* Local Properties */
   /** The number of individuals with properties. */
   int n;
   /** The number of properties of each individual. */
   int* num_properties;
   /** The names of the properties of each individual. */
   char*** property_names;
   /** The values of the properties of each individual. */
   char*** property_values;

} PropertyData;

extern SCIP_RETCODE PR_deallocatePropertyData(SCIP* scip, PropertyData** pd);
extern SCIP_RETCODE PR_copyPropertyData(SCIP* scip, PropertyData* original, PropertyData** duplicate);

extern SCIP_RETCODE PR_writeToFile(SCIP* scip, FILE* file, PropertyData* pd);
extern SCIP_RETCODE PR_parse(SCIP* scip, char* str, PropertyData** pd);

extern void PR_initialise(PropertyData* pd);

extern SCIP_Bool PR_hasGlobalProperty(SCIP* scip, PropertyData* pd, const char* name);
extern char* PR_getGlobalProperty(SCIP* scip, PropertyData* pd, const char* name);
extern SCIP_RETCODE PR_setGlobalProperty(SCIP* scip, PropertyData* pd, const char* name, const char* value);
extern SCIP_RETCODE PR_setGlobalPropertyFromInt(SCIP* scip, PropertyData* pd, const char* name, int value);
extern SCIP_RETCODE PR_setGlobalPropertyFromReal(SCIP* scip, PropertyData* pd, const char* name, SCIP_Real value);
extern SCIP_RETCODE PR_setGlobalPropertyFromBool(SCIP* scip, PropertyData* pd, const char* name, SCIP_Bool value);
extern SCIP_RETCODE PR_setGlobalPropertyFromArray(SCIP* scip, PropertyData* pd, const char* name, const char** value, int length);
extern SCIP_RETCODE PR_setGlobalPropertyFromRealArray(SCIP* scip, PropertyData* pd, const char* name, SCIP_Real* value, int length);

extern SCIP_Bool PR_hasProperty(SCIP* scip, PropertyData* pd, int individual, const char* name);
extern char* PR_getProperty(SCIP* scip, PropertyData* pd, int individual, const char* name);
extern SCIP_RETCODE PR_setProperty(SCIP* scip, PropertyData* pd, int individual, const char* name, const char* value);
extern SCIP_RETCODE PR_setPropertyFromInt(SCIP* scip, PropertyData* pd, int individual, const char* name, int value);
extern SCIP_RETCODE PR_setPropertyFromReal(SCIP* scip, PropertyData* pd, int individual, const char* name, SCIP_Real value);
extern SCIP_RETCODE PR_setPropertyFromBool(SCIP* scip, PropertyData* pd, int individual, const char* name, SCIP_Bool value);
extern SCIP_RETCODE PR_setPropertyFromArray(SCIP* scip, PropertyData* pd, int individual, const char* name, const char** value, int length);
extern SCIP_RETCODE PR_setPropertyFromRealArray(SCIP* scip, PropertyData* pd, int individual, const char* name, SCIP_Real* value, int length);

#endif
