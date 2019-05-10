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

/** @file utils.h
 *  Function declarations for utils.c
 */

#ifndef __UTILS_H__
#define __UTILS_H__

#include <scip/scipdefplugins.h>
#include "parent_set_data.h"

extern SCIP_RETCODE hashtableCreateArrow(SCIP* scip, ParentSetData* psd);
extern SCIP_RETCODE hashtablefreeArrow(SCIP* scip, ParentSetData* psd);
extern SCIP_VAR* get_arrowedge(const ParentSetData* psd, const int i, const int j, SCIP_Bool arrow);
extern SCIP_VAR* get_arrow(const ParentSetData* psd, const int i, const int j);
extern SCIP_RETCODE put_arrow(SCIP* scip, ParentSetData* psd, const int i, const int j, SCIP_VAR* arrow_i_j);
extern SCIP_VAR* get_edge(const ParentSetData* psd, const int i, const int j);
extern SCIP_RETCODE put_edge(SCIP* scip, ParentSetData* psd, const int i, const int j, SCIP_VAR* edge_i_j);

extern char* SCIPstrdup(SCIP* scip, const char* str);

extern SCIP_RETCODE allsubsets(SCIP* scip, const int* set, const int n, int*** subsets_ptr);
extern SCIP_RETCODE free_allsubsets(SCIP* scip, int*** subsets_ptr);


/* Some convenient wrappers for creating new parameters that set many values to sensible defaults */
extern SCIP_RETCODE UT_addBoolParam(SCIP* scip, const char* name, const char* desc, SCIP_Bool value);
extern SCIP_RETCODE UT_addIntParam(SCIP* scip, const char* name, const char* desc, int value, int min, int max);
extern SCIP_RETCODE UT_addLongintParam(SCIP* scip, const char* name, const char* desc, SCIP_Longint value, SCIP_Longint min, SCIP_Longint max);
extern SCIP_RETCODE UT_addRealParam(SCIP* scip, const char* name, const char* desc, SCIP_Real value, SCIP_Real min, SCIP_Real max);
extern SCIP_RETCODE UT_addStringParam(SCIP* scip, const char* name, const char* desc, const char* value);

/* Some convenient wrappers for creating empty linear constraints that set many values to sensible defaults */
extern SCIP_RETCODE UT_createEmptyLinearConstraint(SCIP* scip, SCIP_CONS** cons, const char* name, SCIP_Real lhs, SCIP_Real rhs);
extern SCIP_RETCODE UT_createEmptyGTEConstraint(SCIP* scip, SCIP_CONS** cons, const char* name, SCIP_Real rhs);
extern SCIP_RETCODE UT_createEmptyLTEConstraint(SCIP* scip, SCIP_CONS** cons, const char* name, SCIP_Real lhs);

/* Functions and constants for writing and parsing items */
#define UT_LIST_START   '['
#define UT_LIST_END     ']'
#define UT_LIST_SEP     ','

extern SCIP_RETCODE UT_writeStringArray(SCIP* scip, FILE* file, char** input_array, int length);
extern SCIP_RETCODE UT_writeIntArray(SCIP* scip, FILE* file, int* input_array, int length);
extern SCIP_RETCODE UT_writeCharArray(SCIP* scip, FILE* file, char* input_array, int length);
extern SCIP_RETCODE UT_writeVarArray(SCIP* scip, FILE* file, SCIP_VAR** input_array, int length);
extern SCIP_RETCODE UT_writeIntArrayArray(SCIP* scip, FILE* file, int** input_array, int* lengths, int length);
extern SCIP_RETCODE UT_writeVarArrayArray(SCIP* scip, FILE* file, SCIP_VAR*** input_array, int* lengths, int length);
extern SCIP_RETCODE UT_writeIntArrayArrayArray(SCIP* scip, FILE* file, int*** input_array, int** lengths, int* lengths_of_lengths, int length);

extern SCIP_RETCODE UT_parseArray(char* input_string, char*** output_array);
extern SCIP_RETCODE UT_parseStringArray(char* input_string, char*** output_array, int length);
extern SCIP_RETCODE UT_parseIntArray(char* input_string, int** output_array, int length);
extern SCIP_RETCODE UT_parseCharArray(char* input_string, char** output_array, int length);
extern SCIP_RETCODE UT_parseVarArray(SCIP* scip, char* input_string, SCIP_VAR*** output_array, int length);
extern SCIP_RETCODE UT_parseIntArrayArray(char* input_string, int*** output_array, int* lengths, int length);
extern SCIP_RETCODE UT_parseVarArrayArray(SCIP* scip, char* input_string, SCIP_VAR**** output_array, int* lengths, int length);
extern SCIP_RETCODE UT_parseIntArrayArrayArray(char* input_string, int**** output_array, int** lengths, int* lengths_of_lengths, int length);

/* Functions for reading from a file */
extern SCIP_RETCODE UT_readFileAndSplit(SCIP* scip, FILE* file, char* splitters, int num_splitters, SCIP_Bool multiple_splitters_as_one, char**** lines, int* num_lines, int** line_lengths);

extern int get_index(char* nodeName, ParentSetData* psd);
extern int get_index_names(char* nodeName, int n, char** nodeNames);

extern SCIP_Bool is_dr_feasible(SCIP* scip, const ParentSetData* psd, SCIP_SOL* sol, SCIP_Bool constructing, SCIP_SOL* dr_sol);

#endif
