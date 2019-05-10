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
 *  Function declarations for probdata_bn.c
 */

/*
This file was created by editing the file psd_lop.h that comes with the linear ordering example
in SCIP
*/

#ifndef __BN_probdata_bn__
#define __BN_probdata_bn__

#include "parent_set_data.h"
#include "model_averaging.h"

   
/** Defines the BNSL problem instance */
struct SCIP_ProbData
{
   ParentSetData* psd;         /**< parent set data (includes family, arrow and edge vars) */
   SCIP_Real** scores;         /**< local scores (obj coeff) for family variables */
   SCIP_VAR*** ancestorvars;   /**< ancestor variables (or NULL if they do not exist) */
   SCIP_VAR**** im_vars;       /**< immorality indictator variables (or NULL if they do not exist) 
                                  im_vars[i][j][k]=1 iff there is an i->k<-j immorality */
   SCIP_VAR** posvars;         /**< BN variable position indicators (or NULL if they do not exist) */
   SCIP_VAR*** totalordervars;         /**< BN topological order variables  (or NULL if they do not exist) */
   SCIP_VAR** ispa;            /**< 'Being a parent' indicators (or NULL if they do not exist) */
   SCIP_VAR** pasize;          /**< pasize[i] is the size of the parent set for node i (or NULL if such variables do not exist) */
   SCIP_VAR*** posind;         /**< posind[i][j]==1 if variable i is in position j on 'the' total order (or NULL if they do not exist) */
   SCIP_VAR**** imsetvars;     /**< c-imset vars for subsets of size 3 (or NULL if they do not exist) */
   SCIP_VAR** nch;             /**< nch[i] is the number of children node i has (or NULL if such variables do not exist) */
   SCIP_VAR** isFemale;        /**< For pedigrees: indicators for which individuals are female (or NULL if they do not exist) */
};

extern SCIP_RETCODE BN_setParamaterDefaults(SCIP* scip);
extern SCIP_RETCODE BN_suppresscols(SCIP* scip);
extern SCIP_RETCODE BN_printScores(SCIP* scip);
extern SCIP_RETCODE BN_printProblem(SCIP* scip, int run);
extern SCIP_RETCODE BN_doIterativePrint(SCIP* scip, MA_info* ma_info, int run);
extern SCIP_RETCODE BN_printParameters(SCIP* scip);
extern SCIP_RETCODE BN_printHeader(SCIP* scip);
extern SCIP_RETCODE BN_includePlugins(SCIP* scip);
extern SCIP_RETCODE BN_readProblem(SCIP* scip, char* inputformat, char* frequencyfile, const char* filename);
extern SCIP_RETCODE BN_addNonRepetitionConstraint(SCIP* scip, int run);
extern SCIP_RETCODE BN_addParameters(SCIP* scip);
extern          int BN_getNumberOfRepeats(SCIP* scip);
extern SCIP_RETCODE BN_printcountsols(SCIP* scip, char* filename);
#endif
