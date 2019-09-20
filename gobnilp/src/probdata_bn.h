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
#include "scoring.h"

/** A k-cluster cut or cons. k=1 for 'normal' cluster cuts */
struct Cluster_Cut
{
   SCIP_CONS* cons;    /**< the constraint (or NULL if a cut) */
   SCIP_ROW* row;      /**< the cut (or NULL if a cons) */
   int nelts;          /**< size of the cluster cut */
   int* elts;          /**< BN variables in the cluster */
   int k;              /**< k value for the cluster */
   SCIP_Real dualsol;  /**< dual value for row/cons */
};
typedef struct Cluster_Cut CLUSTER_CUT;

   
/** Defines the BNSL problem instance 
 *
 * Currently incomplete since conss involving optional variables are not included
 */
struct SCIP_ProbData
{
   ParentSetData* psd;           /**< parent set data (includes family, arrow and edge vars) */
   SCIP_VAR*** ancestorvars;     /**< ancestor variables (or NULL if they do not exist) */
   SCIP_VAR**** im_vars;         /**< immorality indictator variables (or NULL if they do not exist) 
                                    im_vars[i][j][k]=1 iff there is an i->k<-j immorality */
   SCIP_VAR** posvars;           /**< BN variable position indicators (or NULL if they do not exist) */
   SCIP_VAR*** totalordervars;   /**< BN topological order variables  (or NULL if they do not exist) */
   SCIP_VAR** ispa;              /**< 'Being a parent' indicators (or NULL if they do not exist) */
   SCIP_VAR** pasize;            /**< pasize[i] is the size of the parent set for node i (or NULL if such variables do not exist) */
   SCIP_VAR*** posind;           /**< posind[i][j]==1 if variable i is in position j on 'the' total order (or NULL if such
                                    variables do not exist) */
   SCIP_VAR**** imsetvars;       /**< c-imset vars for subsets of size 3 (or NULL if they do not exist) */
   SCIP_VAR** nch;               /**< nch[i] is the number of children node i has (or NULL if such variables do not exist) */
   SCIP_VAR** isFemale;          /**< For pedigrees: indicators for which individuals are female (or NULL if they do not exist) */
   SCIP_VAR** arrow_vars;        /**< arrow_vars[i][j] is the indicator variable for i<-j, such variables also accessible via
                                    the hashmap in `psd` (or NULL if this array not used) */
   SCIP_CONS** one_parent_set;   /**< Constraints stating that each BN var has exactly one parent set */
   SCIP_CONS** arrow_conss;      /**< Constraints upper bounding sums of family variables by appropriate arrow variable */
   SCIP_CONS*  dagcluster_cons;  /**< the single dagcluster constraint */
   int nspc_conss;               /**< number of 'set packing' constraints */
   CLUSTER_CUT** spc_conss;      /**< the 'set packing' constraints */
   SCIP_Real** scores;           /**< scores[i][k] is the local score for the kth parent set for variable i */
   DATA_ETC* data_etc;           /**< everything needed to price in new family variables */
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
extern SCIP_RETCODE storeclustercut(SCIP* scip, CLUSTER_CUT** cluster_cut, SCIP_CONS* cons, SCIP_ROW* row, int nelts, int* elts, int k);
#endif
