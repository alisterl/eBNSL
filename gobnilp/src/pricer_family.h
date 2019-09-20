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


/**@file   pricer_family.h
 * @ingroup PRICERS
 * @brief  family variable pricer
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRICER_FAMILY_H__
#define __SCIP_PRICER_FAMILY_H__


#include "scip/scip.h"
#include "cons_dagcluster.h"
#include "probdata_bn.h"
   
struct dualinfo
{
   int n;                           /**< number of BN vars in the problem */
   SCIP_Real* childpenalties;       /**< a dual penalty for each child */
   SCIP_Real* arrowpenalties;       /**< for looking up arrow dual penalties */
   int* nclusters;                  /**< nclusters[i] is the number of 'tight' cluster cuts/conss involving BN variable i */
   CLUSTER_CUT*** clusters;         /**< clusters[i][j] is the jth 'tight' cluster cut/cons involving BN variable i */
};
typedef struct dualinfo DUALINFO;

EXTERN
SCIP_RETCODE SCIPincludePricerFamily(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** added problem specific data to pricer and activates pricer */
extern
SCIP_RETCODE SCIPpricerFamilyActivate(
   SCIP*       scip,                     /**< SCIP data structure */
   int         n,                        /**< number of BN variables */
   SCIP_CONS** one_parent_set_conss,     /**< constraints stating that there is at most one parent set for any child */
   SCIP_CONS** arrow_conss,              /**< arrow_conss[n*i + j] is the cons where i<-j is an upper bound on sum of
                                            relevant family variables */
   SCIP_CONS*  dagcluster_cons,          /**< the single dagcluster constraint */
   int         nspc_conss,               /**< number of set packing constraints */
   CLUSTER_CUT** spc_conss              /**< 'set packing' constraints */
   );


#endif
