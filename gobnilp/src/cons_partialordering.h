/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_partialordering.h
 * @brief  constraint handler for partial ordering constraints
 * @author Marc Pfetsch
 *
 * This constraint a partial order
 * as follows. The rows and columns correspond to the elements of the set to be ordered. A variable
 * x[i][j] is 1 if and only if element i appears before j in the order.
 *
 * In this constraint handler we only add the symmetry inequalities and separate the triangle
 * inequalities yielding a correct IP model.
 *
 * The variables on the diagonal are ignored.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __LOP_CONS_PARTIALORDERING_H__
#define __LOP_CONS_PARTIALORDERING_H__


#include <scip/scip.h>
#include "parent_set_data.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for partial ordering constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrPartialOrdering(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a partial ordering constraint */
extern
SCIP_RETCODE SCIPcreateConsPartialOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of binary variables */
   ParentSetData*        psd,                /**< contains arrow and parent set variables */
   SCIP_Bool             minimal,            /**< whether the partial order must minimal while being consistent with arrows */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,               /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   );


/** creates and captures a partial ordering constraint with all its constraint flags set to their
 *  default values
 *
 */
extern
SCIP_RETCODE SCIPcreateConsBasicPartialOrdering(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   n,                  /**< number of elements */
   SCIP_VAR***           vars,               /**< n x n matrix of binary variables */
   ParentSetData*        psd,                /**< contains arrow and parent set variables */
   SCIP_Bool             minimal             /**< whether the partial order must minimal while being consistent with arrows */
   );

#ifdef __cplusplus
}
#endif

#endif
