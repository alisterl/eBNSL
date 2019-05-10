
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

/**@file   cons_ci.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for ci constraints
 * @author James Cussens
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CI_H__
#define __SCIP_CONS_CI_H__


#include "scip/scip.h"
#include "parent_set_data.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for ci constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrCi(
   SCIP*                 scip                /**< SCIP data structure */
);

/** creates and captures a ci constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
extern
SCIP_RETCODE SCIPcreateConsCi(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   ParentSetData*        psd,                /**< parent set data containing family, edge and arrow variables */
   SCIP_VAR***           ancestorvars,       /**< anectorvars (or NULL if absent) */
   int*                  a,                  /* the set a */
   int                   n_a,                /* size of set a */
   int*                  b,                  /* the set b */
   int                   n_b,                /* size of set b */
   int*                  s,                  /* the set s */
   int                   n_s,                /* size of set s */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
);

/** creates and captures a ci constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
/* extern */
/* SCIP_RETCODE SCIPcreateConsBasicCi( */
/*    SCIP*                 scip,               /\**< SCIP data structure *\/ */
/*    SCIP_CONS**           cons,               /\**< pointer to hold the created constraint *\/ */
/*    const char*           name,               /\**< name of constraint *\/ */
/*    int                   nvars,              /\**< number of variables in the constraint *\/ */
/*    SCIP_VAR**            vars,               /\**< array with variables of constraint entries *\/ */
/*    SCIP_Real*            coefs,              /\**< array with coefficients of constraint entries *\/ */
/*    SCIP_Real             lhs,                /\**< left hand side of constraint *\/ */
/*    SCIP_Real             rhs                 /\**< right hand side of constraint *\/ */
/*    ); */

#ifdef __cplusplus
}
#endif

#endif
