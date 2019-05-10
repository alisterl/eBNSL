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
 *  Provides a simple event handler to detect when the program should
 *  attempt to split dagcluster constraints into their individual
 *  strongly connected components.
 *
 *  The code for the actual splitting is in @link DC_tryToSplit() @endlink , this
 *  file just provides the code to ensure that that function is called at
 *  the correct time.
 */

#include "event_splitdag.h"
#include "cons_dagcluster.h"

/** The name of the event handler. */
#define EVENTHDLR_NAME "splitdag"
/** A description of the event handler. */
#define EVENTHDLR_DESC "event handler for splitting dagcluster constraints"

/** Adds the event catching once the problem is transformed. */
static SCIP_DECL_EVENTINIT(eventInitSplitDAG)
{
   /* Only want to listen to the event if gobnilp/splitdags = TRUE */
   SCIP_Bool splitdags;
   SCIP_CALL( SCIPgetBoolParam(scip, "gobnilp/splitdags", &splitdags) );
   if( splitdags )
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL) );
   return SCIP_OKAY;
}

/** Executes the code to try splitting the dagcluster constraint when the program enters a new node. */
static SCIP_DECL_EVENTEXEC(eventExecSplitDAG)
{
   /* Only want to do anything if this is a normal node and we are at the solving stage */
   if( SCIPnodeGetType(SCIPeventGetNode(event)) == SCIP_NODETYPE_FOCUSNODE && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      SCIP_CALL( DC_tryToSplit(scip) );
   return SCIP_OKAY;
}

/** Creates the event handler and registers it with the program.
 *  @param scip The SCIP instance to contain the handler.
 *  @return SCIP_OKAY if the handler was correctly created, or an error code otherwise.
 */
SCIP_RETCODE SCIPincludeEventHdlrSplitDAG(SCIP* scip)
{
   SCIP_EVENTHDLR* eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSplitDAG, NULL) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSplitDAG) );
   return SCIP_OKAY;
}
