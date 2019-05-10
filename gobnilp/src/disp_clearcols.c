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
 *  Clearer column headings for GOBNILP
 */

#include "disp_clearcols.h"

#define DISP_NAME_SOFAR         "bestsofar"
#define DISP_DESC_SOFAR         "Score of the best network found so far"
#define DISP_HEADER_SOFAR       "Best Network Found So Far"
#define DISP_WIDTH_SOFAR        27
#define DISP_PRIORITY_SOFAR     0
#define DISP_POSITION_SOFAR     100
#define DISP_STRIPLINE_SOFAR    TRUE

#define DISP_NAME_POSSIBLE         "bestpossible"
#define DISP_DESC_POSSIBLE         "Score of the best possible network possible"
#define DISP_HEADER_POSSIBLE       "Best Network Possible"
#define DISP_WIDTH_POSSIBLE        27
#define DISP_PRIORITY_POSSIBLE     0
#define DISP_POSITION_POSSIBLE     200
#define DISP_STRIPLINE_POSSIBLE    TRUE

#define dispCopySoFar NULL
#define dispFreeSoFar NULL
#define dispInitSoFar NULL
#define dispExitSoFar NULL
#define dispInitsolSoFar NULL
#define dispExitsolSoFar NULL

#define dispCopyPossible NULL
#define dispFreePossible NULL
#define dispInitPossible NULL
#define dispExitPossible NULL
#define dispInitsolPossible NULL
#define dispExitsolPossible NULL

static SCIP_DECL_DISPOUTPUT(dispOutputSoFar)
{
   SCIP_Real bestsofar = SCIPgetPrimalbound(scip);
   if( SCIPisInfinity(scip, REALABS(bestsofar)) )
      SCIPinfoMessage(scip, file, "   No Solution Found Yet   ");
   else
      SCIPinfoMessage(scip, file, "       % 10.6e       ", bestsofar);
   return SCIP_OKAY;
}

static SCIP_DECL_DISPOUTPUT(dispOutputPossible)
{
   SCIP_Real bestpossible = SCIPgetDualbound(scip);
   if( SCIPisInfinity(scip, REALABS(bestpossible)) )
      SCIPinfoMessage(scip, file, "  No Upperbound Known Yet  ");
   else
      SCIPinfoMessage(scip, file, "       % 10.6e       ", bestpossible);
   return SCIP_OKAY;
}


SCIP_RETCODE SCIPincludeDispClearCols(SCIP* scip)
{
   SCIP_CALL(SCIPincludeDisp(scip, DISP_NAME_SOFAR, DISP_DESC_SOFAR, DISP_HEADER_SOFAR, SCIP_DISPSTATUS_AUTO,
                             dispCopySoFar, dispFreeSoFar, dispInitSoFar, dispExitSoFar, dispInitsolSoFar, dispExitsolSoFar, dispOutputSoFar,
                             NULL, DISP_WIDTH_SOFAR, DISP_PRIORITY_SOFAR, DISP_POSITION_SOFAR, DISP_STRIPLINE_SOFAR));
   SCIP_CALL(SCIPincludeDisp(scip, DISP_NAME_POSSIBLE, DISP_DESC_POSSIBLE, DISP_HEADER_POSSIBLE, SCIP_DISPSTATUS_AUTO,
                             dispCopyPossible, dispFreePossible, dispInitPossible, dispExitPossible, dispInitsolPossible, dispExitsolPossible, dispOutputPossible,
                             NULL, DISP_WIDTH_POSSIBLE, DISP_PRIORITY_POSSIBLE, DISP_POSITION_POSSIBLE, DISP_STRIPLINE_POSSIBLE));
   return SCIP_OKAY;
}
