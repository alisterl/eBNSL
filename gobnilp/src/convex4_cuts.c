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
 *  Implements the functions needed to add acyclicity constraints based on
 *  the convex hull of the family-variable polytope for 4 nodes.
 */

#include "convex4_cuts.h"

#define DEFAULT_INITCONVEXHULL4BCUTS  0
#define DEFAULT_INITCONVEXHULL4CCUTS  0
#define DEFAULT_INITCONVEXHULL4DCUTS  0
#define DEFAULT_INITCONVEXHULL4ECUTS  0
#define DEFAULT_INITCONVEXHULL4GCUTS  0
#define DEFAULT_INITCONVEXHULL4HCUTS  0
#define DEFAULT_INITCONVEXHULL4ICUTS  0

#define DEFAULT_INITCONVEXHULL4BCUTS_ADDTOPOOL  FALSE
#define DEFAULT_INITCONVEXHULL4CCUTS_ADDTOPOOL  FALSE
#define DEFAULT_INITCONVEXHULL4DCUTS_ADDTOPOOL  FALSE
#define DEFAULT_INITCONVEXHULL4ECUTS_ADDTOPOOL  FALSE
#define DEFAULT_INITCONVEXHULL4GCUTS_ADDTOPOOL  FALSE
#define DEFAULT_INITCONVEXHULL4HCUTS_ADDTOPOOL  FALSE
#define DEFAULT_INITCONVEXHULL4ICUTS_ADDTOPOOL  FALSE

#define DEFAULT_INITCONVEXHULL4BCUTS_ADDCUT  TRUE
#define DEFAULT_INITCONVEXHULL4CCUTS_ADDCUT  TRUE
#define DEFAULT_INITCONVEXHULL4DCUTS_ADDCUT  TRUE
#define DEFAULT_INITCONVEXHULL4ECUTS_ADDCUT  TRUE
#define DEFAULT_INITCONVEXHULL4GCUTS_ADDCUT  TRUE
#define DEFAULT_INITCONVEXHULL4HCUTS_ADDCUT  TRUE
#define DEFAULT_INITCONVEXHULL4ICUTS_ADDCUT  TRUE

#define DEFAULT_CONVEXHULL4BCUTS  TRUE
#define DEFAULT_CONVEXHULL4CCUTS  FALSE
#define DEFAULT_CONVEXHULL4DCUTS  FALSE
#define DEFAULT_CONVEXHULL4ECUTS  FALSE
#define DEFAULT_CONVEXHULL4GCUTS  FALSE
#define DEFAULT_CONVEXHULL4HCUTS  FALSE
#define DEFAULT_CONVEXHULL4ICUTS  FALSE

#define DEFAULT_CONVEXHULL4BCUTS_ADDTOPOOL  TRUE
#define DEFAULT_CONVEXHULL4CCUTS_ADDTOPOOL  TRUE
#define DEFAULT_CONVEXHULL4DCUTS_ADDTOPOOL  TRUE
#define DEFAULT_CONVEXHULL4ECUTS_ADDTOPOOL  TRUE
#define DEFAULT_CONVEXHULL4GCUTS_ADDTOPOOL  TRUE
#define DEFAULT_CONVEXHULL4HCUTS_ADDTOPOOL  TRUE
#define DEFAULT_CONVEXHULL4ICUTS_ADDTOPOOL  TRUE

#define DEFAULT_CONVEXHULL4BCUTS_ADDCUT  TRUE
#define DEFAULT_CONVEXHULL4CCUTS_ADDCUT  TRUE
#define DEFAULT_CONVEXHULL4DCUTS_ADDCUT  TRUE
#define DEFAULT_CONVEXHULL4ECUTS_ADDCUT  TRUE
#define DEFAULT_CONVEXHULL4GCUTS_ADDCUT  TRUE
#define DEFAULT_CONVEXHULL4HCUTS_ADDCUT  TRUE
#define DEFAULT_CONVEXHULL4ICUTS_ADDCUT  TRUE

#define DEFAULT_CONVEXHULL4BCUTS_DELAY  FALSE
#define DEFAULT_CONVEXHULL4CCUTS_DELAY  TRUE
#define DEFAULT_CONVEXHULL4DCUTS_DELAY  TRUE
#define DEFAULT_CONVEXHULL4ECUTS_DELAY  TRUE
#define DEFAULT_CONVEXHULL4GCUTS_DELAY  TRUE
#define DEFAULT_CONVEXHULL4HCUTS_DELAY  TRUE
#define DEFAULT_CONVEXHULL4ICUTS_DELAY  TRUE


/** Add parameters controlling the addition of 'convex4' cuts to the SCIP instance */
SCIP_RETCODE C4_addParams(
   SCIP* scip              /**< SCIP data structure */
   )
{
   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/b_init_limit",
                             "maximum number of 'convexhull4b' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4BCUTS, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/c_init_limit",
                             "maximum number of 'convexhull4c' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4CCUTS, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/d_init_limit",
                             "maximum number of 'convexhull4d' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4DCUTS, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/e_init_limit",
                             "maximum number of 'convexhull4e' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4ECUTS, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/g_init_limit",
                             "maximum number of 'convexhull4g' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4GCUTS, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/h_init_limit",
                             "maximum number of 'convexhull4h' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4HCUTS, -1, INT_MAX, NULL, NULL));

   SCIP_CALL(SCIPaddIntParam(scip,
                             "constraints/convex4_cuts/i_init_limit",
                             "maximum number of 'convexhull4i' cuts to add initially (-1 for no upper bound)",
                             NULL, FALSE, DEFAULT_INITCONVEXHULL4ICUTS, -1, INT_MAX, NULL, NULL));




   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/b_init_addtopool",
                              "whether to add initial 'convexhull4b' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4BCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/c_init_addtopool",
                              "whether to add initial 'convexhull4c' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4CCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/d_init_addtopool",
                              "whether to add initial 'convexhull4d' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4DCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/e_init_addtopool",
                              "whether to add initial 'convexhull4e' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4ECUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/g_init_addtopool",
                              "whether to add initial 'convexhull4g' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4GCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/h_init_addtopool",
                              "whether to add initial 'convexhull4h' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4HCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/i_init_addtopool",
                              "whether to add initial 'convexhull4i' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4ICUTS_ADDTOPOOL, NULL, NULL));




   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/b_init_addcut",
                              "whether to add initial 'convexhull4b' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4BCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/c_init_addcut",
                              "whether to add initial 'convexhull4c' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4CCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/d_init_addcut",
                              "whether to add initial 'convexhull4d' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4DCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/e_init_addcut",
                              "whether to add initial 'convexhull4e' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4ECUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/g_init_addcut",
                              "whether to add initial 'convexhull4g' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4GCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/h_init_addcut",
                              "whether to add initial 'convexhull4h' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4HCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/i_init_addcut",
                              "whether to add initial 'convexhull4i' cuts to the initial LP",
                              NULL, FALSE, DEFAULT_INITCONVEXHULL4ICUTS_ADDCUT, NULL, NULL));




   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/b_enable",
                              "whether to use 'convexhull4b' cuts",
                              NULL, FALSE, DEFAULT_CONVEXHULL4BCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/c_enable",
                              "whether to use 'convexhull4c' cuts",
                              NULL, FALSE, DEFAULT_CONVEXHULL4CCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/d_enable",
                              "whether to use 'convexhull4d' cuts",
                              NULL, FALSE, DEFAULT_CONVEXHULL4DCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/e_enable",
                              "whether to use 'convexhull4e' cuts ",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ECUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/g_enable",
                              "whether to use 'convexhull4g' cuts",
                              NULL, FALSE, DEFAULT_CONVEXHULL4GCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/h_enable",
                              "whether to use 'convexhull4h' cuts ",
                              NULL, FALSE, DEFAULT_CONVEXHULL4HCUTS, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/i_enable",
                              "whether to use 'convexhull4i' cuts ",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ICUTS, NULL, NULL));




   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/b_addtopool",
                              "whether to add separating 'convexhull4b' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4BCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/c_addtopool",
                              "whether to add separating 'convexhull4c' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4CCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/d_addtopool",
                              "whether to add separating 'convexhull4d' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4DCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/e_addtopool",
                              "whether to add separating 'convexhull4e' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ECUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/g_addtopool",
                              "whether to add separating 'convexhull4g' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4GCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/h_addtopool",
                              "whether to add separating 'convexhull4h' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4HCUTS_ADDTOPOOL, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/i_addtopool",
                              "whether to add separating 'convexhull4i' cuts to the global cut pool",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ICUTS_ADDTOPOOL, NULL, NULL));




   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/b_addcut",
                              "whether to add separating 'convexhull4b' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4BCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/c_addcut",
                              "whether to add separating 'convexhull4c' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4CCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/d_addcut",
                              "whether to add separating 'convexhull4d' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4DCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/e_addcut",
                              "whether to add separating 'convexhull4e' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ECUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/g_addcut",
                              "whether to add separating 'convexhull4g' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4GCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/h_addcut",
                              "whether to add separating 'convexhull4h' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4HCUTS_ADDCUT, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/i_addcut",
                              "whether to add separating 'convexhull4i' cuts to the LP",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ICUTS_ADDCUT, NULL, NULL));




   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/b_delay",
                              "whether to only add 'convexhull4b' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4BCUTS_DELAY, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/c_delay",
                              "whether to only add 'convexhull4b' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4CCUTS_DELAY, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/d_delay",
                              "whether to only add 'convexhull4d' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4DCUTS_DELAY, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/e_delay",
                              "whether to only add 'convexhull4e' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ECUTS_DELAY, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/g_delay",
                              "whether to only add 'convexhull4g' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4GCUTS_DELAY, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/h_delay",
                              "whether to only add 'convexhull4h' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4HCUTS_DELAY, NULL, NULL));

   SCIP_CALL(SCIPaddBoolParam(scip,
                              "constraints/convex4_cuts/i_delay",
                              "whether to only add 'convexhull4i' when no cluster cuts can be found",
                              NULL, FALSE, DEFAULT_CONVEXHULL4ICUTS_DELAY, NULL, NULL));

   return SCIP_OKAY;

}


/** (conditionally) add a convex4 cut */
static
SCIP_RETCODE convexhull4x(
   SCIP*                  scip,                    /**< SCIP data structure */
   SCIP_CONSHDLR*         conshdlr,                /**< constraint handler */
   SolutionInfo*         solinfo,                /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_SOL*              sol,                     /**< solution to be separated */
   const int                    cluster[4],              /**< the 4 variables in the convex4 cut to be added e.g. [v0,v1,v2,v3] = [2,5,4,3] */
   const int                    convex4x,                /**< convex4 cut class, e.g 4b=1, 4c=2 */
   const int                    convex4x_ind,            /**< convex4 cut class, e.g 4b=1, 4c=2 but offset to pick out right pattern*/
   SCIP_Bool*             result,                  /**< to return  whether the cut was added */
   const SCIP_Bool              must_separate,           /**< whether to only add if it separates sol, otherwise added unconditionally */
   const SCIP_Bool              add_to_pool,             /**< whether to add cut to global cut pool */
   const SCIP_Bool              addcut,                  /**< whether to add cut to LP ( might just add to cut pool ) */
   const SCIP_Bool              forcecuts,               /**< whether to force cuts to be added */
   SCIP_Bool*             found_efficacious_ptr,    /**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*             cutoff
)
{
   char name[SCIP_MAXSTRLEN];
   SCIP_ROW* row;
   SCIP_VAR** vars;
   SCIP_Real* coeffs;
   int nvars = 0;
   int ci, cj;
   int i;
   int k;
   int size = 0;

   int coeff;
   unsigned int coeffi;
   unsigned int clausei;
   unsigned int liti;

   SCIP_Bool cnf_sat;
   SCIP_Bool clause_sat;
   SCIP_Bool ok;

   SCIP_Real lhs = 0;

   /* explanation of patterns:
      Each pattern is textually arranged below in 4 lines
      specifying when a parent set for cluster[0], cluster[1], cluster[2] and cluster[3]
      can be included in the cut. Each line specifies up to 3 clauses (a CNF). Each clause
      gives a condition which must be satisfied for a parent set variable to be included
      in the cut. For example, the first line of the 4B pattern is:
      { 3,-1,-1, 1, 2,-1,-1,-1,-1,
      3,-1,-1 states that 3 must be in the parent set of 0 for that parent set variable
      to be included (the -1's are just fillers)
      the next clause 1,2,-1 states that one of 1 and 2 must be in the parent set. there is no third clause which is indicated by 3 -1's

      if a cut class can have variables with coefficient 2 then patternshc2 is used otherwise patternshc1 is used
   */

   static const short int patternshc1[5][36] =
   {
      {
         1, 2, 3, -1, -1, -1, -1, -1, -1,
         0, 2, 3, -1, -1, -1, -1, -1, -1,
         0, 1, 3, -1 - 1, -1, -1, -1, -1,
         0, 1, 2, -1, -1, -1, -1, -1, -1
      },

      {
         3, -1, -1, 1, 2, -1, -1, -1, -1,
         0, 2, -1, 2, 3, -1, -1, -1, -1,
         0, 1, -1, 1, 3, -1, -1, -1, -1,
         0, -1, -1, 1, 2, -1, -1, -1, -1
      },

      {
         3, -1, -1, 1, 2, -1, -1, -1, -1,
         0, 2, -1, -1, -1, -1, -1, -1, -1,
         3, -1, -1, 0, 1, -1, -1, -1, -1,
         0, -1, -1, 2, -1, -1, -1, -1, -1
      },

      {
         2, -1, -1, 3, -1, -1, -1, -1, -1,
         2, 3, -1, -1, -1, -1, -1, -1, -1,
         0, -1, -1, 1, 3, -1, -1, -1, -1,
         0, -1, -1, 1, 2, -1, -1, -1, -1
      },

      {
         1, -1, -1, 2, -1, -1, 3, -1, -1,
         0, -1, -1, 2, -1, -1, 3, -1, -1,
         0, -1, 3, 1, -1, -1, 3, -1, -1,
         0, -1, -1, 1, -1, -1, 2, -1, -1
      }
   };

   static const short int patternshc2[5][72] =
   {
      {
         3, -1, -1, -1, -1, -1, -1, -1, -1,
         0, 2, 3, -1, -1, -1, -1, -1, -1,
         3, -1, -1, -1, -1, -1, -1, -1, -1,
         0, 1, 2, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         0, -1, -1, 2, -1, -1, -1, -1, -1
      },

      {
         1, 2, -1, 1, 3, -1, 2, 3, -1,
         0, -1, -1, 2, 3, -1, -1, -1, -1,
         0, -1, -1, 1, 3, -1, -1, -1, -1,
         0, -1, -1, 1, 2, -1, -1, -1, -1,
         1, -1, -1, 2, -1, -1, 3, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1
      },

      {
         1, 2, -1, -1, -1, -1, -1, -1, -1,
         0, 2, 3, -1, -1, -1, -1, -1, -1,
         0, -1, -1, -1, -1, -1, -1, -1, -1,
         0, -1, -1, 1, -1, -1, -1, -1, -1,
         1, -1, -1, 2, -1, -1, 3, -1, -1,
         0, -1, -1, 3, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1
      },

      {
         1, 2, -1, -1, -1, -1, -1, -1, -1,
         0, -1, -1, 3, -1, -1, -1, -1, -1,
         0, -1, -1, 3, -1, -1, -1, -1, -1,
         0, -1, -1, 1, -1, -1, 2, -1, -1,
         1, -1, -1, 2, -1, -1, 3, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1,
         -1, -1, -1, -1, -1, -1, -1, -1, -1
      },

      {
         1, 2, -1, -1, -1, -1, -1, -1, -1,
         0, 2, 3, -1, -1, -1, -1, -1, -1,
         0, 1, 3, -1, -1, -1, -1, -1, -1,
         1, 2, -1, -1, -1, -1, -1, -1, -1,
         1, -1, -1, 2, -1, -1, 3, -1, -1,
         0, -1, -1, 3, -1, -1, -1, -1, -1,
         0, -1, -1, 3, -1, -1, -1, -1, -1,
         0, -1, -1, 1, -1, -1, 2, -1, -1
      }
   };


   static const SCIP_Real rhss[10] = {3, 2, 2, 3, 2, 2, 3, 2, 4, 1};
   static const unsigned short int highest_coeffs[10] = {1, 1, 1, 2, 2, 1, 2, 2, 2, 1};
   static const char ids[] = "abcdefghij";

   const unsigned short int highest_coeff = highest_coeffs[convex4x];
   const SCIP_Real rhs = rhss[convex4x];
   const short int* pattern;
   const char id = ids[convex4x];

   SCIP_Bool* tmp_store;

   assert(highest_coeff == 1 || highest_coeff == 2);
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(solinfo != NULL);
   assert(result != NULL);
   assert(cluster != NULL);
   assert(found_efficacious_ptr != NULL);

   if( highest_coeff == 1 )
      pattern = patternshc1[convex4x_ind];
   else
      pattern = patternshc2[convex4x_ind];

   *result = FALSE;

   /* (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "convex%c(%d,%d,%d,%d)", convex4x->id, cluster[0], cluster[1], cluster[2], cluster[3]); */
   /* SCIP_CALL(  SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, -SCIPinfinity(scip), convex4x->rhs, FALSE, FALSE, TRUE )  ); */
   /* SCIP_CALL(  SCIPcacheRowExtensions(scip, row)  ); */

   for( ci = 0; ci < 4; ++ci )
      size = size + psd->nParentSets[cluster[ci]];

   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, size) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &coeffs, size) );


   /* find variables in cut for each child */
   for( ci = 0; ci < 4; ++ci )
   {

      /* shouldn't need this if cluster intelligently chosen */
      if( must_separate && !(lhs + highest_coeff * (4 - ci) > rhs) )
         break;

      i = cluster[ci];

      assert(i > -1);

      /* check each parent set for child i */
      for( k = 0; k < psd->nParentSets[i]; ++k )
      {

         tmp_store = store[i][k];

         /* each parent set must contain at least one of the other three,
            can just continue if not */
         ok = FALSE;
         for( cj = 0; cj < 4; ++cj )
         {
            if( cj == ci )
               continue;

            if( tmp_store[cluster[cj]] )
            {
               ok = TRUE;
               break;
            }
         }

         if( !ok )
            continue;

         /* look for parent sets with coeff=2 first */
         /* then with coeff=1  */
         for( coeff = highest_coeff; coeff > 0; --coeff )
         {
            assert(coeff == 2 || coeff == 1);
            /* coeffi = (9*ci)+(36*(coeff-1)); */
            /*assert( coeffi == 9*ci || coeffi == 36+(9*ci) ); */
            /* if ( coeff == 2 && convex4x->pattern[coeffi] == -1 ) */
            /*    continue; */
            if( coeff == 2 )
            {
               coeffi = 36 + (9 * ci);
               if( pattern[coeffi] == -1 )
                  continue;
            }
            else
               coeffi = 9 * ci;

            assert(coeffi < 72);

            cnf_sat = TRUE;
            /* check each clause is satisfied */
            for( clausei = coeffi; pattern[clausei] != -1 && clausei < coeffi + 9; clausei += 3 )
            {
               assert(clausei < 72);

               clause_sat = FALSE;
               /* search for a satisfied literal */

               for( liti = clausei; pattern[liti] != -1 && liti < clausei + 3; ++liti )
               {
                  assert(liti < 72);

                  /* if ( store[i][k][cluster[convex4x->pattern[liti]]] ) */
                  if( tmp_store[cluster[pattern[liti]]] )
                  {
                     /* found a satisfying literal, this clause satisfied */
                     clause_sat = TRUE;
                     break;
                  }
               }

               if( !clause_sat )
               {
                  /* clause broken so cnf unsatisfied
                     don't add parent set */
                  cnf_sat = FALSE;
                  break;
               }
            }

            if( cnf_sat )
            {
               vars[nvars] = psd->PaVars[i][k];
               coeffs[nvars++] = coeff;

               /* SCIP_CALL(  SCIPaddVarToRow(scip, row, solinfo->PaVars[i][k], coeff)  ); */
               lhs += coeff * solinfo->lpsolvals[i][k];
               /* never add a variable to cut twice, hence break */
               break;
            }
         }
      }
   }


   if( !must_separate || SCIPisGT(scip, lhs, rhs) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "convex%c(%d,%d,%d,%d)", id, cluster[0], cluster[1], cluster[2], cluster[3]);

      SCIP_CALL(SCIPcreateEmptyRowCons(scip, &row, conshdlr, name, -SCIPinfinity(scip), rhs,
                                       FALSE, FALSE, TRUE));
      SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, vars, coeffs) );

      /* SCIP_CALL(  SCIPflushRowExtensions(scip, row)  ); */
      SCIPdebugMessage(" -> Convex4cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                       name, SCIPgetRowLPActivity(scip, row), SCIProwGetRhs(row), SCIProwGetNorm(row),
                       SCIPgetCutEfficacy(scip, NULL, row),
                       SCIPgetRowMinCoef(scip, row), SCIPgetRowMaxCoef(scip, row),
                       SCIPgetRowMaxCoef(scip, row) / SCIPgetRowMinCoef(scip, row));
      SCIPdebug(SCIP_CALL( SCIPprintRow(scip, row, NULL) ));

      if( addcut )
      {

         SCIP_CALL( SCIPaddRow(scip, row, forcecuts, cutoff) );
         *result = TRUE;
      }
      if( add_to_pool )
         SCIP_CALL( SCIPaddPoolCut(scip, row) );

      if( SCIPisCutEfficacious(scip, sol, row) )
         *found_efficacious_ptr = TRUE;
   }
   else
   {
      SCIPdebugMessage("Failed: convex%c(%d,%d,%d,%d) since lhs=%f and rhs=%f\n", id, cluster[0], cluster[1], cluster[2], cluster[3], lhs, rhs);
   }

   /* SCIP_CALL(  SCIPreleaseRow(scip, &row)  ); */

   SCIPfreeMemoryArray(scip, &vars);
   SCIPfreeMemoryArray(scip, &coeffs);

   return SCIP_OKAY;
}


/**
  Search for 4B constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

  Example convex4b cut:

  \f$
  I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 2)+I(1 \leftarrow 0,2)+I(1 \leftarrow 0,3)+I(1 \leftarrow 2,3)+I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 1)+I(2 \leftarrow 0,1)+I(2 \leftarrow 0,3)+I(2 \leftarrow 1,3)+I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 0,1)+I(3 \leftarrow 0,2)+I(3 \leftarrow 0,1,2)
  \leq 2
  \f$

*/
static
SCIP_RETCODE cut_convexhull4b(
   SCIP*           scip,                  /**< SCIP data structure */
   SolutionInfo*  solinfo,                /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,              /**< constraint handler */
   int*            nGen,                  /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   SCIP_Bool       must_separate,         /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,           /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,                /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                 /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,             /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr, /**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0;
   int i1;
   int i2;
   int i3;
   int i;
   int j;
   int k;
   int ki;
   SCIP_Bool flag;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;

   for( i0 = 0; i0 < psd->n; ++i0 )
   {
      if( done )
         break;

      cluster[0] = i0;
      /* constraint not a cut unless i3 a parent of i0
         so consider all and only positive parents of i0 as potential choices for i3 */
      for( i = 0; i < solinfo->npa[i0]; ++i )
      {
         if( done )
            break;

         i3 = solinfo->pa[i0][i];

         /* swapping i0 and i3 does not alter constraint so only consider i3 > i0
            constraint not a cut unless i0 a parent of i3 */
         if( i3 < i0 || !solinfo->ispa[i3][i0] )
            continue;

         cluster[3] = i3;

         SCIPdebugMessage("Working on convex4b(%d,-,-,%d) since %d < %d and %d and %d are mutual positive parents\n", i0, i3, i0, i3, i0, i3);

         /* i2 must have a parent set containing both i0 and i3
            so at least i2 must be a child of i0 ... */

         for( j = 0; j < solinfo->nch[i0]; ++j )
         {

            if( done )
               break;

            i2 = solinfo->ch[i0][j];

            /* i2 must also be a child of i3 */
            if( !solinfo->ispa[i2][i3] )
               continue;

            assert(i2 != i0 && i2 != i3);

            /* now check there is at least one parent set for i2
               containing both i0 and i3 */

            flag = FALSE;
            for( ki = 0; ki < solinfo->nposvars[i2]; ++ki )
            {
               k = solinfo->posvars[i2][ki];
               if( store[i2][k][i0] && store[i2][k][i3] )
               {
                  flag = TRUE;
                  break;
               }
            }

            if( !flag )
               continue;

            cluster[2] = i2;

            SCIPdebugMessage("Working on convex4b(%d,-,%d,%d) since %d has a parent set containing both %d and %d\n", i0, i2, i3, i2, i0, i3);

            /* swapping i1 and i2 does not alter constraint so only consider i1 > i2 */

            for( i1 = i2 + 1; i1 < psd->n; ++i1 )
            {

               if( i1 == i0 || i1 == i3 )
                  continue;

               /* i1 must be a parent of one of the others */
               /* and some other stuff */
               if(!(
                        (solinfo->ispa[i0][i1] || solinfo->ispa[i2][i1] || solinfo->ispa[i3][i1])
                        && (solinfo->ispa[i0][i2] || solinfo->ispa[i1][i2] || solinfo->ispa[i3][i2])
                        && (solinfo->ispa[i1][i3] || solinfo->ispa[i2][i3])
                        && (solinfo->ispa[i1][i0] || solinfo->ispa[i2][i0])
                     )
                 )
                  continue;

               /* some prospect this is a cut .. */

               cluster[1] = i1;

               SCIPdebugMessage("Proposing convex4b(%d,%d,%d,%d) since: %d > %d; %d and %d both positive parents; %d is a positive parent of %d or %d; %d is a positive parent of %d or %d\n", i0, i1, i2, i3, i1, i2, i1, i2, i3, i1, i2, i0, i1, i2);
               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 1, 1,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if( *cutoff )
                  return SCIP_OKAY;

               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }
   SCIPdebugMessage("Found %d 4b cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}
/**
  Search for 4C constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

  Example 4c cut:

  \f$
  I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 0)+I(1 \leftarrow 2)+I(1 \leftarrow 0,2)+I(1 \leftarrow 0,3)+I(1 \leftarrow 2,3)+I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 0,3)+I(2 \leftarrow 1,3)+I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 0,2)+I(3 \leftarrow 0,1,2)
  \leq 2
  \f$

*/
static
SCIP_RETCODE cut_convexhull4c(
   SCIP*           scip,                 /**< SCIP data structure */
   SolutionInfo*  solinfo,             /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,             /**< constraint handler */
   int*            nGen,                 /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                  /**< solution to be separated */
   SCIP_Bool       must_separate,        /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,          /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,               /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,            /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,/**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0;
   int i1;
   int i2;
   int i3;
   int i;
   int j;
   int k;
   int ki;
   SCIP_Bool flag;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;

   for( i3 = 0; i3 < psd->n; ++i3 )
   {
      if( done )
         break;

      cluster[3] = i3;
      /* constraint not a cut unless i0 and i2 are parents of i3
      so consider all and only positive parent pairs of i3 as potential choices for i0 */
      for( i = 0; i < solinfo->npa[i3]; ++i )
      {
         if( done )
            break;

         i0 = solinfo->pa[i3][i];

         assert(i0 != i3);

         /* i3 must also be a parent of i0 */
         if( !(solinfo->ispa[i0][i3]) )
            continue;

         cluster[0] = i0;

         for( j = i + 1; j < solinfo->npa[i3]; ++j )
         {
            if( done )
               break;

            i2 = solinfo->pa[i3][j];

            assert(i2 != i3 && i2 != i0);

            /* think the following is a mistake */
            /* /\* i2 must also be a parent of i0 *\/ */
            /* if ( !(solinfo->ispa[i0][i2]) ) */
            /*    continue; */

            cluster[2] = i2;

            /* i3 must have a parent set containing both i0 and i2 */
            flag = FALSE;
            for( ki = 0; ki < solinfo->nposvars[i3]; ++ki )
            {
               k = solinfo->posvars[i3][ki];
               if( store[i3][k][i0] && store[i3][k][i2] )
               {
                  flag = TRUE;
                  break;
               }
            }

            if( !flag )
               continue;

            /* now just look for i1 */
            /* perhaps we can do better than just scanning? */
            for( i1 = 0; i1 < psd->n; ++i1 )
            {

               if( i1 == i0 || i1 == i2 || i1 == i3 )
                  continue;

               if( !(solinfo->ispa[i1][i0] || solinfo->ispa[i1][i2]) )
                  continue;

               /* some prospect this is a cut .. */

               cluster[1] = i1;

               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 2, 2,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if( *cutoff )
                  return SCIP_OKAY;
               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }
   SCIPdebugMessage("Found %d 4c cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}

/**
  Search for 4D constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

  \f$
  I(0 \leftarrow 3)+I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 0)+I(1 \leftarrow 2)+I(1 \leftarrow 3)+I(1 \leftarrow 0,2)+I(1 \leftarrow 0,3)+I(1 \leftarrow 2,3)+I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 3)+I(2 \leftarrow 0,3)+I(2 \leftarrow 1,3)+I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 0)+I(3 \leftarrow 1)+I(3 \leftarrow 2)+I(3 \leftarrow 0,1)+2I(3 \leftarrow 0,2)+I(3 \leftarrow 1,2)+2I(3 \leftarrow 0,1,2)
  \leq 3
  \f$

*/
static
SCIP_RETCODE cut_convexhull4d(
   SCIP*           scip,                 /**< SCIP data structure */
   SolutionInfo*  solinfo,             /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,             /**< constraint handler */
   int*            nGen,                 /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                  /**< solution to be separated */
   SCIP_Bool       must_separate,        /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,          /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,               /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,            /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,/**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0;
   int i1;
   int i2;
   int i3;
   int i;
   int j;
   int k;
   int ki;
   SCIP_Bool flag;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;

   for( i3 = 0; i3 < psd->n; ++i3 )
   {
      if( done )
         break;

      cluster[3] = i3;
      /* constraint not a cut unless i0 and i2 are parents of i3
         so consider all and only positive parent pairs of i3 as potential choices for i0 */
      for( i = 0; i < solinfo->npa[i3]; ++i )
      {
         if( done )
            break;

         i0 = solinfo->pa[i3][i];

         assert(i0 != i3);

         /* i3 must also be a parent of i0 */
         if( !(solinfo->ispa[i0][i3]) )
            continue;

         cluster[0] = i0;

         for( j = i + 1; j < solinfo->npa[i3]; ++j )
         {
            if( done )
               break;

            i2 = solinfo->pa[i3][j];

            assert(i2 != i3 && i2 != i0);

            /* i3 must also be a parent of i2 */
            if( !(solinfo->ispa[i2][i3]) )
               continue;

            cluster[2] = i2;

            /* i3 must have a parent set containing both i0 and i2 */
            flag = FALSE;
            for( ki = 0; ki < solinfo->nposvars[i3]; ++ki )
            {
               k = solinfo->posvars[i3][ki];
               if( store[i3][k][i0] && store[i3][k][i2] )
               {
                  flag = TRUE;
                  break;
               }
            }

            if( !flag )
               continue;

            /* now just look for i1 */
            /* perhaps we can do better than just scanning? */
            for( i1 = 0; i1 < psd->n; ++i1 )
            {

               if( i1 == i0 || i1 == i2 || i1 == i3 )
                  continue;

               /* some prospect this is a cut .. */

               cluster[1] = i1;

               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 3, 0,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if( *cutoff )
                  return SCIP_OKAY;
               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }
   SCIPdebugMessage("Found %d 4d cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}


/**
  Search for 4E constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

   Example 4e cut:

  \f$
  I(0 \leftarrow 1,2)+I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+2I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 0,2)+I(1 \leftarrow 0,3)+I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 0,1)+I(2 \leftarrow 0,3)+I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 0,1)+I(3 \leftarrow 0,2)+I(3 \leftarrow 0,1,2)
  \leq 2
  \f$
*/
static
SCIP_RETCODE cut_convexhull4e(
   SCIP*           scip,                  /**< SCIP data structure */
   SolutionInfo*  solinfo,              /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,              /**< constraint handler */
   int*            nGen,                  /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                   /**< solution to be separated */
   SCIP_Bool       must_separate,         /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,           /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,                /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                 /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,             /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr, /**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0, i1, i2, i3;
   int k, ki;
   int j_1, j2;
   SCIP_Bool flag;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;

   for( i0 = 0; i0 < psd->n; ++i0 )
   {

      if( done )
         break;

      /* i0 must have at least two parents */
      if( solinfo->npa[i0] < 2 )
         continue;

      cluster[0] = i0;

      /* due to symmetry breaking, can say i1 must be a parent of i0 */
      for( j_1 = 0; j_1 < solinfo->npa[i0]; ++j_1 )
      {
         if( done )
            break;

         i1 = solinfo->pa[i0][j_1];
         assert(i1 != i0);
         cluster[1] = i1;

         /* due to symmetry breaking, can say i2 must be a parent of i0 */
         for( j2 = j_1 + 1; j2 < solinfo->npa[i0]; ++j2 )
         {
            if( done )
               break;

            i2 = solinfo->pa[i0][j2];

            assert(i2 != i0  && i2 != i1);

            /* i0 must have a parent set containing both i1 and i2 */

            flag = FALSE;
            for( ki = 0; ki < solinfo->nposvars[i0]; ++ki )
            {
               k = solinfo->posvars[i0][ki];
               if( store[i0][k][i1] && store[i0][k][i2] )
               {
                  flag = TRUE;
                  break;
               }
            }

            if( !flag )
               continue;


            cluster[2] = i2;

            /* just try all possibilities for i3 */

            for( i3 = 0; i3 < psd->n; ++i3 )
            {

               if( i3 == i0  || i3 == i1 || i3 == i2 )
                  continue;

               if( !(solinfo->ispa[i0][i3] || solinfo->ispa[i1][i3] || solinfo->ispa[i2][i3]) )
                  continue;

               /* some prospect this is a cut .. */

               cluster[3] = i3;

               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 4, 1,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if (*cutoff )
                  return SCIP_OKAY;
               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }

            }
         }
      }
   }
   SCIPdebugMessage("Found %d 4e cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}


/**
  Search for 4G constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

  \f$
  I(0 \leftarrow 1)+I(0 \leftarrow 2)+I(0 \leftarrow 1,2)+I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+2I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 0)+I(1 \leftarrow 2)+I(1 \leftarrow 3)+I(1 \leftarrow 0,2)+2I(1 \leftarrow 0,3)+I(1 \leftarrow 2,3)+2I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 0)+I(2 \leftarrow 0,1)+I(2 \leftarrow 0,3)+I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 0,1)+I(3 \leftarrow 0,1,2)
\leq 3
\f$

*/
static
SCIP_RETCODE cut_convexhull4g(
   SCIP*           scip,                   /**< SCIP data structure */
   SolutionInfo*  solinfo,               /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,               /**< constraint handler */
   int*            nGen,                   /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                    /**< solution to be separated */
   SCIP_Bool       must_separate,          /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,            /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,                 /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                  /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,              /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,  /**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0, i1, i2, i3;
   int k, ki;
   int j_1, j_0, j2;
   SCIP_Bool flag;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;

   for( i3 = 0; i3 < psd->n; ++i3 )
   {
      /* both i0 and i1 must be parents of i3 */
      if( solinfo->npa[i3] < 2 )
         continue;

      cluster[3] = i3;

      for( j_1 = 0; j_1 < solinfo->npa[i3]; ++j_1 )
      {

         if( done )
            break;

         i1 = solinfo->pa[i3][j_1];
         assert(i1 != i3);

         /* i3 must be a parent of i1 */
         /* not sure that it must be so ..
         if ( !solinfo->ispa[i1][i3] )
            continue;
         */
         cluster[1] = i1;

         for( j_0 = 0; j_0 < solinfo->npa[i3]; ++j_0 )
         {
            if( done )
               break;

            /* no symmetry between i0 and i1 but can't be the same */
            if( j_1 == j_0 )
               continue;

            i0 = solinfo->pa[i3][j_0];
            assert(i0 != i3  && i0 != i1);

            /* i1 must be a parent of i0 */

            if( !solinfo->ispa[i0][i1] )
               continue;

            /* i3 must have at least one parent set containing both of i0 and i1 */

            flag = FALSE;
            for( ki = 0; ki < solinfo->nposvars[i3]; ++ki )
            {
               k = solinfo->posvars[i3][ki];
               if( store[i3][k][i0] && store[i3][k][i1] )
               {
                  flag = TRUE;
                  break;
               }
            }

            if( !flag )
               continue;

            cluster[0] = i0;

            /* i2 must be a child of i0 */

            for( j2 = 0; j2 < solinfo->nch[i0]; ++j2 )
            {
               i2 = solinfo->ch[i0][j2];

               assert(i2 != i0);
               if( i2 == i1 || i2 == i3 )
                  continue;

               /* some prospect this is a cut .. */

               cluster[2] = i2;

               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 6, 2,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if( *cutoff )
                  return SCIP_OKAY;

               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }
   SCIPdebugMessage("Found %d 4g cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}
/**
  Search for 4H constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

   Example 4h cut:

  \f$
  I(0 \leftarrow 1)+I(0 \leftarrow 2)+I(0 \leftarrow 1,2)+I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+2I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 0,3)+I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 0,3)+I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 0,1,2)
\leq 2
\f$

*/
static
SCIP_RETCODE cut_convexhull4h(
   SCIP*           scip,                   /**< SCIP data structure */
   SolutionInfo*  solinfo,               /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,               /**< constraint handler */
   int*            nGen,                   /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                    /**< solution to be separated */
   SCIP_Bool       must_separate,          /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,            /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,                 /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                  /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,              /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,  /**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0, i1, i2, i3;
   int k, ki;
   int j_1, j3, j2;
   SCIP_Bool flag;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;


   /* either have i0<-i1,i2,i3
      or i3<-i0,i1,i2
      treat two alternatives separately
   */


   for( i0 = 0; i0 < psd->n; ++i0 )
   {

      if( done )
         break;

      /* i1, i2 and i3 must be parents of i0 */
      if( solinfo->npa[i0] < 3 )
         continue;

      cluster[0] = i0;

      for( j_1 = 0; j_1 < solinfo->npa[i0]; ++j_1 )
      {

         if( done )
            break;

         i1 = solinfo->pa[i0][j_1];
         assert(i1 != i0);
         cluster[1] = i1;

         /* i1 and i2 are symmetric */
         for( j2 = j_1 + 1; j2 < solinfo->npa[i0]; ++j2 )
         {
            if( done )
               break;

            i2 = solinfo->pa[i0][j2];
            assert(i2 != i0 && i2 != i1);
            cluster[2] = i2;

            for( j3 = 0; j3 < solinfo->npa[i0]; ++j3 )
            {
               if( j3 == j_1 || j3 == j2 )
                  continue;

               i3 = solinfo->pa[i0][j3];
               assert(i3 != i0 && i3 != i1 && i3 != i2);

               /* i0 must have at least one parent set containing all of i1, i2 and i3 */

               flag = FALSE;
               for( ki = 0; ki < solinfo->nposvars[i0]; ++ki )
               {
                  k = solinfo->posvars[i0][ki];
                  if( store[i0][k][i1] && store[i0][k][i2] && store[i0][k][i3] )
                  {
                     flag = TRUE;
                     break;
                  }
               }

               if( !flag )
                  continue;


               /* some prospect this is a cut .. */

               cluster[3] = i3;

               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 7, 3,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if( *cutoff )
                  return SCIP_OKAY;
               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }

   SCIPdebugMessage("Found %d 4h cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}
/**
  Search for 4I constraints (derived from convex hull of DAG polytope with 4 nodes )
  Can be use to add initial cuts or cuts for a specific LP solution

  greedy search incomplete algorithm

   Example 4i cut:

  \f$
  I(0 \leftarrow 1)+I(0 \leftarrow 2)+I(0 \leftarrow 1,2)+I(0 \leftarrow 1,3)+I(0 \leftarrow 2,3)+2I(0 \leftarrow 1,2,3)+
  I(1 \leftarrow 0)+I(1 \leftarrow 2)+I(1 \leftarrow 3)+I(1 \leftarrow 0,2)+2I(1 \leftarrow 0,3)+I(1 \leftarrow 2,3)+2I(1 \leftarrow 0,2,3)+
  I(2 \leftarrow 0)+I(2 \leftarrow 1)+I(2 \leftarrow 3)+I(2 \leftarrow 0,1)+2I(2 \leftarrow 0,3)+I(2 \leftarrow 1,3)+2I(2 \leftarrow 0,1,3)+
  I(3 \leftarrow 1)+I(3 \leftarrow 2)+I(3 \leftarrow 0,1)+I(3 \leftarrow 0,2)+I(3 \leftarrow 1,2)+2I(3 \leftarrow 0,1,2)
  \leq 4
  \f$

*/

static
SCIP_RETCODE cut_convexhull4i(
   SCIP*           scip,                   /**< SCIP data structure */
   SolutionInfo*  solinfo,                 /**< solution information */
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR*  conshdlr,               /**< constraint handler */
   int*            nGen,                   /**< *nGen is how many cuts added */
   SCIP_SOL*       sol,                    /**< solution to be separated */
   SCIP_Bool       must_separate,          /**< if true cut is only added if it separates sol, otherwise added unconditionally */
   SCIP_Bool       add_to_pool,            /**< whether to add cut to global cut pool */
   SCIP_Bool       addcut,                 /**< whether to add cut ( might just add to cut pool ) */
   int             limit,                  /**< limit on number of cuts to add  (-1 for no limit ) */
   SCIP_Bool       forcecuts,              /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,  /**< to return whether an efficacious cutting plane was found */
   SCIP_Bool*      cutoff
)
{
   int i0, i1, i2, i3;
   int j_0, j3;
   int cluster[4];
   SCIP_Bool result;

   int ncuts = 0;
   SCIP_Bool done;

   if( limit == 0 )
      done = TRUE;
   else
      done = FALSE;

   /* say that i1 must have i0 and i3 as parents */

   for( i1 = 0; i1 < psd->n; ++i1 )
   {

      if( done )
         break;

      /* i0 and i3 must be parents of i1 */
      if( solinfo->npa[i1] < 2 )
         continue;

      cluster[1] = i1;

      for( j_0 = 0; j_0 < solinfo->npa[i1]; ++j_0 )
      {

         if( done )
            break;

         i0 = solinfo->pa[i1][j_0];
         assert(i0 != i1);
         cluster[0] = i0;

         /* say: i0 and i3 are symmetric */
         for( j3 = j_0 + 1; j3 < solinfo->npa[i1]; ++j3 )
         {
            if( done )
               break;

            i3 = solinfo->pa[i1][j3];
            assert(i3 != i0 && i3 != i1);
            cluster[3] = i3;

            /* say i0 and i3  must also be parents of i2 */

            for( i2 = 0; i2 < psd->n; ++i2 )
            {
               if( i2 == i0 || i2 == i1 || i2 == i3 )
                  continue;

               if( !(solinfo->ispa[i2][i0] || solinfo->ispa[i2][i3]) )
                  continue;

               /* some prospect this is a cut .. */

               cluster[2] = i2;

               SCIP_CALL( convexhull4x(scip, conshdlr, solinfo, psd, store, sol, cluster, 8, 4,
                     &result, must_separate, add_to_pool, addcut, forcecuts, found_efficacious_ptr, cutoff) );
               if( *cutoff )
                  return SCIP_OKAY;
               if( result )
               {
                  (*nGen)++;

                  if( ++ncuts == limit )
                  {
                     done = TRUE;
                     break;
                  }
               }
            }
         }
      }
   }

   SCIPdebugMessage("Found %d 4i cuts in separation method.\n", ncuts);
   return SCIP_OKAY;
}

/**
  Search for (potentially) all initial 'convex4' constraints (derived from convex hull of DAG polytope with 4 nodes )

  (Not all will be searched for; it depends on the parameter settings.)
*/
SCIP_RETCODE C4_add_initconvexhull_constraints(
   SCIP* scip,
   SCIP_CONSHDLR* conshdlr,
   SolutionInfo* solinfo,
   ParentSetData* psd,
   SCIP_Bool*** store
   )
{
   int nGen;
   SCIP_Bool found_efficacious;
   SCIP_Bool cutoff = FALSE;

   SCIP_Bool initconvexhull4xcuts_addtopool;
   SCIP_Bool initconvexhull4xcuts_addcut;
   int initconvexhull4xcuts;

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/b_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/b_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/b_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4b(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/c_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/c_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/c_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4c(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/d_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/d_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/d_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4d(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/e_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/e_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/e_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4e(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );
   
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/g_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/g_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/g_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4g(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/h_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/h_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/h_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4h(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/i_init_addtopool", &initconvexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/i_init_addcut", &initconvexhull4xcuts_addcut) );
   SCIP_CALL( SCIPgetIntParam(scip, "constraints/convex4_cuts/i_init_limit", &initconvexhull4xcuts) );
   SCIP_CALL( cut_convexhull4i(scip, solinfo, psd, store, conshdlr, &nGen, NULL, FALSE, initconvexhull4xcuts_addtopool, initconvexhull4xcuts_addcut, initconvexhull4xcuts, TRUE, &found_efficacious, &cutoff) );

   return SCIP_OKAY;
}

/**
  Search for (potentially) all non-initial 'convex4' constraints (derived from convex hull of DAG polytope with 4 nodes )

  (Not all will be searched for; it depends on the parameter settings.)
*/
SCIP_RETCODE C4_add_convexhull_constraints(
   SCIP* scip,
   SolutionInfo* solinfo,
   ParentSetData* psd,
   SCIP_Bool*** store,
   SCIP_CONSHDLR* conshdlr,
   SCIP_SOL* sol,
   int* nGen,
   SCIP_Bool forcecuts,
   SCIP_Bool* found_efficacious,
   SCIP_Bool* cutoff
   )
{

   SCIP_Bool convexhull4xcuts;
   SCIP_Bool convexhull4xcuts_delay;
   SCIP_Bool convexhull4xcuts_addtopool;
   SCIP_Bool convexhull4xcuts_addcut;

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/b_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/b_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/b_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/b_addcut", &convexhull4xcuts_addcut) );
   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4b(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/c_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/c_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/c_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/c_addcut", &convexhull4xcuts_addcut) );
   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4c(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
   }
   
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/d_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/d_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/d_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/d_addcut", &convexhull4xcuts_addcut) );
   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4d(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
   }
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/e_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/e_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/e_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/e_addcut", &convexhull4xcuts_addcut) );

   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4e(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
   }
   
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/g_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/g_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/g_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/g_addcut", &convexhull4xcuts_addcut) );
   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4g(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
   }
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/h_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/h_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/h_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/h_addcut", &convexhull4xcuts_addcut) );
   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4h(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
      if( *cutoff )
         return SCIP_OKAY;
   }
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/i_enable", &convexhull4xcuts) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/i_delay", &convexhull4xcuts_delay) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/i_addtopool", &convexhull4xcuts_addtopool) );
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/convex4_cuts/i_addcut", &convexhull4xcuts_addcut) );
   if( convexhull4xcuts && (!(*found_efficacious) || !(convexhull4xcuts_delay)) )
   {
      SCIP_CALL( cut_convexhull4i(scip, solinfo, psd, store, conshdlr, nGen, sol, TRUE,
            convexhull4xcuts_addtopool, convexhull4xcuts_addcut, -1, forcecuts, found_efficacious, cutoff) );
   }
   return SCIP_OKAY;
}
