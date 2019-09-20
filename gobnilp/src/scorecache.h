/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2019 James Cussens, Mark Bartlett        *
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

/** @file scorecache.h
 *  Function declarations for scorecache.c
 */

#ifndef __SCORECACHE_H__
#define __SCORECACHE_H__

#include "scip/scip.h"
#include "scoring.h"

struct scorecache
{
   int nvarscachelimit;           /**< subsets must have size below this to be cached. */
   int cachesizelimit;            /**< the maximum number of scores and count 
                                     values to cache (limit is common to both). */
   int cacheblocksize;            /**< how much to increase the size of the cache for scores 
                                     and counts when it is too small (common to both).  */
   int cachesize;                 /**< the current size of the cache for scores
                                     and counts (common to both) */
   int* nsubsets;                 /**< nsubsets[nvariables] is how many susbsets have size strictly less than nvariables */
   SCORE* scores;                 /**< scores[r] is the score for (the data projected onto ) 
                                     the unique subset of variables with rank r */
   COUNT* counts;                 /**< counts[r] is a count for (the data projected onto ) 
                                     the unique subset of variables with rank r 
                                     A value of zero for counts[r] indicates that no score (or count)
                                     has been previously stored for the subset with rank r */
                                  
};
typedef struct scorecache SCORECACHE;

extern void delete_scorecache(SCIP* scip, SCORECACHE* scorecache);
extern SCIP_RETCODE make_scorecache(SCIP* scip, SCORECACHE* scorecache, int nvars,
   int nvarscachelimit, int cachesizelimit,int cacheblocksize);
extern SCIP_RETCODE set_score_count_from_rank(SCIP* scip, SCORECACHE* scorecache, int rank,
   SCORE score, COUNT count);
extern SCORE get_score_count_rank(const SCORECACHE* scorecache, const VARIABLE* variables,
   int nvariables, int* rank, COUNT* count);

#endif
