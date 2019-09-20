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

/** @file
 *  Functions for caching scores
 */

#include "scorecache.h"
#include <limits.h>
#include <stdio.h>

/** Map a subset of (BN) variables to a unique index.
 * If A is a subset of B then rank(A) < rank(B)
 * Note that the set of which A and B are subsets is not given explicitly as input
 * but it determines the value of `nsubsets`.
 * @return The rank of the subset
 */
static int rank_subset(
   const VARIABLE *variables,  /**< the subset of the variables */
   int nvariables,             /**< the size of the subset of the variables */
   int nsubsets                /**< nsubsets is how many subsets have size strictly less than nvariables */
)
{

   /* the first subset of size nvariables has rank nsubsets
    it is how many subsets have size strictly less than nvariables */
   int rank;
   int i;
   int j;
   int v_choose_iplusone;
   VARIABLE v;
   VARIABLE v2;
   VARIABLE v3;
   VARIABLE v4;

   assert(nvariables == 0 || variables != NULL);
   assert(nsubsets >= 0);

   rank = nsubsets;
   
   switch(nvariables)
   {
   case 0 :
      rank = 0;
      break;
   case 1 :
      rank = 1 + variables[0];
      break;
   case 2 :
      v2 = variables[1];
      assert(v2 > variables[0]);
      /* (C(v1,1) + C(v2,2)) */
      rank += variables[0] + v2 * (v2 - 1) / 2;
      break;
   case 3 :
      v2 = variables[1];
      v3 = variables[2];
      assert(v2 > variables[0]);
      assert(v3 > v2);
      /* (C(v1,1) + C(v2,2)) + C(v3,3) */
      rank += variables[0] + v2 * (v2 - 1) / 2  + v3 * (v3 - 1) * (v3 - 2) / 6;
      break;
   case 4 :
      v2 = variables[1];
      v3 = variables[2];
      v4 = variables[3];
      assert(v2 > variables[0]);
      assert(v3 > v2);
      assert(v4 > v3);
      /* (C(v1,1) + C(v2,2)) + C(v3,3) + C(v4,4) */
      rank += variables[0] + v2 * (v2 - 1) / 2  + v3 * (v3 - 1) * (v3 - 2) / 6  + v4 * (v4 - 1) * (v4 - 2) * (v4 - 3) / 24;
      break;
   default :
      for( i = 0; i < nvariables; ++i )
      {
         assert(i == 0 || variables[i] > v);
         v = variables[i];
         /* compute C(v,i+1)*/
         v_choose_iplusone = 1;
         for( j = 0; j < i + 1; ++j )
         {
            v_choose_iplusone *= (v - j);
            v_choose_iplusone /= (j + 1);
         }
         rank += v_choose_iplusone;
      }
      break;
   }

   assert(rank > -1);

   return rank;
}

/** Free the memory associated with a score cache **/
void delete_scorecache(
   SCIP* scip,              /**< SCIP instance */
   SCORECACHE* scorecache   /**< Score cache */
   )
{
   assert( scip !=  NULL);
   assert( scorecache !=  NULL);
   assert( scorecache->scores !=  NULL);
   assert( scorecache->counts !=  NULL);
   assert( scorecache->nsubsets !=  NULL);
   
   SCIPfreeMemoryArray(scip, &(scorecache->scores));
   SCIPfreeMemoryArray(scip, &(scorecache->counts));
   SCIPfreeMemoryArray(scip, &(scorecache->nsubsets));   

}

/** Initialise a score cache
 * @return SCIP_OKAY if all is well
 */
SCIP_RETCODE make_scorecache(
   SCIP* scip,                    /**< SCIP instance */
   SCORECACHE* scorecache,        /**< Score cache (to initialise) */
   int nvars,                     /**< size of set containing subsets */
   int nvarscachelimit,           /**< subsets must have size below this to be cached, 
                                     but not all subsets up to this limit will be cached if their
                                     rank (ie index) is too high */
   int cachesizelimit,            /**< the maximum number of scores and count 
                                     values to cache (limit is common to both). */
   int cacheblocksize             /**< how much to increase the size of the cache for scores 
                                     and counts when it is too small (common to both).  */
   )
{

   int accumulator;
   SCIP_Bool overflow;
   int i;
   int k;
   int nvars_choose_i;

   assert(scip != NULL);
   assert(scorecache != NULL);
   assert(nvars >= 0);
   assert(nvarscachelimit >= 0);
   assert(cachesizelimit >= 0);
   assert(cacheblocksize > 0);
   
   SCIP_CALL( SCIPallocMemoryArray(scip, &(scorecache->nsubsets), MIN(nvars + 1,nvarscachelimit)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(scorecache->scores), cacheblocksize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(scorecache->counts), cacheblocksize) );
   /* indicate that cache is empty by setting all counts to zero */
   for( i = 0; i < cacheblocksize; ++i )
      (scorecache->counts)[i] = 0;

   /* compute the number of subsets less than a given size, for all possible sizes */
   accumulator = 0;
   overflow = FALSE;
   for( i = 0; i <= MIN(nvars, nvarscachelimit - 1); ++i )
      {
         /* there are 'accumulator' subsets of cardinality less than i */
         (scorecache->nsubsets)[i] = accumulator;
         
         nvars_choose_i = 1;
         for( k = 0; k < i; ++k )
         {
            if( nvars_choose_i > INT_MAX / (nvars - k) )
            {
               overflow = TRUE;
               break;
            }
            nvars_choose_i *= (nvars - k);
            nvars_choose_i /= (k + 1);
         }
         if( INT_MAX - nvars_choose_i < accumulator )
            overflow = TRUE;
         else
            accumulator += nvars_choose_i;

         assert(accumulator > 0);

         /* accumulator is the highest rank of subsets of size i
            if accumulator is above INT_MAX then it is unsafe to compute
            ranks for subsets of this size, so set nvarscachelimit to
            stop this happening
         */
         if( overflow )
         {
            nvarscachelimit = i;
            break;
         }

         /* if accumulator has reached cachesizelimit then no point in ranking
            subsets of size greater than i
         */
         if( accumulator >= cachesizelimit )
         {
            nvarscachelimit = i+1;
            break;
         }
      }
   
   scorecache->nvarscachelimit = nvarscachelimit;
   scorecache->cachesizelimit = cachesizelimit;
   scorecache->cacheblocksize = cacheblocksize;
   scorecache->cachesize = cacheblocksize;

   return SCIP_OKAY;
}

/** Set the score and count associated with a subset using the rank of that subset
 * Any existing values associated with the given rank are silently over-written.
 * @return SCIP_OKAY if all is well, otherwise an error value
 */
SCIP_RETCODE set_score_count_from_rank(
   SCIP* scip,                   /**< SCIP instance */ 
   SCORECACHE* scorecache,       /**< scorecache */
   int rank,                     /**< rank (ie index for some subset) */
   SCORE score,                  /**< score associated with rank */
   COUNT count                   /**< count associated with rank */
)
{
   int i;

   assert(scip != NULL);
   assert(scorecache != NULL);
   assert(rank >= 0);

   if( rank < 0 )
   {
      SCIPerrorMessage("Rank is %d, but should not be negative\n", rank);
      return SCIP_ERROR;
   }

   if( count <= 0 )
   {
      SCIPerrorMessage("Count is %d, but should be positive\n", count);
      return SCIP_ERROR;
   }

   if( rank >= scorecache->cachesizelimit )
   {
      SCIPerrorMessage("Rank %d exceeds allowed size %d of scorecache\n", rank, scorecache->cachesizelimit);
      return SCIP_ERROR;
   }

   while( rank >= scorecache->cachesize )
   {
      scorecache->cachesize += scorecache->cacheblocksize;
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(scorecache->scores), scorecache->cachesize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(scorecache->counts), scorecache->cachesize) );

      /* flag newly allocated space as empty by setting count to 0, no need to initialise score */
      for( i = scorecache->cachesize - scorecache->cacheblocksize; i < scorecache->cachesize; ++i )
         (scorecache->counts)[i] = 0;
   }

   (scorecache->scores)[rank] = score;
   (scorecache->counts)[rank] = count;

   return SCIP_OKAY;
}

/** Return score, count and rank for the given subset, the latter two via pointers
 *  If `*rank==-1` then the subset has too high a rank to have a score/count stored for it.
 *  If `*count==0` then the score and count for this subset have yet to be computed
 * @return The score associated with the subset which has the given rank or 0.0 if there is no such subset.
 */
SCORE get_score_count_rank(
   const SCORECACHE* scorecache, /**< scorecache */
   const VARIABLE *variables,    /**< the subset of the variables */
   int nvariables,               /**< the size of the subset of the variables */
   int* rank,                    /**< rank (ie index for the subset), -1 indicates that 
                                    the actual rank would be too large */
   COUNT* count                  /**< count associated with rank, 0 indicates no score or count
                                    has yet been computed for the subset with this rank */
)
{
   int actualrank;

   assert(scorecache != NULL);
   assert(scorecache->nsubsets != NULL);
   assert(nvariables == 0 || variables != NULL);
   assert(rank != NULL);
   assert(count != NULL);
   assert(nvariables >= 0);
   
   /* 'default' values */
   *rank = -1;
   *count = 0;
   
   if( nvariables < scorecache->nvarscachelimit )
   {
      actualrank = rank_subset(variables, nvariables, (scorecache->nsubsets)[nvariables]);
      
      if( actualrank < scorecache->cachesizelimit )
         *rank = actualrank;
   }
   /* if rank too high, return dummy value for score */
   if( *rank == -1 )
      return 0.0;

   assert(*rank == actualrank);
   
   if( actualrank < scorecache->cachesize )
   {
      *count = (scorecache->counts)[actualrank];
      return (scorecache->scores)[actualrank];
   }
   else
   {
      /* if score is missing, return dummy value for score */
      return 0.0;
   }
}


