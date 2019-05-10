/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2015 James Cussens, Mark Bartlett        *
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
 * Function declarations for bge_score.c
 */

#include "bge_matrix.h"

#ifndef __BGE_SCORE_H__
#define __BGE_SCORE_H__

extern double LogBgeScore(int node, unsigned int * family, int no_parents, int alpha_mu, int alpha_omega,
              double log_prefactor,double* log_gamma_ratio_table, Bge_Matrix* prior_matrix,
              Bge_Matrix* posterior_matrix, Bge_Matrix* data);

extern double log_gamma_ratio(int samples, int alpha_omega, int total_no_nodes, int no_nodes);

extern double log_prefactor(int no_samples, int alpha_mu);

extern double * create_log_gamma_ratio_table(int samples, int alpha_omega, int no_vars);

extern double LogBgeScoreWithoutParents( int node, int alpha_mu, int alpha_omega, double log_prefactor,
              double* log_gamma_ratio_table, Bge_Matrix* prior_matrix, Bge_Matrix* posterior_matrix, Bge_Matrix* data);

extern double LogBgeScoreWithParents(int node, unsigned int *family, int no_parents,
              int alpha_mu, int alpha_omega, double log_prefactor,
              double *log_gamma_ratio_table, Bge_Matrix* prior_matrix,
              Bge_Matrix* posterior_matrix, Bge_Matrix* data);

#endif
