/**   LICENSE
* Copyright (c) 2014-2019 Genome Research Ltd.
*
* Author: Cancer Genome Project cgpit@sanger.ac.uk
*
* This file is part of CaVEMan.
*
* CaVEMan is free software: you can redistribute it and/or modify it under
* the terms of the GNU Affero General Public License as published by the Free
* Software Foundation; either version 3 of the License, or (at your option) any
* later version.
*
* This program is distributed in the hope that it will be useful, but WITHOUT
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
* details.
*
* You should have received a copy of the GNU Affero General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*
*    1. The usage of a range of years within a copyright statement contained within
*    this distribution should be interpreted as being equivalent to a list of years
*    including the first and last year specified and all consecutive years between
*    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
*    2009, 2011-2012’ should be interpreted as being identical to a statement that
*    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
*    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
*    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
*    2009, 2010, 2011, 2012’."
*
*/

#ifndef _algos_h
#define _algos_h

#include <alg_bean.h>
#include <genotype.h>
#include <bam_access.h>
#include <stdint.h>
#include "zlib.h"

typedef struct estep_position_t{
	uint8_t norm_cn;
	uint8_t tum_cn;
	int8_t ref_base_idx;
	uint32_t ref_pos;
	uint32_t total_cvg_tum;
	uint32_t total_cvg_norm;
	long double base_norm_contam;
	long double total_snp_prob;
	long double total_mut_prob;
    long double total_mut_allele_prob;
	genotype_t *norm_fwd_cvg;
	genotype_t *norm_rev_cvg;
	genotype_t *tum_fwd_cvg;
	genotype_t *tum_rev_cvg;
	combined_genotype_t *top_geno;
	combined_genotype_t *sec_geno;
	genotype_store_t *genos;
	char *ref_base;
} estep_position_t;

int algos_mstep_read_position(alg_bean_t *alg,uint64_t ********covs, char *chr_name, uint32_t from, uint32_t to, char *ref_base, int split_size);
int algos_estep_read_position(alg_bean_t *alg,long double ********prob_arr, char *chr_name, uint32_t from, uint32_t to, char *ref_base,
												char *norm_cn, char *tum_cn, gzFile snp_out, gzFile tum_out, gzFile dbg, int split_size);

int algos_check_var_position_alleles(estep_position_t *pos, char *chr_name, char *type, uint8_t warnings);
long double algos_calculate_per_base_normal_contamination(uint8_t norm_copy_no,uint8_t tum_copy_no);
void finalise_probabilities_and_find_top_prob(estep_position_t *pos,long double norm_factor);
int algos_run_per_read_estep_maths(genotype_store_t *genos,read_pos_t *read, int8_t ref_base_idx, long double base_norm_contam);
void algos_run_per_position_estep_maths(estep_position_t *pos);

void set_snp_warnings();
void set_min_mut_prob(float f);
void set_min_snp_prob(float f);
void set_norm_contam(float f);
void set_ref_bias(float f);
void set_prior_mut_prob(float f);
void set_prior_snp_prob(float f);
void set_max_tum_cvg(int i);
void set_min_tum_cvg(int i);
void set_min_norm_cvg(int i);
int get_normal_cn();
void set_normal_cn(int i);
int get_tumour_cn();
void set_tumour_cn(int i);

float get_min_mut_prob();
float get_min_snp_prob();
float get_norm_contam();
float get_ref_bias();
float get_prior_mut_prob();
float get_prior_snp_prob();
int get_min_tum_cvg();
int get_min_norm_cvg();

#endif
