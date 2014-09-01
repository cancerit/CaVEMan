/**   LICENSE
* Copyright (c) 2014 Genome Research Ltd.
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
*/

#ifndef _algos_h
#define _algos_h

#include <alg_bean.h>
#include <genotype.h>
#include <bam_access.h>

typedef struct estep_position_t{
	int ref_pos;
	genotype_store_t *genos;
	int norm_cn;
	int tum_cn;
	int ref_base_idx;
	char *ref_base;
	long double base_norm_contam;
	genotype_t *norm_fwd_cvg;
	genotype_t *norm_rev_cvg;
	genotype_t *tum_fwd_cvg;
	genotype_t *tum_rev_cvg;
	int total_cvg_tum;
	int total_cvg_norm;
	long double total_snp_prob;
	long double total_mut_prob;
	combined_genotype_t *top_geno;
	combined_genotype_t *sec_geno;
} estep_position_t;

int algos_mstep_read_position(alg_bean_t *alg,int ********covs, char *chr_name, int from, int to, char *ref_base, int split_size);
int algos_estep_read_position(alg_bean_t *alg,long double ********prob_arr, char *chr_name, int from, int to, char *ref_base,
												char *norm_cn, char *tum_cn, FILE *snp_out, FILE *tum_out, FILE *dbg, int split_size);

int algos_check_var_position_alleles(estep_position_t *pos, char *chr_name, char *type);
inline long double algos_calculate_per_base_normal_contamination(int norm_copy_no,int tum_copy_no);
inline void finalise_probabilities_and_find_top_prob(estep_position_t *pos,long double norm_factor);
inline int algos_run_per_read_estep_maths(genotype_store_t *genos,read_pos_t *read, int ref_base_idx, long double base_norm_contam);
inline void algos_run_per_position_estep_maths(estep_position_t *pos);

void set_min_mut_prob(float f);
void set_min_snp_prob(float f);
void set_norm_contam(float f);
void set_ref_bias(float f);
void set_prior_mut_prob(float f);
void set_prior_snp_prob(float f);
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
