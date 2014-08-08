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

#include <assert.h>
#include <algos.h>
#include <cn_access.h>
#include <ignore_reg_access.h>
#include <output.h>
#include <ctype.h>
#include <dbg.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <fenv.h>

float min_m_prob = 0.8;
float min_s_prob = 0.95;
float norm_cont = 0.1;
float ref_b = 0.95;
float prior_m_prob = 0.000006;
float prior_s_prob = 0.0001;
int normal_cn = 2;
int tumour_cn = 2;
int min_t_cvg = 1;
int min_n_cvg = 1;

char bases_char[4] = {'A','C','G','T'};
char *bases_list_str[4] = {"A","C","G","T"};

float get_min_mut_prob(){
	return min_m_prob;
}

float get_min_snp_prob(){
	return min_s_prob;
}

float get_norm_contam(){
	return norm_cont;
}

float get_ref_bias(){
	return ref_b;
}

float get_prior_mut_prob(){
	return prior_m_prob;
}

float get_prior_snp_prob(){
	return prior_s_prob;
}

int get_min_tum_cvg(){
	return min_t_cvg;
}

int get_normal_cn(){
  return normal_cn;
}

void set_normal_cn(int i){
  normal_cn = i;
}

int get_tumour_cn(){
  return tumour_cn;
}

void set_tumour_cn(int i){
  tumour_cn = i;
}

int get_min_norm_cvg(){
	return min_n_cvg;
}

void set_min_mut_prob(float f){
	min_m_prob = f;
}

void set_min_snp_prob(float f){
	min_s_prob = f;
}

void set_norm_contam(float f){
	norm_cont = f;
}

void set_ref_bias(float f){
	ref_b = f;
}

void set_prior_mut_prob(float f){
	prior_m_prob = f;
}

void set_prior_snp_prob(float f){
	prior_s_prob = f;
}

void set_min_tum_cvg(int i){
	min_t_cvg = i;
}

void set_min_norm_cvg(int i){
	min_n_cvg = i;
}

int algos_mstep_read_position(alg_bean_t *alg,int ********covs, char *chr_name, int from, int to, char *ref_base, int split_size){
	//Fetch all reads for this pos? Or a struct for a single read at that position?
	List *reads;
	char *ref_b = malloc(sizeof(char) *2);
	check_mem(ref_b);
	char *cbase = malloc(sizeof(char)*2) ;
	check_mem(cbase);

	int start = from;
	int stop;
	/*time_t rawtime;
	time(&rawtime);
	struct tm * timeinfo;*/
	//estep_position_t **positions;

	for(start=from;start<=to;start+=(split_size+1)){
		stop = start+split_size;


	//for(start=from;start<=to;start++){
		//stop = start;
		if(stop>to){
			stop = to;
		}
	//estep_position_t **positions;

  		/*time (&rawtime);
  		timeinfo = localtime (&rawtime);
		//printf ("Fetching reads %d-%d\n", start,stop);*/
		reads = bam_access_get_reads_at_this_pos(chr_name,start,stop,0,alg);
		/*time (&rawtime);
  		timeinfo = localtime (&rawtime);
		printf ("Reads retrieved: %s", asctime(timeinfo));*/
		check(reads >= 0,"Error retrieving read positions for section %s:%d-%d.",chr_name,start,stop);
		LIST_FOREACH(reads, first, next, cur){
			read_pos_t *pos_t = cur->value;
			if(pos_t->ref_pos >= start && pos_t->ref_pos <= stop){
				//char ref_b_up = toupper(ref_base[((pos_t->ref_pos)-from)]);
				char ref_b_up = toupper(ref_base[((pos_t->ref_pos)-from)]);
				if((ref_b_up != 'N') && ((ref_b_up == 'A') || (ref_b_up == 'C') || (ref_b_up == 'G') || (ref_b_up = 'T')) ){
					//int lane_i = alg_bean_get_index_for_str_arr(alg->lane,pos_t->lane);
					//check(lane_i>=0,"Error calculating lane index.");
					int rpos_i = alg_bean_get_index_for_read_pos_prop_arr(alg->rd_pos,pos_t->rd_pos,pos_t->rd_len);
					check(rpos_i>=0,"Error calculating read position index.");
					int mq_i = alg_bean_get_index_for_intrange_arr(alg->map_qual,pos_t->map_qual);
					check(mq_i>=0,"Error calculating map qual index.");
					int bq_i = alg_bean_get_index_for_intrange_arr(alg->base_qual,pos_t->base_qual);
					check(bq_i>=0,"Error calculating base qual index.");
					sprintf(ref_b,"%c",ref_b_up);
					int rbase_i = alg_bean_get_index_for_char_arr(alg->ref_base,ref_b);
					check(rbase_i>=0,"Error calculating ref base index for base %s.",ref_b);
					sprintf(cbase,"%c",pos_t->called_base);
					int callbase_i = alg_bean_get_index_for_char_arr(alg->call_base,cbase);
					check(callbase_i>=0,"Error calculating called base index.");
					covs[pos_t->read_order]
								[pos_t->strand]
								[pos_t->lane_i]
								[rpos_i]
								[mq_i]
								[bq_i]
								[rbase_i]
								[callbase_i]++;
				}
			}
			//free(pos_t->lane);
		}//End of LIST_FOREACH
		List_clear_destroy(reads);
	}//End of each section.

	free(cbase);
	free(ref_b);
	return 0;
error:
	List_clear_destroy(reads);
	if(cbase) free(cbase);
	if(ref_b) free(ref_b);
	return -1;
}

void destroy_positions_array(estep_position_t **pos_arr, int size){
	int i=0;
	if(pos_arr != NULL){
		for(i=0;i<size;i++){
			if(pos_arr[i] != NULL){
				if(pos_arr[i]->genos != NULL){
					genotype_destroy_genotype_store(pos_arr[i]->genos);
				}
				if(pos_arr[i]->ref_base) free(pos_arr[i]->ref_base);
				if(pos_arr[i]->norm_fwd_cvg != NULL) free(pos_arr[i]->norm_fwd_cvg);
				if(pos_arr[i]->norm_rev_cvg != NULL) free(pos_arr[i]->norm_rev_cvg);
				if(pos_arr[i]->tum_fwd_cvg != NULL) free(pos_arr[i]->tum_fwd_cvg);
				if(pos_arr[i]->tum_rev_cvg != NULL) free(pos_arr[i]->tum_rev_cvg);
				free(pos_arr[i]);
			}
		}
		free(pos_arr);
	}
	return;
}

void destroy_position(estep_position_t *pos_arr){
	if(pos_arr->genos != NULL){
		genotype_destroy_genotype_store(pos_arr->genos);
	}
	if(pos_arr->ref_base) free(pos_arr->ref_base);
	if(pos_arr->norm_fwd_cvg != NULL) free(pos_arr->norm_fwd_cvg);
	if(pos_arr->norm_rev_cvg != NULL) free(pos_arr->norm_rev_cvg);
	if(pos_arr->tum_fwd_cvg != NULL) free(pos_arr->tum_fwd_cvg);
	if(pos_arr->tum_rev_cvg != NULL) free(pos_arr->tum_rev_cvg);
	free(pos_arr);
	return;
}

inline long double algos_calculate_per_base_normal_contamination(int norm_copy_no,int tum_copy_no){
	long double per_base_norm = ((long double)norm_cont * (long double)norm_copy_no) /
								(
									(
										(
											(long double)1 - (long double)norm_cont
										)
										*
										(long double)tum_copy_no
									)
									+
									(
										(long double)norm_cont * (long double)norm_copy_no
									)

								);
	return per_base_norm;
}

int algos_run_per_read_estep_maths(genotype_store_t *genos,read_pos_t *read, int ref_base_idx, long double base_norm_contam){
	//Normal read
	if(read->normal==1){

		//Calculate initial normal probabilities.
		long double res = genos->ref_geno_norm_prob + read->ref_base_probs[ref_base_idx];
		genos->ref_geno_norm_prob = res;

		//iterate through from 0 to highest available of genos->hom_count, genos->het_norm_count,
		int iter=0;
		for(iter=0;iter<genos->norm_max;iter++){
			//Hom snps

			if(iter<genos->hom_count){
				long double ans = genos->hom_snp_genotypes[iter]->prob + read->ref_base_probs[(genos->hom_snp_genotypes[iter]->norm_geno->var_base_idx)];
				genos->hom_snp_genotypes[iter]->prob = ans;
			}//End of iteration through hom snps

			//Het snps
			if(iter<genos->het_norm_count){
				long double norm_var_base_prop = genos->het_snp_norm_genotypes[iter]->norm_geno->var_base_prop;
				long double nu_fx = (long double)ref_b * norm_var_base_prop;
				long double tmp_psi_var_prob_norm = nu_fx / (nu_fx + ((long double)1 - norm_var_base_prop));
				long double ans = genos->het_snp_norm_genotypes[iter]->prob +
												(
													logl(
														expl( read->ref_base_probs[ref_base_idx] + logl((long double)1 - tmp_psi_var_prob_norm) )
															+
														expl( read->ref_base_probs[(genos->het_snp_norm_genotypes[iter]->norm_geno->var_base_idx)] + logl(tmp_psi_var_prob_norm) )
													)
												);
				genos->het_snp_norm_genotypes[iter]->prob = ans;
			}//End of iteration through het snps
		}

	}else if(read->normal==0){//A tumour read
		long double res = genos->ref_geno_tum_prob + read->ref_base_probs[ref_base_idx];
		genos->ref_geno_tum_prob = res;

		//iterate through from 0 to highest available of genos->somatic_count, genos->hom_count, genos->het_count
		int iter=0;
		for(iter=0;iter<genos->tum_max;iter++){
			//Somatics

			if(iter<genos->somatic_count){
				long double nu_xf = (long double)ref_b *  genos->somatic_genotypes[iter]->tum_geno->var_base_prop;
				long double tmp_psi_var_prob = nu_xf * ((long double)1 - base_norm_contam) / (nu_xf + ((long double)1 - genos->somatic_genotypes[iter]->tum_geno->var_base_prop));
				//Calculate read component probability.
				long double log_tmp_psi_var_prob = logl(tmp_psi_var_prob);
				long double log_1_minus_tmp_psi_var_prob = logl((long double)1 - tmp_psi_var_prob);

				long double ans = genos->somatic_genotypes[iter]->prob +
																	logl(
																		(
																			expl( (read->ref_base_probs[ref_base_idx] + log_1_minus_tmp_psi_var_prob) )
																				+
																			expl( (read->ref_base_probs[genos->somatic_genotypes[iter]->tum_geno->var_base_idx] + log_tmp_psi_var_prob) )
																		)
																	);
				genos->somatic_genotypes[iter]->prob = ans;

			}//End of somatics

			//Hom snps
			if(iter<genos->hom_count){
				long double ans = genos->hom_snp_genotypes[iter]->prob + read->ref_base_probs[(genos->hom_snp_genotypes[iter]->tum_geno->var_base_idx)];
				genos->hom_snp_genotypes[iter]->prob = ans;
			}//End of hom snps

			//Het snps
			if(iter<genos->het_count){
				long double nu_fx = (long double)ref_b * genos->het_snp_genotypes[iter]->norm_geno->var_base_prop;
				long double tmp_psi_var_prob_norm = nu_fx / (nu_fx + ((long double)1 - genos->het_snp_genotypes[iter]->norm_geno->var_base_prop));

				long double nu_fx_two = (long double)ref_b * genos->het_snp_genotypes[iter]->tum_geno->var_base_prop;
				long double tmp_psi_var_prob_tum = ((nu_fx_two * ((long double)1 - (long double)base_norm_contam))
															/
														(nu_fx_two + ((long double)1 - genos->het_snp_genotypes[iter]->tum_geno->var_base_prop)))
														+
														( (long double)base_norm_contam * tmp_psi_var_prob_norm);
				long double tum_res = genos->het_snp_genotypes[iter]->prob +
												(
													logl(
														expl( read->ref_base_probs[ref_base_idx] + logl((long double)1 - tmp_psi_var_prob_tum) )
															+
														expl( read->ref_base_probs[(genos->het_snp_genotypes[iter]->norm_geno->var_base_idx)] + logl(tmp_psi_var_prob_tum) )
													)
												);

				genos->het_snp_genotypes[iter]->prob = tum_res;
			}//End of het snps

		}//End of iteration.



	}//End of it this is normal or tumour read.
	return 0;
}

inline long double calculateLogSumExpNormFactor(estep_position_t *pos,long double norm_factor_max){
	long double sum_tot = 0.0;
   //Reference
   sum_tot += expl(pos->genos->ref_genotype->prob - norm_factor_max);
  // printf("%s/%s = %Le\n",genotype_get_genotype_t_as_string(pos->genos->ref_genotype->norm_geno),
   	//			genotype_get_genotype_t_as_string(pos->genos->ref_genotype->tum_geno),pos->genos->ref_genotype->prob);
  //
   //Somatic
   int somatic_i = 0;
   for(somatic_i = 0;somatic_i<pos->genos->somatic_count;somatic_i++){
   	sum_tot += expl(pos->genos->somatic_genotypes[somatic_i]->prob - norm_factor_max);
   //	printf("%s/%s = %Le\n",genotype_get_genotype_t_as_string(pos->genos->somatic_genotypes[somatic_i]->norm_geno),
   //				genotype_get_genotype_t_as_string(pos->genos->somatic_genotypes[somatic_i]->tum_geno),pos->genos->somatic_genotypes[somatic_i]->prob);
   }

   //het
   int het_snp_it=0;
   for(het_snp_it=0;het_snp_it<pos->genos->het_count;het_snp_it++){
   	sum_tot += expl(pos->genos->het_snp_genotypes[het_snp_it]->prob - norm_factor_max);
   //	printf("%s/%s = %Le\n",genotype_get_genotype_t_as_string(pos->genos->het_snp_genotypes[het_snp_it]->norm_geno),
   //				genotype_get_genotype_t_as_string(pos->genos->het_snp_genotypes[het_snp_it]->tum_geno),pos->genos->het_snp_genotypes[het_snp_it]->prob);
   }

   //Hom
   int hom_snp_i=0;
	for(hom_snp_i=0;hom_snp_i<pos->genos->hom_count;hom_snp_i++){
		sum_tot += expl(pos->genos->hom_snp_genotypes[hom_snp_i]->prob - norm_factor_max);
	//	printf("%s/%s = %Le\n",genotype_get_genotype_t_as_string(pos->genos->hom_snp_genotypes[hom_snp_i]->norm_geno),
	//				genotype_get_genotype_t_as_string(pos->genos->hom_snp_genotypes[hom_snp_i]->tum_geno),pos->genos->hom_snp_genotypes[hom_snp_i]->prob);
	}

  	//Return final logged value.
  	long double norm_factor = norm_factor_max + logl(sum_tot);
  	return norm_factor;
}

inline void finalise_probabilities_and_find_top_prob(estep_position_t *pos,long double norm_factor_max){

	combined_genotype_t *top_geno = NULL;
	combined_genotype_t *sec_geno = NULL;

	//We find the normalisation factor via the log-sum-exp method to avoid underflow, then we can calculate the final probabilities
	long double norm_factor = calculateLogSumExpNormFactor(pos,norm_factor_max);
	//printf("NORM FACTOR: %5.1Le\n",norm_factor);
	//Iterate through every genotype and normalise it. Adding it to a list.
	int iter = 0;
	long double somatic_sum = 0.0;
	long double snp_sum = 0.0;
	for(iter=0;iter<pos->genos->total_max;iter++){
		//somatics
		if(iter<pos->genos->somatic_count){
			pos->genos->somatic_genotypes[iter]->prob = expl(pos->genos->somatic_genotypes[iter]->prob - norm_factor);

			somatic_sum += pos->genos->somatic_genotypes[iter]->prob;
			if(top_geno==NULL){
				top_geno = pos->genos->somatic_genotypes[iter];
			}else if(pos->genos->somatic_genotypes[iter]->prob > top_geno->prob){
				sec_geno = top_geno;
				top_geno = pos->genos->somatic_genotypes[iter];
			}else if(sec_geno==NULL || pos->genos->somatic_genotypes[iter]->prob > sec_geno->prob){
				sec_geno = pos->genos->somatic_genotypes[iter];
			}
		}

		//het snps
		if(iter<pos->genos->het_count){
			pos->genos->het_snp_genotypes[iter]->prob = expl(pos->genos->het_snp_genotypes[iter]->prob - norm_factor);

			snp_sum += pos->genos->het_snp_genotypes[iter]->prob;
			if(pos->genos->het_snp_genotypes[iter]->prob > top_geno->prob){
				sec_geno = top_geno;
				top_geno = pos->genos->het_snp_genotypes[iter];
			}else if(sec_geno==NULL || pos->genos->het_snp_genotypes[iter]->prob > sec_geno->prob){
				sec_geno = pos->genos->het_snp_genotypes[iter];
			}
		}

		//hom snps
		if(iter<pos->genos->hom_count){
			pos->genos->hom_snp_genotypes[iter]->prob = expl(pos->genos->hom_snp_genotypes[iter]->prob - norm_factor);

			snp_sum += pos->genos->hom_snp_genotypes[iter]->prob;
			if(pos->genos->hom_snp_genotypes[iter]->prob > top_geno->prob){
				sec_geno = top_geno;
				top_geno = pos->genos->hom_snp_genotypes[iter];
			}else if(sec_geno==NULL || pos->genos->hom_snp_genotypes[iter]->prob > sec_geno->prob){
				sec_geno = pos->genos->hom_snp_genotypes[iter];
			}
		}
	}


	//normalise the reference genotype probability.
	pos->genos->ref_genotype->prob = expl(pos->genos->ref_genotype->prob - norm_factor);

	if(pos->genos->ref_genotype->prob > top_geno->prob){
		sec_geno = top_geno;
		top_geno = pos->genos->ref_genotype;
	}else if(pos->genos->ref_genotype->prob > sec_geno->prob){
		sec_geno = pos->genos->ref_genotype;
	}

	//Sort the list then set the two top genotypes.

	pos->top_geno = top_geno;
	pos->sec_geno = sec_geno;

	//Round the total probabilities to make sure we're not losing ones that round to the cutoff value.s
	pos->total_snp_prob = snp_sum;
	pos->total_mut_prob = somatic_sum;
	return;
}

void algos_run_per_position_estep_maths(estep_position_t *pos){
	//Put together the reference genotype probability.
	long double norm_factor_max = -DBL_MAX;
	pos->genos->ref_genotype->prob = 0.0;
	long double ref_prob = logl(((long double)1 - (long double)prior_s_prob) - (long double)prior_m_prob)
							+ pos->genos->ref_geno_norm_prob + pos->genos->ref_geno_tum_prob;

	pos->genos->ref_genotype->prob = ref_prob;
	if(pos->genos->ref_genotype->prob > norm_factor_max){
		norm_factor_max = pos->genos->ref_genotype->prob;
	}

	//Iterate through the 3 arrays of genotypes at once
	int iter=0;
	for(iter=0;  iter<pos->genos->total_max; iter++){
		//Somatics
		if(iter<pos->genos->somatic_count){
			long double res = logl((long double)prior_m_prob) - logl((long double)pos->genos->somatic_count)
							+ pos->genos->ref_geno_norm_prob + pos->genos->somatic_genotypes[iter]->prob;

			pos->genos->somatic_genotypes[iter]->prob = res;

			if(res > norm_factor_max){
				norm_factor_max = res;
			}
		}

		//Het snps
		if(iter<pos->genos->het_count){
			int het_norm_i=0;
			for(het_norm_i=0;het_norm_i<pos->genos->het_norm_count;het_norm_i++){
				//find the correct het normal genotype
				if(genotype_equals(pos->genos->het_snp_norm_genotypes[het_norm_i]->norm_geno,pos->genos->het_snp_genotypes[iter]->norm_geno)==1){
					long double res = logl((long double)prior_s_prob) - logl((long double)(pos->genos->het_count + pos->genos->hom_count))
								+ pos->genos->het_snp_norm_genotypes[het_norm_i]->prob + pos->genos->het_snp_genotypes[iter]->prob;

					pos->genos->het_snp_genotypes[iter]->prob = res;

					if(pos->genos->het_snp_genotypes[iter]->prob > norm_factor_max){
						norm_factor_max = pos->genos->het_snp_genotypes[iter]->prob;
					}
				}
			}
		}


		//Hom snps
		if(iter<pos->genos->hom_count){
			//pos->genos->hom_snp_genotypes[iter]->prob = pos->genos->hom_snp_genotypes[iter]->prob;
		 	if(pos->genos->hom_snp_genotypes[iter]->prob > norm_factor_max){
				norm_factor_max = pos->genos->hom_snp_genotypes[iter]->prob;
			}
		}
	}

	finalise_probabilities_and_find_top_prob(pos, norm_factor_max);

}

inline int algos_get_read_specific_indices(alg_bean_t *alg, read_pos_t *pos_t, int *rpos_i, int *mq_i, int *bq_i, int *callbase_i){

	*rpos_i = alg_bean_get_index_for_read_pos_prop_arr(alg->rd_pos,pos_t->rd_pos,pos_t->rd_len);
	check(*rpos_i>=0,"Error calculating read position index.");
	*mq_i = alg_bean_get_index_for_intrange_arr(alg->map_qual,pos_t->map_qual);
	check(*mq_i>=0,"Error calculating map qual index.");
	*bq_i = alg_bean_get_index_for_intrange_arr(alg->base_qual,pos_t->base_qual);
	check(*bq_i>=0,"Error calculating base qual index.");
	char *cbase = malloc(sizeof(char)*2) ;
	//Build a list of probabilities for each possible ref base at this position.
	check_mem(cbase);
	sprintf(cbase,"%c",pos_t->called_base);
	*callbase_i = alg_bean_get_index_for_char_arr(alg->call_base,cbase);
	free(cbase);
	check(*callbase_i>=0,"Error calculating called base index.");

	return 0;
error:
	return -1;
}

int algos_estep_read_position(alg_bean_t *alg,long double ********prob_arr, char *chr_name, int from, int to,
						char *ref_base, char *norm_cn, char *tum_cn, FILE *snp_out, FILE *tum_out, FILE *debug_output, int split_size){

	List *reads;
	//If this region is larger than n bases, split it into regions of length n.
	int start = from;
	int stop;
	estep_position_t *pos;

	//estep_position_t **positions;
	for(start=from;start<=to;start+=(split_size+1)){
		stop = start+split_size;
		if(stop>to){
			stop = to;
		}
		int arr_size = (stop-start)+1;
		int *pos_cvg = malloc(sizeof(int) * arr_size);
		check_mem(pos_cvg);
		int i=0;
		for(i=0; i<arr_size; i++){
			pos_cvg[i]=0;
		}
		pos = NULL;
		//Fetch all reads for this region - returned sorted by position.
		reads = bam_access_get_reads_at_this_pos(chr_name,start,stop,1,alg);
		check(reads >= 0,"Error retrieving read positions for section %s:%d-%d.",chr_name,start,stop);
		//Iterate through each read once to generate position beans.
		//Change position when we reach a new position as list is sorted
		if(List_count(reads) > 0){
			LIST_FOREACH(reads, first, next, cur){
				read_pos_t *pos_t = cur->value;
				if(pos == NULL || pos->ref_pos < pos_t->ref_pos){
					if(pos != NULL){
						//Finish the last position we ran over
						if(pos->ref_base_idx >= 0 && pos->norm_cn > 0 && pos->tum_cn > 0 && pos->total_cvg_norm > 0 && pos->total_cvg_tum > 0){
							algos_run_per_position_estep_maths(pos);
							//printf("Total_mut");
							int rounded_mut = (int)(round(floorl(pos->total_mut_prob*100 + 0.5)));
							int rounded_snp = (int)(round(floorl(pos->total_snp_prob*100 + 0.5)));
							if(rounded_mut >= ((int)(min_m_prob * 100))){
							//if((((floorl(pos->total_mut_prob*100 + 0.5)/100) - (long double)min_m_prob)) >= (long double)0){
								int write = output_vcf_variant_position(pos,tum_out,chr_name);
								check(write==0,"Error writing mutation to file.");
							}else if(rounded_snp >= (int)(min_s_prob * 100)){
								int write = output_vcf_variant_position(pos,snp_out,chr_name);
								check(write==0,"Error writing germline to file.");
							}
							if(debug_output != NULL){
								int write = output_vcf_variant_position(pos,debug_output,chr_name);
								check(write==0,"Error writing debug output.");
							}
							pos_cvg[pos->ref_pos-start] = 1;
						}
						destroy_position(pos);
					}
					pos = malloc(sizeof(estep_position_t));
					check_mem(pos);
					pos->ref_pos = pos_t->ref_pos;
					pos->norm_cn = cn_access_get_copy_number_for_location(norm_cn,chr_name,pos_t->ref_pos,1);
					check(pos->norm_cn != -1,"Error fetching normal copy number for %s:%d\n",chr_name,pos_t->ref_pos);
					if(pos->norm_cn < 2){
            pos->norm_cn = normal_cn;
					}
					pos->tum_cn = cn_access_get_copy_number_for_location(tum_cn,chr_name,pos_t->ref_pos,0);
					check(pos->tum_cn!=-1,"Error fetching tumour copy number for %s:%d\n",chr_name,pos_t->ref_pos);
					if(pos->tum_cn < 2){
					  pos->tum_cn = tumour_cn;
					}
					char *ref_b = malloc(sizeof(char) *2);
					check_mem(ref_b);
					char ref_b_up = toupper(ref_base[((pos_t->ref_pos)-from)]);
					int pr_chk = sprintf(ref_b,"%c",ref_b_up);
					check(pr_chk==1,"Error writing uppercase ref base to string.");
					//Uppercase ref base to ensure all comparisons are like for like
					pos->ref_base = ref_b;
					if((strcmp(ref_b,"A") == 0 || strcmp(ref_b,"C") == 0 || strcmp(ref_b,"G") == 0 || strcmp(ref_b,"T") == 0)
																				&& pos->norm_cn > 0 && pos->tum_cn > 0){
						pos->ref_base_idx = alg_bean_get_index_for_char_arr(alg->ref_base,ref_b);
						//build up genotype stores for genotype_store_t *genos
						pos->genos = genotype_generate_genotype_list_for_cn_and_ref_base(pos->norm_cn, pos->tum_cn, ref_b);

						//Setup initial hom snp probabilities.
						long double init_hom_snp = logl((long double)prior_s_prob) - logl((long double) (pos->genos->het_count + pos->genos->hom_count));
						int hom_snp_i = 0;
						for(hom_snp_i=0;hom_snp_i<pos->genos->hom_count;hom_snp_i++){
							pos->genos->hom_snp_genotypes[hom_snp_i]->prob = init_hom_snp;
						}

						check(pos->genos != NULL,"Error generating genotypes.");
						//The coverage stats genotypes.
						pos->norm_fwd_cvg = genotype_init_genotype();
						pos->norm_rev_cvg = genotype_init_genotype();
						pos->tum_fwd_cvg = genotype_init_genotype();
						pos->tum_rev_cvg = genotype_init_genotype();
						pos->total_cvg_norm = 0;
						pos->total_cvg_tum = 0;

						//Per base normal contamination value...
						pos->base_norm_contam = algos_calculate_per_base_normal_contamination(pos->norm_cn,pos->tum_cn);
					}else{
						pos->ref_base_idx = -1;
						pos->norm_fwd_cvg = NULL;
						pos->norm_rev_cvg = NULL;
						pos->tum_fwd_cvg = NULL;
						pos->tum_rev_cvg = NULL;
						pos->genos = NULL;

					}
				}
				//read_pos_t struct has copy number and ref base is a usable base.
				if(pos->ref_base_idx >= 0 && pos->norm_cn > 0 && pos->tum_cn > 0){
					//Adding to the counts at this position.
					if(pos_t->normal == 1){
						if(pos_t->strand == 1){
							genotype_add_base_to_count(pos->norm_rev_cvg,pos_t->called_base);
						}else{
							genotype_add_base_to_count(pos->norm_fwd_cvg,pos_t->called_base);
						}
						pos->total_cvg_norm++;
					}else if(pos_t->normal == 0){
						if(pos_t->strand == 1){
							genotype_add_base_to_count(pos->tum_rev_cvg,pos_t->called_base);
						}else{
							genotype_add_base_to_count(pos->tum_fwd_cvg,pos_t->called_base);
						}
						pos->total_cvg_tum++;
					}
					//Get the index for things we'll use during analysis
					int rpos_i,mq_i,bq_i,callbase_i;
					int ok = algos_get_read_specific_indices(alg,pos_t,&rpos_i,&mq_i,&bq_i,&callbase_i);
					check(ok==0,"Error fetching index information for read specific info.");

					//Assign the probability for each reference base, so we only look it up once per read.
					int ref_base_iter=0;
					for(ref_base_iter=0;ref_base_iter<4;ref_base_iter++){
						pos_t->ref_base_probs[ref_base_iter] = 0.0;
						int tmp_index = alg_bean_get_index_for_char_arr(alg->ref_base,bases_list_str[ref_base_iter]);

						pos_t->ref_base_probs[ref_base_iter] = prob_arr[pos_t->read_order]
								[pos_t->strand]
								[pos_t->lane_i]
								[rpos_i]
								[mq_i]
								[bq_i]
								[tmp_index]
								[callbase_i];
					}

					//With all the indexes we need for now, and the genotypes calculated we can run the mathematical part...
					algos_run_per_read_estep_maths(pos->genos,pos_t,pos->ref_base_idx,pos->base_norm_contam);
				}
				//free(pos_t->lane);
			}//end of iteration through reads.
			if(pos->ref_base_idx >= 0 && pos->norm_cn > 0 && pos->tum_cn > 0 && pos->total_cvg_norm > 0 && pos->total_cvg_tum > 0){
				algos_run_per_position_estep_maths(pos);
				int rounded_mut = (int)(round(floorl(pos->total_mut_prob*100 + 0.5)));
				int rounded_snp = (int)(round(floorl(pos->total_snp_prob*100 + 0.5)));
				if(rounded_mut >= ((int)(min_m_prob * 100))){
				//if((((floorl(pos->total_mut_prob*100 + 0.5)/100) - (long double)min_m_prob)) >= (long double)0){
					int write = output_vcf_variant_position(pos,tum_out,chr_name);
					check(write==0,"Error writing mutation to file.");
				}else if(rounded_snp >= (int)(min_s_prob * 100)){
					int write = output_vcf_variant_position(pos,snp_out,chr_name);
					check(write==0,"Error writing germline to file.");
				}
				if(debug_output != NULL){
					int write = output_vcf_variant_position(pos,debug_output,chr_name);
					check(write==0,"Error writing debug to file.");
				}
				pos_cvg[pos->ref_pos-start] = 1;
			}
			destroy_position(pos);
		}//End of if there are reads.
		List_clear_destroy(reads);
		i=0;
		for(i=0; i<arr_size; i++){
			if(pos_cvg[i]==0){
				int sect_start = start+i;
				int sect_stop = start+i;
				while(i<arr_size && pos_cvg[i]==0 ){
					sect_stop = start+i;
					i++;
				}
				output_append_position_to_no_analysis(chr_name, sect_start, sect_stop);
			}
		}
		free(pos_cvg);
	}

	//Code to destroy positions array
	clear_copy_number_store();
	genotype_clear_genotype_cache();
	return 0;
error:
	if(pos) destroy_position(pos);
	genotype_clear_genotype_cache();
	clear_copy_number_store();
	List_clear_destroy(reads);
	return -1;
}
