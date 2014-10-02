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

#include "minunit.h"
#include <algos.h>
#include <covs_access.h>
#include <bam_access.h>
#include <alg_bean.h>
#include <output.h>
#include <genotype.h>
#include <dbg.h>
#include <math.h>

char *norm = "testData/wt.bam";
char *tum = "testData/mt.bam";
char *mut_norm = "testData/testing_wt.bam";
char *mut_tum = "testData/testing_mt.bam";
char *mstep_mut = "testData/mstep_test_mt.bam";
char *mstep_norm = "testData/mstep_test_wt.bam";
char *mut_probs = "testData/test_mut_probs_array";
char *mut_alg = "testData/test_mut_alg";
char *mut_mt_cn = "testData/mc.cave.cn";
char *mut_wt_cn = "testData/wc.cave.cn";
char *test_snp_out = "testData/snp.vcf";
char *test_mut_out = "testData/mut.vcf";
char *test_dbg_out = "testData/dbg.vcf";
char *test_no_anal_out = "testData/no_analysis.bed";

char *test_algos_mstep_read_position(){
	//algos_mstep_read_position(alg_bean_t *alg,int ********covs, char *chr_name, int from, int to, char *ref_base);
	alg_bean_t *alg = alg_bean_generate_default_alg_bean(norm,tum);
	mu_assert(bam_access_openbams(norm, tum)==0,"Bams not opened.\n");
	int ********arr = covs_access_generate_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	mu_assert(arr != NULL,"Array not properly created.\n");
	char *chr = "22";
	int from = 17619559;
	int to = 17619559;
	char *ref_base = "A";
	algos_mstep_read_position(alg, arr, chr, from, to, ref_base, 50000);
	//check we have a count of 2 in the expected place...
	mu_assert(arr[0][1][1][4][1][4][0][3] == 1, "Incorrect data returned.\n");
	covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),arr);
	alg_bean_destroy(alg);
	return NULL;
}

char *test_algos_mstep_read_position_two(){
	//2	38000243	.	G	A	.	.	DP=178;MP=1.7e-84;GP=1.0e+00;TG=AG/AAGG;TP=1.0e+00;SG=AG/AGGG;SP=1.9e-03	GT:AF:CF:GF:TF:AR:CR:GR:TR:PM	0|1:17:0:28:0:19:0:17:0:4.4e-01	0|1:25:0:30:0:18:0:24:0:4.4e-01
	//char *mstep_mut = "tests/mstep_test_mt.bam";
	//char *mstep_norm = "tests/mstep_test_wt.bam";
	alg_bean_t *alg = alg_bean_generate_default_alg_bean(mstep_norm,mstep_mut);
	mu_assert(bam_access_openbams(mstep_norm,mstep_mut)==0,"Bams not opened.\n");
	int ********arr = covs_access_generate_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	mu_assert(arr != NULL,"Array not properly created.\n");
	char *chr = "2";
	int from = 38000243;
	int to = 38000243;
	char *ref_base = "G";
	algos_mstep_read_position(alg, arr, chr, from, to, ref_base, 50000);

	int sum=0;
	int i,j,k,m,n,p,r,s;
	for(i=0;i<List_count(alg->read_order);i++){
		for(j=0;j<List_count(alg->strand);j++){
			for(k=0;k<List_count(alg->lane);k++){
				for(m=0;m<List_count(alg->rd_pos);m++){
					for(n=0;n<List_count(alg->map_qual);n++){
						for(p=0;p<List_count(alg->base_qual);p++){
							for(r=0;r<List_count(alg->ref_base);r++){
								for(s=0;s<List_count(alg->call_base);s++){
									sum+=arr[i][j][k][m][n][p][r][s];
								}
							}
						}
					}
				}
			}
		}
	}
	mu_assert(sum==178,"Incorrect total count in cov array.");
	covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),arr);
	alg_bean_destroy(alg);
	return NULL;
}

int estep_no_analysis(){
	FILE *alg_file = fopen(mut_alg,"r");
	check(bam_access_openbams(mut_norm, mut_tum)==0,"Bams not opened.");
	alg_bean_t *alg = alg_bean_read_file(alg_file);
	char *ref_base = malloc(sizeof(char)*51);
	memset(ref_base,'C',51*sizeof(char));
	fclose(alg_file);
	long double ********probs = covs_access_read_probs_from_file(mut_probs,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	FILE *snp_out = fopen(test_snp_out,"w");
	check(snp_out!=NULL,"Error opening file.");
	FILE *mut_out = fopen(test_mut_out,"w");
	check(mut_out!=NULL,"Error opening file.");
	FILE *dbg_out = fopen(test_dbg_out,"w");
	check(dbg_out!=NULL,"Error opening file.");
	FILE *no_anal_out = fopen(test_no_anal_out,"w");
	check(no_anal_out!=NULL,"Error opening file.");
	output_set_no_analysis_file(no_anal_out);
	output_set_no_analysis_section_list(List_create());
	int estep = algos_estep_read_position(alg, probs,"1", 192462250, 192462300, ref_base, mut_wt_cn, mut_mt_cn, snp_out, mut_out, dbg_out, 50000);
	output_flush_no_analysis("1");
	check(estep==0,"Error running estep.");
	fclose(snp_out);
	fclose(mut_out);
	fclose(dbg_out);
	fclose(no_anal_out);

	bam_access_closebams();
	covs_access_free_prob_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),probs);
				alg_bean_destroy(alg);


	no_anal_out = fopen(test_no_anal_out,"r");
	char line[5000];
	int count = 0;
	while ( fgets(line,sizeof(line),no_anal_out) != NULL ){
		check(strncmp("1\t192462249\t192462263\n",line,(sizeof(char) * 22))==0,"Incorrect no analysis output in file.");
		count++;
	}
	check(count==1,"Wrong number of lines found in no analysis file.");
	fclose(no_anal_out);
	free(ref_base);
	unlink(test_snp_out);
	unlink(test_mut_out);
	unlink(test_dbg_out);
	unlink(test_no_anal_out);
	return 0;
error:
	return -1;
}

char *test_algos_estep_read_position_real_data_no_analysis(){
	mu_assert(estep_no_analysis()==0,"Error testing no analysis estep.");
	return NULL;
}

//1	192462357				.	C	A	.		.	DP=77;MP=1.0e+00;GP=1.3e-06;TG=CC/AAA;TP=1.0e+00;SG=CC/AAC;SP=8.3e-05	GT:AF:CF:GF:TF:AR:CR:GR:TR:PM	0|0:0:21:0:0:0:4:0:0:0.0e+00	1|1:14:7:0:0:23:8:0:0:7.1e-01
//1	192462357	958902214	C	A	.	PASS	DP=77;MP=1.0e+00;GP=2.2e-06;TG=CC/AAC;TP=5.9e-01;SG=CC/AAA;SP=4.1e-01;VT=Sub	GT:AA:CA:GA:TA:PM	0/0:0:25:0:0:0.0e+00	0/1:37:15:0:0:7.1e-01

char *test_algos_estep_read_position(){
	FILE *alg_file = fopen(mut_alg,"r");
	mu_assert(bam_access_openbams(mut_norm, mut_tum)==0,"Bams not opened.\n");
	alg_bean_t *alg = alg_bean_read_file(alg_file);
	fclose(alg_file);
	long double ********probs = covs_access_read_probs_from_file(mut_probs,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	FILE *snp_out = fopen(test_snp_out,"w");
	FILE *mut_out = fopen(test_mut_out,"w");
	FILE *dbg_out = fopen(test_dbg_out,"w");
	FILE *no_anal_out = fopen(test_no_anal_out,"w");
	output_set_no_analysis_file(no_anal_out);
	output_set_no_analysis_section_list(List_create());
	int estep = algos_estep_read_position(alg, probs,"1", 192462357, 192462357, "C", mut_wt_cn, mut_mt_cn, snp_out, mut_out, dbg_out, 50000);
	mu_assert(estep==0,"Error running estep.");
	fclose(snp_out);
	fclose(mut_out);
	fclose(dbg_out);
	fclose(no_anal_out);

	bam_access_closebams();
	covs_access_free_prob_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),probs);
				alg_bean_destroy(alg);

	mut_out = fopen(test_mut_out,"r");
	char line[5000];
	int count = 0;
	while ( fgets(line,sizeof(line),mut_out) != NULL ){
		mu_assert(strncmp("1\t192462357",line,(sizeof(char) * 11))==0,"Incorrect mutation output in file.");
		//printf("%s",line);
		//mu_assert(strcmp(
			//"1\t192462357\t.\tC\tA\t.\t.\tDP=77;MP=1.0e+00;GP=2.2e-06;TG=CC/AAC;TP=5.9e-01;SG=CC/AAA;SP=4.1e-01\tGT:AF:CF:GF:TF:AR:CR:GR:TR:PM\t0|0:0:21:0:0:0:4:0:0:0.0e+00\t0|1:14:7:0:0:23:8:0:0:7.1e-01\n",
			//line)==0,"Incorrect mutation line in file.");
		//printf("%s\n",line);
		count++;
	}
	mu_assert(count==1,"Wrong number of mutations output in the mutant file.");
	fclose(mut_out);
	unlink(test_snp_out);
	unlink(test_mut_out);
	unlink(test_dbg_out);
	unlink(test_no_anal_out);
	return NULL;
}

char *test_algos_calculate_per_base_normal_contamination(){
	float norm = 0.5;
	int norm_copy_no = 4;
	int tum_copy_no = 6;
	set_norm_contam(norm);
	mu_assert(get_norm_contam()==norm,"Wrong normal contamination set.");
	long double got = algos_calculate_per_base_normal_contamination(norm_copy_no,tum_copy_no);
	long double exp = (long double)2/(long double)5;
	mu_assert(got==exp,"Wrong normal contamination calculated.");
	return NULL;
}

int test_finalise(){
	genotype_t *som_tum2 = NULL;
	combined_genotype_t *ref = NULL;
	combined_genotype_t *het = NULL;
	combined_genotype_t *hom = NULL;
	combined_genotype_t *som = NULL;
	combined_genotype_t *som2 = NULL;

	genotype_t *ref_norm = NULL;
	genotype_t *ref_tum = NULL;
	genotype_t *som_tum = NULL;
	genotype_t *het_norm = NULL;
	genotype_t *het_tum = NULL;
	genotype_t *hom_norm = NULL;
	genotype_t *hom_tum = NULL;

	combined_genotype_t **het_snp_genotypes = NULL;
	combined_genotype_t **hom_snp_genotypes = NULL;
	combined_genotype_t **somatic_genotypes = NULL;

	estep_position_t *pos = malloc(sizeof(estep_position_t));
	genotype_store_t *genos = malloc(sizeof(genotype_store_t));
	check_mem(pos);
	check_mem(pos);

	char ref_base = 'C';
	char mut_base = 'T';

	int het_count = 1;
	int hom_count = 1;
	int somatic_count = 1;
	ref = malloc(sizeof(combined_genotype_t));
	check_mem(ref);
	het = malloc(sizeof(combined_genotype_t));
	check_mem(het);
	hom = malloc(sizeof(combined_genotype_t));
	check_mem(hom);
	som = malloc(sizeof(combined_genotype_t));
	check_mem(som);
	som2 = malloc(sizeof(combined_genotype_t));
	check_mem(som2);

	ref_norm = genotype_init_genotype();
	genotype_set_base_count(ref_norm, ref_base, 2);
	ref_tum = genotype_init_genotype();
	genotype_set_base_count(ref_tum, ref_base, 2);
	som_tum = genotype_init_genotype();
	genotype_set_base_count(som_tum, mut_base, 2);
	het_norm = genotype_init_genotype();
	genotype_set_base_count(het_norm, ref_base, 1);
	genotype_set_base_count(het_norm, mut_base, 1);
	het_tum = genotype_init_genotype();
	genotype_set_base_count(het_tum, ref_base, 1);
	genotype_set_base_count(het_tum, mut_base, 1);
	hom_norm = genotype_init_genotype();
	genotype_set_base_count(hom_norm, mut_base, 2);
	hom_tum = genotype_init_genotype();
	genotype_set_base_count(hom_tum, mut_base, 2);

	ref->norm_geno = ref_norm;
	ref->tum_geno = ref_tum;
	ref->prob = logl(0.002);

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = logl(0.5);

	het->norm_geno = het_norm;
	het->tum_geno = het_tum;
	het->prob = logl(0.00000125);

	hom->norm_geno = hom_norm;
	hom->tum_geno = hom_tum;
	hom->prob = logl(0.0000125);

	het_snp_genotypes = malloc(sizeof(combined_genotype_t *) * het_count);
	check_mem(het_snp_genotypes);
	het_snp_genotypes[0] = het;
	hom_snp_genotypes = malloc(sizeof(combined_genotype_t *) * hom_count);
	check_mem(hom_snp_genotypes);
	hom_snp_genotypes[0] = hom;
	somatic_genotypes = malloc(sizeof(combined_genotype_t *) * somatic_count);
	check_mem(somatic_genotypes);
	somatic_genotypes[0] = som;

	genos->ref_genotype = ref;
	genos->het_snp_genotypes = het_snp_genotypes;
	genos->hom_snp_genotypes = hom_snp_genotypes;
	genos->somatic_genotypes = somatic_genotypes;
	genos->het_count = het_count;
	genos->hom_count = hom_count;
	genos->somatic_count = somatic_count;
	genos->tum_max = 1;
	genos->norm_max = 1;
	genos->total_max = 1;

	pos->genos = genos;

	long double norm_factor_max = logl(0.5);

	finalise_probabilities_and_find_top_prob(pos,norm_factor_max);

	check(pos->top_geno==som,"Wrong top geno.");
	check(genotype_equals(pos->top_geno->norm_geno,som->norm_geno) == 1,"Top geno has incorrect normal genotype");
	check(genotype_equals(pos->top_geno->tum_geno,som->tum_geno) == 1,"Top geno has incorrect tumour genotype");
	check(pos->top_geno->prob==som->prob,"Wrong top geno prob.");
	check(pos->sec_geno==ref,"Wrong sec geno.");
	check(genotype_equals(pos->sec_geno->norm_geno,ref->norm_geno) == 1,"Sec geno has incorrect normal genotype");
	check(genotype_equals(pos->sec_geno->tum_geno,ref->tum_geno) == 1,"Sec geno has incorrect tumour genotype");
	check(pos->sec_geno->prob==ref->prob,"Wrong sec geno prob.");

	//Now with 2 somatics to ensure correct sum etc.
	somatic_count = 2;

	som_tum2 = genotype_init_genotype();
	genotype_set_base_count(som_tum2, mut_base, 1);
	genotype_set_base_count(som_tum2, ref_base, 1);

	ref->norm_geno = ref_norm;
	ref->tum_geno = ref_tum;
	ref->prob = logl(0.00000125);

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = logl(0.5);

	som2->norm_geno = ref_norm;
	som2->tum_geno = som_tum;
	som2->prob = logl(0.4);

	het->norm_geno = het_norm;
	het->tum_geno = het_tum;
	het->prob = logl(0.00000125);

	hom->norm_geno = hom_norm;
	hom->tum_geno = hom_tum;
	hom->prob = logl(0.00000125);

	het_snp_genotypes[0] = het;
	hom_snp_genotypes[0] = hom;
	free(somatic_genotypes);
	somatic_genotypes = malloc(sizeof(combined_genotype_t *) * somatic_count);
	check_mem(somatic_genotypes);
	somatic_genotypes[0] = som;
	somatic_genotypes[1] = som2;

	genos->ref_genotype = ref;
	genos->het_snp_genotypes = het_snp_genotypes;
	genos->hom_snp_genotypes = hom_snp_genotypes;
	genos->somatic_genotypes = somatic_genotypes;
	genos->het_count = het_count;
	genos->hom_count = hom_count;
	genos->somatic_count = somatic_count;
	genos->tum_max = 2;
	genos->norm_max = 1;
	genos->total_max = 2;

	pos->genos = genos;

	norm_factor_max = logl(0.5);

	pos->total_snp_prob=0;
	pos->total_mut_prob=0;

	finalise_probabilities_and_find_top_prob(pos,norm_factor_max);

	check(pos->top_geno==som,"Wrong top geno.");
	check(genotype_equals(pos->top_geno->norm_geno,som->norm_geno) == 1,"Top geno has incorrect normal genotype");
	check(genotype_equals(pos->top_geno->tum_geno,som->tum_geno) == 1,"Top geno has incorrect tumour genotype");
	check(pos->top_geno->prob==som->prob,"Wrong top geno prob.");
	check(pos->sec_geno==som2,"Wrong sec geno.");
	check(genotype_equals(pos->sec_geno->norm_geno,som2->norm_geno) == 1,"Sec geno has incorrect normal genotype");
	check(genotype_equals(pos->sec_geno->tum_geno,som2->tum_geno) == 1,"Sec geno has incorrect tumour genotype");
	check(pos->sec_geno->prob==som2->prob,"Wrong sec geno prob.");

	long double exp_snp = 0.6;
	long double exp_mut = 1.88;

	check(abs(pos->total_snp_prob-exp_snp) == 0,"Wrong total SNP prob %Le != %Le\n",pos->total_snp_prob,exp_snp);
	check(abs(pos->total_mut_prob-exp_mut) == 0,"Wrong total mut prob %Le != %Le\n",pos->total_mut_prob,exp_mut);

	if(het_snp_genotypes) free(het_snp_genotypes);
	if(hom_snp_genotypes) free(hom_snp_genotypes);
	if(somatic_genotypes) free(somatic_genotypes);
	if(ref_norm) free(ref_norm);
	if(ref_tum) free(ref_tum);
	if(som_tum) free(som_tum);
	if(som_tum2) free(som_tum2);
	if(het_norm) free(het_norm);
	if(het_tum) free(het_tum);
	if(hom_norm) free(hom_norm);
	if(hom_tum) free(hom_tum);
	if(ref) free(ref);
	if(het) free(het);
	if(hom) free(hom);
	if(som) free(som);
	if(som2) free(som2);
	if(genos) free(genos);
	if(pos) free(pos);
	return 0;
error:
	if(het_snp_genotypes) free(het_snp_genotypes);
	if(hom_snp_genotypes) free(hom_snp_genotypes);
	if(somatic_genotypes) free(somatic_genotypes);
	if(ref_norm) free(ref_norm);
	if(ref_tum) free(ref_tum);
	if(som_tum) free(som_tum);
	if(som_tum2) free(som_tum2);
	if(het_norm) free(het_norm);
	if(het_tum) free(het_tum);
	if(hom_norm) free(hom_norm);
	if(hom_tum) free(hom_tum);
	if(ref) free(ref);
	if(het) free(het);
	if(hom) free(hom);
	if(som) free(som);
	if(som2) free(som2);
	if(genos) free(genos);
	if(pos) free(pos);
	return -1;
}

char *test_finalise_probabilities_and_find_top_prob(){
	mu_assert(test_finalise() == 0,"Error testing finalize method");
	return NULL;
}

int test_per_read_estep(){
	int ref_base_idx = 1;
	unsigned char cbase = 'T';
	long double base_norm_contam = 0.4;

	long double ldone = (long double) logl(0.000552201);
	long double ldtwo = (long double) logl(0.000404538);

	combined_genotype_t *ref = NULL;
	combined_genotype_t *het = NULL;
	combined_genotype_t *hom = NULL;
	combined_genotype_t *som = NULL;
	combined_genotype_t *het_norms = NULL;

	genotype_t *ref_norm = NULL;
	genotype_t *ref_tum = NULL;
	genotype_t *som_tum = NULL;
	genotype_t *het_norm = NULL;
	genotype_t *het_tum = NULL;
	genotype_t *hom_norm = NULL;
	genotype_t *hom_tum = NULL;

	combined_genotype_t **het_snp_genotypes = NULL;
	combined_genotype_t **hom_snp_genotypes = NULL;
	combined_genotype_t **somatic_genotypes = NULL;
	combined_genotype_t **het_norm_genotypes = NULL;

	read_pos_t *norm_read = malloc(sizeof(read_pos_t));
	norm_read->ref_base_probs[0] = &ldone;
	norm_read->ref_base_probs[1] = &ldone;
	norm_read->ref_base_probs[2] = &ldtwo;
	norm_read->ref_base_probs[3] = &ldone;
	norm_read->called_base = bam_nt16_table[cbase];
	norm_read->normal = 1;
	check_mem(norm_read);
	read_pos_t *tum_read = malloc(sizeof(read_pos_t));
	tum_read->called_base = bam_nt16_table[cbase];
	tum_read->normal = 0;
	tum_read->ref_base_probs[0] = &ldone;
	tum_read->ref_base_probs[1] = &ldone;
	tum_read->ref_base_probs[2] = &ldtwo;
	tum_read->ref_base_probs[3] = &ldone;

	genotype_store_t *genos = malloc(sizeof(genotype_store_t));
	check_mem(genos);

	char ref_base = 'C';
	char mut_base = 'T';

	int het_count = 1;
	int hom_count = 1;
	int somatic_count = 1;
	int het_norm_count = 1;
	ref = malloc(sizeof(combined_genotype_t));
	check_mem(ref);
	het = malloc(sizeof(combined_genotype_t));
	check_mem(het);
	hom = malloc(sizeof(combined_genotype_t));
	check_mem(hom);
	som = malloc(sizeof(combined_genotype_t));
	check_mem(som);
	het_norms = malloc(sizeof(combined_genotype_t));
	check_mem(het_norms);

	ref_norm = genotype_init_genotype();
	genotype_set_base_count(ref_norm, ref_base, 2);
	ref_norm->var_base_idx = 1;
	ref_norm->var_base_prop = 0;
	ref_tum = genotype_init_genotype();
	genotype_set_base_count(ref_tum, ref_base, 2);
	ref_tum->var_base_idx = 1;
	ref_tum->var_base_prop = 0;
	som_tum = genotype_init_genotype();
	genotype_set_base_count(som_tum, mut_base, 2);
	som_tum->var_base_idx = 3;
	som_tum->var_base_prop = 1;
	het_norm = genotype_init_genotype();
	genotype_set_base_count(het_norm, ref_base, 1);
	genotype_set_base_count(het_norm, mut_base, 1);
	het_norm->var_base_idx = 3;
	het_norm->var_base_prop = 0.5;
	het_tum = genotype_init_genotype();
	genotype_set_base_count(het_tum, ref_base, 1);
	genotype_set_base_count(het_tum, mut_base, 1);
	het_tum->var_base_idx = 3;
	het_tum->var_base_prop = 0.5;
	hom_norm = genotype_init_genotype();
	genotype_set_base_count(hom_norm, mut_base, 2);
	hom_norm->var_base_idx = 3;
	hom_norm->var_base_prop = 1;
	hom_tum = genotype_init_genotype();
	genotype_set_base_count(hom_tum, mut_base, 2);
	hom_tum->var_base_idx = 3;
	hom_tum->var_base_prop = 1;

	ref->norm_geno = ref_norm;
	ref->tum_geno = ref_tum;
	ref->prob = 0.0;

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = 0.0;

	het->norm_geno = het_norm;
	het->tum_geno = het_tum;
	het->prob = 0.0;

	hom->norm_geno = hom_norm;
	hom->tum_geno = hom_tum;
	hom->prob = 0.0;

	het_norms->norm_geno = het_norm;
	het_norms->prob = 0.0;

	het_snp_genotypes = malloc(sizeof(combined_genotype_t *) * het_count);
	check_mem(het_snp_genotypes);
	het_snp_genotypes[0] = het;
	hom_snp_genotypes = malloc(sizeof(combined_genotype_t *) * hom_count);
	check_mem(hom_snp_genotypes);
	hom_snp_genotypes[0] = hom;
	somatic_genotypes = malloc(sizeof(combined_genotype_t *) * somatic_count);
	check_mem(somatic_genotypes);
	somatic_genotypes[0] = som;
	het_norm_genotypes = malloc(sizeof(combined_genotype_t *) * het_norm_count);
	check_mem(het_norm_genotypes);
	het_norm_genotypes[0] = het_norms;

	genos->ref_genotype = ref;
	genos->het_snp_genotypes = het_snp_genotypes;
	genos->het_snp_norm_genotypes = het_norm_genotypes;
	genos->hom_snp_genotypes = hom_snp_genotypes;
	genos->somatic_genotypes = somatic_genotypes;
	genos->het_count = het_count;
	genos->hom_count = hom_count;
	genos->somatic_count = somatic_count;
	genos->het_norm_count = het_norm_count;
	genos->ref_geno_norm_prob = 0.0;
	genos->ref_geno_tum_prob = 0.0;
	genos->tum_max = 1;
	genos->norm_max = 1;
	genos->total_max = 1;

	check_mem(tum_read);
	int chk = algos_run_per_read_estep_maths(genos, norm_read, ref_base_idx, base_norm_contam);
	check(chk==0,"Normal read failed read estep maths.");
	chk = 0;
	chk = algos_run_per_read_estep_maths(genos, tum_read, ref_base_idx, base_norm_contam);
	check(chk==0,"Tumour read failed read estep maths.");

	check(abs(som->prob - (-7.501598e+00))==0,"Wrong somatic probability.");
	check(abs(het->prob - (-7.64117748830039))==0,"Wrong het_tum probability.");
	check(abs(het_norms->prob - (-7.64117748830039))==0,"Wrong het_norm probability.");
	check(abs(hom->prob - (logl(0.00040453)*2))==0,"Wrong homozygous probability.");
	check(abs(genos->ref_geno_norm_prob - logl(0.000552201)==0),"Incorrect ref_geno_norm_prob.");
	check(abs(genos->ref_geno_tum_prob - logl(0.00040453))==0,"Incorrect ref_geno_tum_prob.");

	if(het_snp_genotypes) free(het_snp_genotypes);
	if(hom_snp_genotypes) free(hom_snp_genotypes);
	if(somatic_genotypes) free(somatic_genotypes);
	if(ref_norm) free(ref_norm);
	if(ref_tum) free(ref_tum);
	if(som_tum) free(som_tum);
	if(het_norm) free(het_norm);
	if(het_tum) free(het_tum);
	if(hom_norm) free(hom_norm);
	if(hom_tum) free(hom_tum);
	if(ref) free(ref);
	if(het) free(het);
	if(hom) free(hom);
	if(som) free(som);
	if(genos) free(genos);
	if(norm_read) free(norm_read);
	if(tum_read) free(tum_read);
	return 0;
error:
	if(het_snp_genotypes) free(het_snp_genotypes);
	if(hom_snp_genotypes) free(hom_snp_genotypes);
	if(somatic_genotypes) free(somatic_genotypes);
	if(ref_norm) free(ref_norm);
	if(ref_tum) free(ref_tum);
	if(som_tum) free(som_tum);
	if(het_norm) free(het_norm);
	if(het_tum) free(het_tum);
	if(hom_norm) free(hom_norm);
	if(hom_tum) free(hom_tum);
	if(ref) free(ref);
	if(het) free(het);
	if(hom) free(hom);
	if(som) free(som);
	if(genos) free(genos);
	if(norm_read) free(norm_read);
	if(tum_read) free(tum_read);
	return -1;
}

char *test_algos_run_per_read_estep_maths(){
	float ref_bias = 0.95;
	set_ref_bias(ref_bias);
	mu_assert(get_ref_bias()==ref_bias,"Wrong reference bias set.");
	mu_assert(test_per_read_estep()==0,"Error testing per read estep.");
	return NULL;
}

int test_estep_pos(){
	int ref_base_idx = 1;
	unsigned char cbase = 'T';
	long double base_norm_contam = 0.4;

	long double ldone =(long double) logl(0.000552201);
	long double ldtwo =(long double) logl(0.000404538);

	genotype_store_t *genos = NULL;
	estep_position_t *pos = NULL;

	combined_genotype_t *ref = NULL;
	combined_genotype_t *het = NULL;
	combined_genotype_t *hom = NULL;
	combined_genotype_t *som = NULL;
	combined_genotype_t *het_norms = NULL;

	genotype_t *ref_norm = NULL;
	genotype_t *ref_tum = NULL;
	genotype_t *som_tum = NULL;
	genotype_t *het_norm = NULL;
	genotype_t *het_tum = NULL;
	genotype_t *hom_norm = NULL;
	genotype_t *hom_tum = NULL;

	combined_genotype_t **het_snp_genotypes = NULL;
	combined_genotype_t **hom_snp_genotypes = NULL;
	combined_genotype_t **somatic_genotypes = NULL;
	combined_genotype_t **het_norm_genotypes = NULL;


	read_pos_t *norm_read = malloc(sizeof(read_pos_t));
	read_pos_t *tum_read = malloc(sizeof(read_pos_t));
	check_mem(tum_read);
	check_mem(norm_read);
	norm_read->ref_base_probs[0] = &ldone;
	norm_read->ref_base_probs[1] = &ldone;
	norm_read->ref_base_probs[2] = &ldtwo;
	norm_read->ref_base_probs[3] = &ldone;
	norm_read->called_base = bam_nt16_table[cbase];
	norm_read->normal = 1;

	tum_read->called_base = bam_nt16_table[cbase];
	tum_read->normal = 0;
	tum_read->ref_base_probs[0] = &ldone;
	tum_read->ref_base_probs[1] = &ldone;
	tum_read->ref_base_probs[2] = &ldtwo;
	tum_read->ref_base_probs[3] = &ldone;

	genos = malloc(sizeof(genotype_store_t));
	check_mem(genos);

	pos = malloc(sizeof(estep_position_t));
	check_mem(pos);


	char ref_base = 'C';
	char mut_base = 'T';

	int het_count = 1;
	int hom_count = 1;
	int somatic_count = 1;
	int het_norm_count = 1;
	ref = malloc(sizeof(combined_genotype_t));
	check_mem(ref);
	het = malloc(sizeof(combined_genotype_t));
	check_mem(het);
	hom = malloc(sizeof(combined_genotype_t));
	check_mem(hom);
	som = malloc(sizeof(combined_genotype_t));
	check_mem(som);
	het_norms = malloc(sizeof(combined_genotype_t));
	check_mem(het_norms);

	ref_norm = genotype_init_genotype();
	genotype_set_base_count(ref_norm, ref_base, 2);
	ref_norm->var_base_idx = 1;
	ref_norm->var_base_prop = 0;
	ref_tum = genotype_init_genotype();
	genotype_set_base_count(ref_tum, ref_base, 2);
	ref_tum->var_base_idx = 1;
	ref_tum->var_base_prop = 0;
	som_tum = genotype_init_genotype();
	genotype_set_base_count(som_tum, mut_base, 2);
	som_tum->var_base_idx = 3;
	som_tum->var_base_prop = 1;
	het_norm = genotype_init_genotype();
	genotype_set_base_count(het_norm, ref_base, 1);
	genotype_set_base_count(het_norm, mut_base, 1);
	het_norm->var_base_idx = 3;
	het_norm->var_base_prop = 0.5;
	het_tum = genotype_init_genotype();
	genotype_set_base_count(het_tum, ref_base, 1);
	genotype_set_base_count(het_tum, mut_base, 1);
	het_tum->var_base_idx = 3;
	het_tum->var_base_prop = 0.5;
	hom_norm = genotype_init_genotype();
	genotype_set_base_count(hom_norm, mut_base, 2);
	hom_norm->var_base_idx = 3;
	hom_norm->var_base_prop = 1;
	hom_tum = genotype_init_genotype();
	genotype_set_base_count(hom_tum, mut_base, 2);
	hom_tum->var_base_idx = 3;
	hom_tum->var_base_prop = 1;

	ref->norm_geno = ref_norm;
	ref->tum_geno = ref_tum;
	ref->prob = 0.0;

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = 0.0;

	het->norm_geno = het_norm;
	het->tum_geno = het_tum;
	het->prob = 0.0;

	hom->norm_geno = hom_norm;
	hom->tum_geno = hom_tum;
	hom->prob = 0.0;

	het_norms->norm_geno = het_norm;
	het_norms->prob = 0.0;

	het_snp_genotypes = malloc(sizeof(combined_genotype_t *) * het_count);
	check_mem(het_snp_genotypes);
	het_snp_genotypes[0] = het;
	hom_snp_genotypes = malloc(sizeof(combined_genotype_t *) * hom_count);
	check_mem(hom_snp_genotypes);
	hom_snp_genotypes[0] = hom;
	somatic_genotypes = malloc(sizeof(combined_genotype_t *) * somatic_count);
	check_mem(somatic_genotypes);
	somatic_genotypes[0] = som;
	het_norm_genotypes = malloc(sizeof(combined_genotype_t *) * het_norm_count);
	check_mem(het_norm_genotypes);
	het_norm_genotypes[0] = het_norms;

	genos->ref_genotype = ref;
	genos->het_snp_genotypes = het_snp_genotypes;
	genos->het_snp_norm_genotypes = het_norm_genotypes;
	genos->hom_snp_genotypes = hom_snp_genotypes;
	genos->somatic_genotypes = somatic_genotypes;
	genos->het_count = het_count;
	genos->hom_count = hom_count;
	genos->somatic_count = somatic_count;
	genos->het_norm_count = het_norm_count;
	genos->ref_geno_norm_prob = 0.0;
	genos->ref_geno_tum_prob = 0.0;

	genos->tum_max = 1;
	genos->norm_max = 1;
	genos->total_max = 1;

	pos->genos = genos;

	check_mem(tum_read);
	int chk = algos_run_per_read_estep_maths(genos, norm_read, ref_base_idx, base_norm_contam);
	check(chk==0,"Normal read failed read estep maths.");
	chk = 0;
	chk = algos_run_per_read_estep_maths(genos, tum_read, ref_base_idx, base_norm_contam);
	check(chk==0,"Tumour read failed read estep maths.");

	check(abs(som->prob - (-7.501598e+00))==0,"Wrong somatic probability.");
	check(abs(het->prob - (-7.64117748830039))==0,"Wrong het_tum probability.");
	check(abs(het_norms->prob - (-7.64117748830039))==0,"Wrong het_norm probability.");
	check(abs(hom->prob - (logl(0.00040453)*2))==0,"Wrong homozygous probability.");
	check(abs(genos->ref_geno_norm_prob - logl(0.000552201)==0),"Incorrect ref_geno_norm_prob.");
	check(abs(genos->ref_geno_tum_prob - logl(0.00040453))==0,"Incorrect ref_geno_tum_prob.");

	algos_run_per_position_estep_maths(pos);

	long double norm_fact = exp((logl(0.000552201) + logl(0.00040453)) + logl(1-0.001-0.00006))
									+ exp(logl(0.00006) + logl(0.000552201) + (-7.501598e+00))
									+ exp((logl(0.00040453)*2))
									+ exp((logl(0.001) - log(2)) + (-7.64117748830039) + (-7.64117748830039));

	check(abs(ref->prob - ((exp((logl(0.000552201) + logl(0.00040453)) + logl(1-0.001-0.00006)))/norm_fact))==0,"Wrong reference probability.");
	check(abs(som->prob - (exp(logl(0.00006) + logl(0.000552201) + (-7.501598e+00))/norm_fact))==0,"Wrong somatic probability.");
	check(abs(het->prob - (exp((logl(0.001) - log(2)) + (-7.64117748830039) + (-7.64117748830039))/norm_fact))==0,"Wrong het probability.");
	check(abs(hom->prob - (exp((logl(0.00040453)*2))/norm_fact))==0,"Wrong hom probability.");

	check(genotype_equals(hom->norm_geno,pos->top_geno->norm_geno)==1,"Wrong first genotype.");
	check(genotype_equals(hom->tum_geno,pos->top_geno->tum_geno)==1,"Wrong first genotype.");
	check(genotype_equals(ref->norm_geno,pos->sec_geno->norm_geno)==1,"Wrong second genotype.");
	check(genotype_equals(ref->tum_geno,pos->sec_geno->tum_geno)==1,"Wrong second genotype.");

	check(abs(pos->total_snp_prob - (long double)(2.500625e-04 + 5.001250e-01))== 0,"Wrong total SNP probability.");
	check(abs(pos->total_mut_prob-(long double)3.000750e-05)==0,"Wrong total mutation probability.");

	if(het_snp_genotypes) free(het_snp_genotypes);
	if(hom_snp_genotypes) free(hom_snp_genotypes);
	if(somatic_genotypes) free(somatic_genotypes);
	if(ref_norm) free(ref_norm);
	if(ref_tum) free(ref_tum);
	if(som_tum) free(som_tum);
	if(het_norm) free(het_norm);
	if(het_tum) free(het_tum);
	if(hom_norm) free(hom_norm);
	if(hom_tum) free(hom_tum);
	if(ref) free(ref);
	if(het) free(het);
	if(hom) free(hom);
	if(som) free(som);
	if(genos) free(genos);
	if(norm_read) free(norm_read);
	if(tum_read) free(tum_read);
	if(pos) free(pos);
	return 0;
error:
	if(het_snp_genotypes) free(het_snp_genotypes);
	if(hom_snp_genotypes) free(hom_snp_genotypes);
	if(somatic_genotypes) free(somatic_genotypes);
	if(ref_norm) free(ref_norm);
	if(ref_tum) free(ref_tum);
	if(som_tum) free(som_tum);
	if(het_norm) free(het_norm);
	if(het_tum) free(het_tum);
	if(hom_norm) free(hom_norm);
	if(hom_tum) free(hom_tum);
	if(ref) free(ref);
	if(het) free(het);
	if(hom) free(hom);
	if(som) free(som);
	if(genos) free(genos);
	if(norm_read) free(norm_read);
	if(tum_read) free(tum_read);
	if(pos) free(pos);
	return -1;
}

char *test_algos_run_per_position_estep_maths(){
	float ref_bias = 0.95;
	float snp_prob = 0.001;
	float mut_prob = 0.00006;
	set_ref_bias(ref_bias);
	mu_assert(get_ref_bias()==ref_bias,"Wrong reference bias set.");
	set_prior_mut_prob(mut_prob);
	mu_assert(get_prior_mut_prob()==mut_prob,"Wrong mut prior set.");
	set_prior_snp_prob(snp_prob);
	mu_assert(get_prior_snp_prob()==snp_prob,"Wrong SNP prior set.");
	mu_assert(test_estep_pos()==0,"Error testing per position estep.");
	return NULL;
}

char *test_algos_check_var_position_alleles(){
	char *type = "TEST_SOMATIC";
	char *chr_name = "TEST_CHR";
	int ref_pos = 10;

	combined_genotype_t *ref = NULL;
	combined_genotype_t *som = NULL;

	genotype_t *ref_norm = NULL;
	genotype_t *ref_tum = NULL;
	genotype_t *som_tum = NULL;
	genotype_t *alt_som_tum = NULL;

	genotype_store_t *genos = malloc(sizeof(genotype_store_t));
	estep_position_t *pos = malloc(sizeof(estep_position_t));
	check_mem(genos);
	check_mem(pos);


	char ref_base = 'C';
	char mut_base = 'T';
	char wrong_mut_base = 'G';

	double top_prob = 0.95;
	double sec_prob = 0.05;

	pos->ref_pos = ref_pos;
	pos->ref_base = "C";

	ref = malloc(sizeof(combined_genotype_t));
	check_mem(ref);
	som = malloc(sizeof(combined_genotype_t));
	check_mem(som);

	ref_norm = genotype_init_genotype();
	genotype_set_base_count(ref_norm, ref_base, 2);
	ref_norm->var_base_idx = 1;
	ref_norm->var_base = ref_base;
	ref_norm->var_base_prop = 0;
	ref_tum = genotype_init_genotype();
	genotype_set_base_count(ref_tum, ref_base, 2);
	ref_tum->var_base_idx = 1;
	ref_tum->var_base = ref_base;
	ref_tum->var_base_prop = 0;
	som_tum = genotype_init_genotype();
	genotype_set_base_count(som_tum, mut_base, 2);
	som_tum->var_base_idx = 3;
	som_tum->var_base = mut_base;
	som_tum->var_base_prop = 1;

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = top_prob;

	ref->norm_geno = ref_norm;
	ref->tum_geno = ref_tum;
	ref->prob	= sec_prob;


	pos->norm_fwd_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->norm_fwd_cvg, 'A', 0);
	genotype_set_base_count(pos->norm_fwd_cvg, 'C', 20);
	genotype_set_base_count(pos->norm_fwd_cvg, 'G', 0);
	genotype_set_base_count(pos->norm_fwd_cvg, 'T', 0);
	pos->norm_rev_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->norm_rev_cvg, 'A', 0);
	genotype_set_base_count(pos->norm_rev_cvg, 'C', 6);
	genotype_set_base_count(pos->norm_rev_cvg, 'G', 0);
	genotype_set_base_count(pos->norm_rev_cvg, 'T', 0);
	pos->tum_fwd_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->tum_fwd_cvg, 'A', 0);
	genotype_set_base_count(pos->tum_fwd_cvg, 'C', 0);
	genotype_set_base_count(pos->tum_fwd_cvg, 'G', 0);
	genotype_set_base_count(pos->tum_fwd_cvg, 'T', 1);
	pos->tum_rev_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->tum_rev_cvg, 'A', 0);
	genotype_set_base_count(pos->tum_rev_cvg, 'C', 0);
	genotype_set_base_count(pos->tum_rev_cvg, 'G', 0);
	genotype_set_base_count(pos->tum_rev_cvg, 'T', 2);

	pos->top_geno = som;
	pos->sec_geno = ref;

	int result = algos_check_var_position_alleles(pos, chr_name, type,0);
	mu_assert(result == 1,"Somatic, homozygous counts check 1");
	mu_assert(pos->top_geno->prob==top_prob,"Top probability unchanged");
	mu_assert(pos->sec_geno->prob==sec_prob,"Second probability unchanged");


	//We just reverse top and second genos to check top genotype is normal.
	pos->top_geno = ref;
	pos->sec_geno = som;
	result = algos_check_var_position_alleles(pos, chr_name, type,0);
	mu_assert(result == 1,"Somatic, homozygous counts check 2");
	mu_assert(pos->top_geno->prob==top_prob,"Top probability changed as normal top genotype");
	mu_assert(pos->sec_geno->prob==sec_prob,"Second probability changed as normal top genotype");


	//TODO top geno and second geno swap
	alt_som_tum = genotype_init_genotype();
	genotype_set_base_count(alt_som_tum, wrong_mut_base, 2);
	alt_som_tum->var_base_idx = 2;
	alt_som_tum->var_base = wrong_mut_base;
	alt_som_tum->var_base_prop = 1;

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = sec_prob;

	ref->norm_geno = ref_norm;
	ref->tum_geno = alt_som_tum;
	ref->prob	= top_prob;

	pos->top_geno = ref;
	pos->sec_geno = som;

	result = algos_check_var_position_alleles(pos, chr_name, type,0);

	mu_assert(result == 1,"Somatic, homozygous counts check 3");
	mu_assert(pos->top_geno->prob==sec_prob,"Top probability changed as erroneous top genotype");
	mu_assert(pos->sec_geno->prob==top_prob,"Second probability changed as erroneous top genotype");

	//TODO complete failure.
	pos->top_geno = ref;
	pos->sec_geno = ref;

	result = algos_check_var_position_alleles(pos, chr_name, type,0);
	mu_assert(result == 0,"Genotype check fail");
	return NULL;
error:
	return "Failed mem check";
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_algos_mstep_read_position);
   mu_run_test(test_algos_mstep_read_position_two);
   mu_run_test(test_finalise_probabilities_and_find_top_prob);
   mu_run_test(test_algos_calculate_per_base_normal_contamination);
   mu_run_test(test_algos_run_per_read_estep_maths);
   mu_run_test(test_algos_run_per_position_estep_maths);
   mu_run_test(test_algos_estep_read_position_real_data_no_analysis);
   mu_run_test(test_algos_estep_read_position);
   mu_run_test(test_algos_check_var_position_alleles);
   return NULL;
}

RUN_TESTS(all_tests);
