/**   LICENSE
* Copyright (c) 2014-2015 Genome Research Ltd.
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

#include "minunit.h"
#include <output.h>
#include <string.h>
#include <stdio.h>
#include <algos.h>
#include <covs_access.h>
#include <bam_access.h>
#include <alg_bean.h>
#include <genotype.h>
#include <unistd.h>
#include <time.h>

char *mut_norm = "testData/testing_wt.bam";
char *mut_tum = "testData/testing_mt.bam";
char *mut_norm_cram = "testData/testing_wt.cram";
char *mut_tum_cram = "testData/testing_mt.cram";
char *out_test_vcf = "testData/test_out.vcf";
char *test_fai_out = TEST_REF;

char *test_output_generate_info_lines(){
	char *result = output_generate_info_lines();
	char exp[2500];
	strcpy(exp,"");
	strcat(exp,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
	strcat(exp,"##INFO=<ID=MP,Number=1,Type=Float,Description=\"Sum of CaVEMan somatic genotype probabilities\">\n");
	strcat(exp,"##INFO=<ID=GP,Number=1,Type=Float,Description=\"Sum of CaVEMan germline genotype probabilities\">\n");
	strcat(exp,"##INFO=<ID=TG,Number=1,Type=String,Description=\"Most probable genotype as called by CaVEMan\">\n");
	strcat(exp,"##INFO=<ID=TP,Number=1,Type=Float,Description=\"Probability of most probable genotype as called by CaVEMan\">\n");
	strcat(exp,"##INFO=<ID=SG,Number=1,Type=String,Description=\"2nd most probable genotype as called by CaVEMan\">\n");
	strcat(exp,"##INFO=<ID=SP,Number=1,Type=Float,Description=\"Probability of 2nd most probable genotype as called by CaVEMan\">\n");
	strcat(exp,"##INFO=<ID=DS,Number=.,Type=String,Description=\"DBSnp ID of known SNP\">\n");
	mu_assert(strcmp(result,exp)==0,"Wrong info format string.");
	return NULL;
}

char *test_output_generate_format_lines(){
	char *result = output_generate_format_lines();
	char exp[2500];
	strcpy(exp,"");
	strcat(exp,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
	strcat(exp,"##FORMAT=<ID=FAZ,Number=1,Type=Integer,Description=\"Reads presenting a A for this position, forward strand\">\n");
	strcat(exp,"##FORMAT=<ID=FCZ,Number=1,Type=Integer,Description=\"Reads presenting a C for this position, forward strand\">\n");
	strcat(exp,"##FORMAT=<ID=FGZ,Number=1,Type=Integer,Description=\"Reads presenting a G for this position, forward strand\">\n");
	strcat(exp,"##FORMAT=<ID=FTZ,Number=1,Type=Integer,Description=\"Reads presenting a T for this position, forward strand\">\n");
	strcat(exp,"##FORMAT=<ID=RAZ,Number=1,Type=Integer,Description=\"Reads presenting a A for this position, reverse strand\">\n");
	strcat(exp,"##FORMAT=<ID=RCZ,Number=1,Type=Integer,Description=\"Reads presenting a C for this position, reverse strand\">\n");
	strcat(exp,"##FORMAT=<ID=RGZ,Number=1,Type=Integer,Description=\"Reads presenting a G for this position, reverse strand\">\n");
	strcat(exp,"##FORMAT=<ID=RTZ,Number=1,Type=Integer,Description=\"Reads presenting a T for this position, reverse strand\">\n");
	strcat(exp,"##FORMAT=<ID=PM,Number=1,Type=Float,Description=\"Proportion of mut allele\">\n");
	mu_assert(strcmp(result,exp)==0,"Wrong format format string.");
	return NULL;
}

int test_output_to_file_int(){
	int ref_pos = 10;

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

	FILE *out = NULL;

	genotype_store_t *genos = malloc(sizeof(genotype_store_t));
	estep_position_t *pos = malloc(sizeof(estep_position_t));
	check_mem(genos);
	check_mem(pos);


	char ref_base = 'C';
	char mut_base = 'T';
	char *chrom = "Y";

	int het_count = 1;
	int hom_count = 1;
	int somatic_count = 1;
	int het_norm_count = 1;

	pos->ref_pos = ref_pos;
	pos->ref_base = "C";

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
	ref_norm->var_base = 'C';
	ref_norm->var_base_prop = 0;
	ref_tum = genotype_init_genotype();
	genotype_set_base_count(ref_tum, ref_base, 2);
	ref_tum->var_base_idx = 1;
	ref_tum->var_base = 'C';
	ref_tum->var_base_prop = 0;
	som_tum = genotype_init_genotype();
	genotype_set_base_count(som_tum, mut_base, 2);
	som_tum->var_base_idx = 3;
	som_tum->var_base = 'T';
	som_tum->var_base_prop = 1;
	het_norm = genotype_init_genotype();
	genotype_set_base_count(het_norm, ref_base, 1);
	genotype_set_base_count(het_norm, mut_base, 1);
	het_norm->var_base_idx = 3;
	het_norm->var_base = 'T';
	het_norm->var_base_prop = 0.5;
	het_tum = genotype_init_genotype();
	genotype_set_base_count(het_tum, ref_base, 1);
	genotype_set_base_count(het_tum, mut_base, 1);
	het_tum->var_base_idx = 3;
	het_tum->var_base = 'T';
	het_tum->var_base_prop = 0.5;
	hom_norm = genotype_init_genotype();
	genotype_set_base_count(hom_norm, mut_base, 2);
	hom_norm->var_base_idx = 3;
	hom_norm->var_base = 'T';
	hom_norm->var_base_prop = 1;
	hom_tum = genotype_init_genotype();
	genotype_set_base_count(hom_tum, mut_base, 2);
	hom_tum->var_base_idx = 3;
	hom_tum->var_base = 'T';
	hom_tum->var_base_prop = 1;

	pos->norm_fwd_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->norm_fwd_cvg, 'A', 1);
	genotype_set_base_count(pos->norm_fwd_cvg, 'C', 2);
	genotype_set_base_count(pos->norm_fwd_cvg, 'G', 3);
	genotype_set_base_count(pos->norm_fwd_cvg, 'T', 4);
	pos->norm_rev_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->norm_rev_cvg, 'A', 5);
	genotype_set_base_count(pos->norm_rev_cvg, 'C', 6);
	genotype_set_base_count(pos->norm_rev_cvg, 'G', 7);
	genotype_set_base_count(pos->norm_rev_cvg, 'T', 8);
	pos->tum_fwd_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->tum_fwd_cvg, 'A', 9);
	genotype_set_base_count(pos->tum_fwd_cvg, 'C', 10);
	genotype_set_base_count(pos->tum_fwd_cvg, 'G', 11);
	genotype_set_base_count(pos->tum_fwd_cvg, 'T', 12);
	pos->tum_rev_cvg = genotype_init_genotype();
	genotype_set_base_count(pos->tum_rev_cvg, 'A', 13);
	genotype_set_base_count(pos->tum_rev_cvg, 'C', 14);
	genotype_set_base_count(pos->tum_rev_cvg, 'G', 15);
	genotype_set_base_count(pos->tum_rev_cvg, 'T', 16);

	ref->norm_geno = ref_norm;
	ref->tum_geno = ref_tum;
	ref->prob = 0.04789;

	som->norm_geno = ref_norm;
	som->tum_geno = som_tum;
	som->prob = 0.948;

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

	pos->genos = genos;

	pos->total_snp_prob = 0.076589;
	pos->total_mut_prob = 0.93678;
	pos->top_geno = som;
	pos->sec_geno = ref;

	out = fopen(out_test_vcf,"w");

	int chk = output_vcf_variant_position(pos, out, chrom);
	check(chk==0,"Error running output_vcf_variant_position.");

	fclose(out);

	char *exp = "Y\t10\t.\tC\tT\t.\t.\tDP=136;MP=9.4e-01;GP=7.7e-02;TG=CC/TT;TP=9.5e-01;SG=CC/CC;SP=4.8e-02\tGT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM\t0|0:1:2:3:4:5:6:7:8:3.3e-01\t1|1:9:10:11:12:13:14:15:16:2.8e-01\n";
	out = fopen(out_test_vcf,"r");
	char line[5000];
	int count = 0;
	while ( fgets(line,sizeof(line),out) != NULL ){
		check(count==0,"Too many lines output to mutation file.");
		check(strcmp(exp,line)==0,"Incorrect variant output in file.");
		count++;
	}
	fclose(out);
	unlink(out_test_vcf);


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
	if(pos) free(pos);
	return -1;
}

char *test_output_to_file(){
	float ref_bias = 0.95;
	float snp_prob = 0.001;
	float mut_prob = 0.00006;
	set_ref_bias(ref_bias);
	mu_assert(get_ref_bias()==ref_bias,"Wrong reference bias set.");
	set_prior_mut_prob(mut_prob);
	mu_assert(get_prior_mut_prob()==mut_prob,"Wrong mut prior set.");
	set_prior_snp_prob(snp_prob);
	mu_assert(get_prior_snp_prob()==snp_prob,"Wrong SNP prior set.");
	mu_assert(test_output_to_file_int()==0,"Error testing output to file.");
	return NULL;
}

char *test_output_header_to_file(){
	FILE *out = fopen(out_test_vcf,"w");
	char *norm_protocol = "WGS";
	char *tum_protocol = "WXS";
	char *norm_plat,*tum_plat;
	norm_plat = malloc(sizeof(char)*50);
	strcpy(norm_plat,".");
	tum_plat = malloc(sizeof(char)*50);
	strcpy(tum_plat,".");
	int chk = output_vcf_header(out, mut_tum, mut_norm, test_fai_out,
																		NULL, NULL, norm_protocol, tum_protocol, norm_plat, tum_plat);
	mu_assert(chk==0,"Error running output header method.");

	fclose(out);

	out = fopen(out_test_vcf,"r");
	char exp[20000];
	strcpy(exp,"");
	strcat(exp,"##fileformat=VCFv4.1\n");
	//fileDate=20120104
	char date[50];
	time_t t = time(NULL);
	strftime(date,sizeof(date),"%Y%m%d",localtime(&t));
	char tmp[2000];
	sprintf(tmp,"##fileDate=%s\n",date);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp,"##reference=%s\n",test_fai_out);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp, "##vcfProcessLog=<InputVCF=<.>,InputVCFSource=<CaVEMan>,InpuVCFVer=<\"%s\">,InputVCFParam=<NORMAL_CONTAMINATION=%g,REF_BIAS=%g,PRIOR_MUT_RATE=%g,PRIOR_SNP_RATE=%g,SNP_CUTOFF=%g,MUT_CUTOFF=%g>>\n",CAVEMAN_VERSION,get_norm_contam(),get_ref_bias(),get_prior_mut_prob(),get_prior_snp_prob(),get_min_snp_prob(),get_min_mut_prob());
	strcat(exp,tmp);
	sprintf(tmp,"##cavemanVersion=%s\n",CAVEMAN_VERSION);
	strcat(exp,tmp);
	strcat(exp,	"##contig=<ID=1,length=249250621,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=2,length=243199373,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=3,length=198022430,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=4,length=191154276,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=5,length=180915260,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=6,length=171115067,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=7,length=159138663,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=8,length=146364022,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=9,length=141213431,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=10,length=135534747,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=11,length=135006516,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=12,length=133851895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=13,length=115169878,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=14,length=107349540,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=15,length=102531392,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=16,length=90354753,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=17,length=81195210,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=18,length=78077248,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=19,length=59128983,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=20,length=63025520,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=21,length=48129895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=22,length=51304566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=X,length=155270560,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=Y,length=59373566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=MT,length=16569,assembly=37,species=HUMAN>\n");

	strcat(exp,output_generate_info_lines());
	strcat(exp,output_generate_format_lines());

	strcat(exp,"##SAMPLE=<ID=NORMAL,Description=\"Normal\",Accession=.,Platform=HiSeq,Protocol=WGS,SampleName=NORMALb,Source=.>\n");
	strcat(exp,"##SAMPLE=<ID=TUMOUR,Description=\"Tumour\",Accession=.,Platform=HiSeq,Protocol=WXS,SampleName=TUMOURa,Source=.>\n");
	strcat(exp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR\n");

	char line[1000];
	int count = 0;
	int exp_lines = 52;
	char got[20000];
	strcpy(got,"");
	while ( fgets(line,sizeof(line),out) != NULL ){
		mu_assert(count<=exp_lines,"Too many header lines.");
		strcat(got,line);
		count++;
	}
	if(strcmp(exp,got)!=0){
		printf("exp:\n%s\n\n\ngot:\n%s\n\n",exp,got);
	}
	mu_assert(strcmp(exp,got)==0,"Header in file doesn't match.");
	fclose(out);
	unlink(out_test_vcf);

	chk=0;

	out = fopen(out_test_vcf,"w");

	char *norm_plat2 = malloc(sizeof(char)*50);
	char *tum_plat2 = malloc(sizeof(char)*50);

	strcpy(norm_plat2,"TEST");
	strcpy(tum_plat2,"TEST");
	chk = output_vcf_header(out, mut_tum, mut_norm, test_fai_out,
																		NULL, NULL, norm_protocol, tum_protocol, norm_plat2, tum_plat2);
	mu_assert(chk==0,"Error running output header method.");

	fclose(out);

	out = fopen(out_test_vcf,"r");
	strcpy(exp,"");
	strcat(exp,"##fileformat=VCFv4.1\n");
	//fileDate=20120104
	t = time(NULL);
	strftime(date,sizeof(date),"%Y%m%d",localtime(&t));
	sprintf(tmp,"##fileDate=%s\n",date);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp,"##reference=%s\n",test_fai_out);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp, "##vcfProcessLog=<InputVCF=<.>,InputVCFSource=<CaVEMan>,InpuVCFVer=<\"%s\">,InputVCFParam=<NORMAL_CONTAMINATION=%g,REF_BIAS=%g,PRIOR_MUT_RATE=%g,PRIOR_SNP_RATE=%g,SNP_CUTOFF=%g,MUT_CUTOFF=%g>>\n",CAVEMAN_VERSION,get_norm_contam(),get_ref_bias(),get_prior_mut_prob(),get_prior_snp_prob(),get_min_snp_prob(),get_min_mut_prob());
	strcat(exp,tmp);
	sprintf(tmp,"##cavemanVersion=%s\n",CAVEMAN_VERSION);
	strcat(exp,tmp);
	strcat(exp,	"##contig=<ID=1,length=249250621,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=2,length=243199373,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=3,length=198022430,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=4,length=191154276,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=5,length=180915260,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=6,length=171115067,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=7,length=159138663,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=8,length=146364022,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=9,length=141213431,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=10,length=135534747,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=11,length=135006516,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=12,length=133851895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=13,length=115169878,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=14,length=107349540,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=15,length=102531392,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=16,length=90354753,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=17,length=81195210,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=18,length=78077248,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=19,length=59128983,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=20,length=63025520,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=21,length=48129895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=22,length=51304566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=X,length=155270560,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=Y,length=59373566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=MT,length=16569,assembly=37,species=HUMAN>\n");

	strcat(exp,output_generate_info_lines());
	strcat(exp,output_generate_format_lines());

	strcat(exp,"##SAMPLE=<ID=NORMAL,Description=\"Normal\",Accession=.,Platform=TEST,Protocol=WGS,SampleName=NORMALb,Source=.>\n");
	strcat(exp,"##SAMPLE=<ID=TUMOUR,Description=\"Tumour\",Accession=.,Platform=TEST,Protocol=WXS,SampleName=TUMOURa,Source=.>\n");
	strcat(exp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR\n");

	count = 0;
	exp_lines = 52;
	strcpy(got,"");
	while ( fgets(line,sizeof(line),out) != NULL ){
		mu_assert(count<=exp_lines,"Too many header lines.");
		strcat(got,line);
		count++;
	}
	if(strcmp(exp,got)!=0){
		printf("exp:\n%s\n\n\ngot:\n%s\n\n",exp,got);
	}
	mu_assert(strcmp(exp,got)==0,"Header in file doesn't match.");
	fclose(out);
	unlink(out_test_vcf);

	return NULL;
}

char *test_output_header_to_file_cram(){
	FILE *out = fopen(out_test_vcf,"w");
	char *norm_protocol = "WGS";
	char *tum_protocol = "WXS";
	char *norm_plat,*tum_plat;
	norm_plat = malloc(sizeof(char)*50);
	strcpy(norm_plat,".");
	tum_plat = malloc(sizeof(char)*50);
	strcpy(tum_plat,".");
	int chk = output_vcf_header(out, mut_tum_cram, mut_norm_cram, test_fai_out,
																		NULL, NULL, norm_protocol, tum_protocol, norm_plat, tum_plat);
	mu_assert(chk==0,"Error running output header method.");

	fclose(out);

	out = fopen(out_test_vcf,"r");
	char exp[20000];
	strcpy(exp,"");
	strcat(exp,"##fileformat=VCFv4.1\n");
	//fileDate=20120104
	char date[50];
	time_t t = time(NULL);
	strftime(date,sizeof(date),"%Y%m%d",localtime(&t));
	char tmp[2000];
	sprintf(tmp,"##fileDate=%s\n",date);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp,"##reference=%s\n",test_fai_out);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp, "##vcfProcessLog=<InputVCF=<.>,InputVCFSource=<CaVEMan>,InpuVCFVer=<\"%s\">,InputVCFParam=<NORMAL_CONTAMINATION=%g,REF_BIAS=%g,PRIOR_MUT_RATE=%g,PRIOR_SNP_RATE=%g,SNP_CUTOFF=%g,MUT_CUTOFF=%g>>\n",CAVEMAN_VERSION,get_norm_contam(),get_ref_bias(),get_prior_mut_prob(),get_prior_snp_prob(),get_min_snp_prob(),get_min_mut_prob());
	strcat(exp,tmp);
	sprintf(tmp,"##cavemanVersion=%s\n",CAVEMAN_VERSION);
	strcat(exp,tmp);
	strcat(exp,	"##contig=<ID=1,length=249250621,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=2,length=243199373,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=3,length=198022430,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=4,length=191154276,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=5,length=180915260,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=6,length=171115067,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=7,length=159138663,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=8,length=146364022,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=9,length=141213431,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=10,length=135534747,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=11,length=135006516,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=12,length=133851895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=13,length=115169878,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=14,length=107349540,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=15,length=102531392,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=16,length=90354753,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=17,length=81195210,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=18,length=78077248,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=19,length=59128983,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=20,length=63025520,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=21,length=48129895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=22,length=51304566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=X,length=155270560,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=Y,length=59373566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=MT,length=16569,assembly=37,species=HUMAN>\n");

	strcat(exp,output_generate_info_lines());
	strcat(exp,output_generate_format_lines());

	strcat(exp,"##SAMPLE=<ID=NORMAL,Description=\"Normal\",Accession=.,Platform=HiSeq,Protocol=WGS,SampleName=NORMALb,Source=.>\n");
	strcat(exp,"##SAMPLE=<ID=TUMOUR,Description=\"Tumour\",Accession=.,Platform=HiSeq,Protocol=WXS,SampleName=TUMOURa,Source=.>\n");
	strcat(exp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR\n");

	char line[1000];
	int count = 0;
	int exp_lines = 52;
	char got[20000];
	strcpy(got,"");
	while ( fgets(line,sizeof(line),out) != NULL ){
		mu_assert(count<=exp_lines,"Too many header lines.");
		strcat(got,line);
		count++;
	}
	if(strcmp(exp,got)!=0){
		printf("exp:\n%s\n\n\ngot:\n%s\n\n",exp,got);
	}
	mu_assert(strcmp(exp,got)==0,"Header in file doesn't match.");
	fclose(out);
	unlink(out_test_vcf);

	chk=0;

	out = fopen(out_test_vcf,"w");

	char *norm_plat2 = malloc(sizeof(char)*50);
	char *tum_plat2 = malloc(sizeof(char)*50);

	strcpy(norm_plat2,"TEST");
	strcpy(tum_plat2,"TEST");
	chk = output_vcf_header(out, mut_tum_cram, mut_norm_cram, test_fai_out,
																		NULL, NULL, norm_protocol, tum_protocol, norm_plat2, tum_plat2);
	mu_assert(chk==0,"Error running output header method.");

	fclose(out);

	out = fopen(out_test_vcf,"r");
	strcpy(exp,"");
	strcat(exp,"##fileformat=VCFv4.1\n");
	//fileDate=20120104
	t = time(NULL);
	strftime(date,sizeof(date),"%Y%m%d",localtime(&t));
	sprintf(tmp,"##fileDate=%s\n",date);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp,"##reference=%s\n",test_fai_out);
	strcat(exp,tmp);
	strcpy(tmp,"");
	sprintf(tmp, "##vcfProcessLog=<InputVCF=<.>,InputVCFSource=<CaVEMan>,InpuVCFVer=<\"%s\">,InputVCFParam=<NORMAL_CONTAMINATION=%g,REF_BIAS=%g,PRIOR_MUT_RATE=%g,PRIOR_SNP_RATE=%g,SNP_CUTOFF=%g,MUT_CUTOFF=%g>>\n",CAVEMAN_VERSION,get_norm_contam(),get_ref_bias(),get_prior_mut_prob(),get_prior_snp_prob(),get_min_snp_prob(),get_min_mut_prob());
	strcat(exp,tmp);
	sprintf(tmp,"##cavemanVersion=%s\n",CAVEMAN_VERSION);
	strcat(exp,tmp);
	strcat(exp,	"##contig=<ID=1,length=249250621,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=2,length=243199373,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=3,length=198022430,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=4,length=191154276,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=5,length=180915260,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=6,length=171115067,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=7,length=159138663,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=8,length=146364022,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=9,length=141213431,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=10,length=135534747,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=11,length=135006516,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=12,length=133851895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=13,length=115169878,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=14,length=107349540,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=15,length=102531392,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=16,length=90354753,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=17,length=81195210,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=18,length=78077248,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=19,length=59128983,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=20,length=63025520,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=21,length=48129895,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=22,length=51304566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=X,length=155270560,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=Y,length=59373566,assembly=37,species=HUMAN>\n");
	strcat(exp,	"##contig=<ID=MT,length=16569,assembly=37,species=HUMAN>\n");

	strcat(exp,output_generate_info_lines());
	strcat(exp,output_generate_format_lines());

	strcat(exp,"##SAMPLE=<ID=NORMAL,Description=\"Normal\",Accession=.,Platform=TEST,Protocol=WGS,SampleName=NORMALb,Source=.>\n");
	strcat(exp,"##SAMPLE=<ID=TUMOUR,Description=\"Tumour\",Accession=.,Platform=TEST,Protocol=WXS,SampleName=TUMOURa,Source=.>\n");
	strcat(exp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOUR\n");

	count = 0;
	exp_lines = 52;
	strcpy(got,"");
	while ( fgets(line,sizeof(line),out) != NULL ){
		mu_assert(count<=exp_lines,"Too many header lines.");
		strcat(got,line);
		count++;
	}
	if(strcmp(exp,got)!=0){
		printf("exp:\n%s\n\n\ngot:\n%s\n\n",exp,got);
	}
	mu_assert(strcmp(exp,got)==0,"Header in file doesn't match.");
	fclose(out);
	unlink(out_test_vcf);

	return NULL;
}


char *all_tests() {
   mu_suite_start();
   mu_run_test(test_output_generate_info_lines);
   mu_run_test(test_output_generate_format_lines);
   mu_run_test(test_output_to_file);
   mu_run_test(test_output_header_to_file);
   mu_run_test(test_output_header_to_file_cram);
	//int output_vcf_header(FILE *out, char *tum_bam, char *norm_bam, char *ref_seq_loc);
   return NULL;
}

RUN_TESTS(all_tests);
