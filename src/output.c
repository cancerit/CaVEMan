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

#include <stdio.h>
#include <assert.h>
#include <output.h>
#include <string.h>
#include <bam_access.h>
#include <ignore_reg_access.h>
#include <genotype.h>
#include <stdio.h>
#include <time.h>
#include <dbg.h>

static const char *VCF_VERSION_KEY = "fileformat";
static const char *VCF_VERSION_VALUE = "VCFv4.1";
static const char *VCF_FILEDATE_KEY = "fileDate";
static const char *VCF_REF_SEQ_KEY = "reference";
static const char *VCF_SAMPLE_KEY = "SAMPLE";
static const char *VCF_NORMAL_NAME = "NORMAL";
static const char *VCF_DESCRIPTION_KEY = "Description";
static const char *VCF_ACCESSION_KEY = "Accession";
static const char *VCF_PLATFORM_KEY = "Platform";
static const char *VCF_PROTOCOL_KEY = "Protocol";
static const char *VCF_INDIVIDUAL_KEY = "SampleName";
static const char *VCF_SOURCE_KEY = "Source";
static const char *VCF_TUMOUR_NAME = "TUMOUR";
static const char *VCF_PROCESSLOG_KEY = "vcfProcessLog";
static const char *VCF_INPUTVCF_KEY = "InputVCF";
static const char *VCF_INPUTSOURCE_KEY = "InputVCFSource";
static const char *VCF_INPUT_SOURCE_VALUE = "CaVEMan";
static const char *VCF_INPUTVERSION_KEY = "InpuVCFVer";
static const char *VCF_INPUTPARAMS_KEY = "InputVCFParam";
static const char *VCF_FORMAT_KEY = "FORMAT";
static const char *VCF_INFO_KEY = "INFO";
static const char *VCF_VAR_FORMAT_LINE = "GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM";
static const char *VCF_FORMAT_A_FWD_COUNT = "ID=FAZ,Number=1,Type=Integer,Description=\"Reads presenting a A for this position, forward strand\"";
static const char *VCF_FORMAT_C_FWD_COUNT = "ID=FCZ,Number=1,Type=Integer,Description=\"Reads presenting a C for this position, forward strand\"";
static const char *VCF_FORMAT_G_FWD_COUNT = "ID=FGZ,Number=1,Type=Integer,Description=\"Reads presenting a G for this position, forward strand\"";
static const char *VCF_FORMAT_T_FWD_COUNT = "ID=FTZ,Number=1,Type=Integer,Description=\"Reads presenting a T for this position, forward strand\"";
static const char *VCF_FORMAT_A_REV_COUNT = "ID=RAZ,Number=1,Type=Integer,Description=\"Reads presenting a A for this position, reverse strand\"";
static const char *VCF_FORMAT_C_REV_COUNT = "ID=RCZ,Number=1,Type=Integer,Description=\"Reads presenting a C for this position, reverse strand\"";
static const char *VCF_FORMAT_G_REV_COUNT = "ID=RGZ,Number=1,Type=Integer,Description=\"Reads presenting a G for this position, reverse strand\"";
static const char *VCF_FORMAT_T_REV_COUNT = "ID=RTZ,Number=1,Type=Integer,Description=\"Reads presenting a T for this position, reverse strand\"";
static const char *VCF_FORMAT_MUT_ALL_PROP = "ID=PM,Number=1,Type=Float,Description=\"Proportion of mut allele\"";
static const char *VCF_FORMAT_GENOTYPE = "ID=GT,Number=1,Type=String,Description=\"Genotype\"";
static const char *VCF_INFO_TOTAL_DEPTH = "ID=DP,Number=1,Type=Integer,Description=\"Total Depth\"";
static const char *VCF_INFO_MUT_PROB_SUM = "ID=MP,Number=1,Type=Float,Description=\"Sum of CaVEMan somatic genotype probabilities\"";
static const char *VCF_INFO_SNP_PROB_SUM = "ID=GP,Number=1,Type=Float,Description=\"Sum of CaVEMan germline genotype probabilities\"";
static const char *VCF_INFO_VAR_ALL_PROB_SUM = "ID=SP,Number=1,Type=Float,Description=\"Sum of CaVEMan somatic genotype containing the called mutant allele probabilities\"";
static const char *VCF_INFO_TOP_GENO = "ID=TG,Number=1,Type=String,Description=\"Most probable genotype as called by CaVEMan\"";
static const char *VCF_INFO_TOP_GENO_PROB = "ID=TP,Number=1,Type=Float,Description=\"Probability of most probable genotype as called by CaVEMan\"";
static const char *VCF_INFO_SEC_GENO = "ID=SG,Number=1,Type=String,Description=\"2nd most probable genotype as called by CaVEMan\"";
static const char *VCF_INFO_SEC_GENO_PROB = "ID=SP,Number=1,Type=Float,Description=\"Probability of 2nd most probable genotype as called by CaVEMan\"";
static const char *VCF_INFO_DBNSP_ID = "ID=DS,Number=.,Type=String,Description=\"DBSnp ID of known SNP\"";

FILE *no_analysis = NULL;
int no_analysis_cache_size = 500;
List *no_analysis_sects = NULL;

char algos_bases[4] = {'A','C','G','T'};

void output_set_genotype_representations_for_genotype_t(genotype_t *geno,char ref_base, int *norm, int *tum){
	if(genotype_get_base_count(geno,ref_base) == 0){
		*norm = 1;
		*tum = 1;
		return;
	}else{
		*norm = 0;
		*tum = 0;
		int i=0;
		for(i=0;i<4;i++){
			if(algos_bases[i]!=ref_base && genotype_get_base_count(geno, algos_bases[i]) > 0){
				*tum = 1;
				return;
			}
		}
		return;
	}
}

int output_vcf_variant_position(estep_position_t *pos, gzFile out, char *chrom){
	assert(pos != NULL);
	assert(out != NULL);
	assert(chrom != NULL);

	char *top_ref = NULL;
	char *top_tum = NULL;
	char *sec_ref = NULL;
	char *sec_tum = NULL;

	int write = gzprintf(out,"%s\t%d\t.\t%s\t%c\t.\t.\t",chrom,pos->ref_pos,pos->ref_base,pos->top_geno->tum_geno->var_base);
	check(write>0,"Error writing initial vcf variant line, top mutant genotype.");
	//total depth seen by CaVEMan
	write = gzprintf(out,"DP=%d;",(genotype_get_total_base_count(pos->norm_fwd_cvg)+genotype_get_total_base_count(pos->norm_rev_cvg)
								+genotype_get_total_base_count(pos->tum_fwd_cvg)+genotype_get_total_base_count(pos->tum_rev_cvg)));
	check(write>0,"Error writing info line.");
	//Total somatic prob
	write = gzprintf(out,"MP=%5.1Le;",pos->total_mut_prob);
	check(write>0,"Error writing info line, mut_prob.");
	//total germline prob
	write = gzprintf(out,"GP=%5.1Le;",pos->total_snp_prob);
	check(write>0,"Error writing info line snp_prob.");
    //Total somatic probability for variant allele changes only
    write = gzprintf(out,"SP=%5.1Le;",pos->total_mut_allele_prob);
	check(write>0,"Error writing info line var allele total prob.");
	//Top genotype and second genotype info
	top_ref = genotype_get_genotype_t_as_string(pos->top_geno->norm_geno);
	top_tum = genotype_get_genotype_t_as_string(pos->top_geno->tum_geno);
	sec_ref = genotype_get_genotype_t_as_string(pos->sec_geno->norm_geno);
	sec_tum = genotype_get_genotype_t_as_string(pos->sec_geno->tum_geno);
	write = gzprintf(out,"TG=%s/%s;TP=%5.1Le;SG=%s/%s;SP=%5.1Le\t",
						top_ref,top_tum,pos->top_geno->prob,
						sec_ref,sec_tum,pos->sec_geno->prob);

	check(write>0,"Error writing info line extra probs.");
	//FORMAT line
	write = gzprintf(out,"%s\t",VCF_VAR_FORMAT_LINE);
	check(write>0,"Error writing vcf variant format line.");
	//NORMAL FORMAT
	int norm_geno_rep = 0;
	int tum_geno_rep = 0;
	int mut_count = genotype_get_base_count(pos->norm_fwd_cvg,pos->top_geno->tum_geno->var_base)+genotype_get_base_count(pos->norm_rev_cvg,pos->top_geno->tum_geno->var_base);
	double mut_prop = 0.0;
	if(mut_count > 0){
		mut_prop = (double) mut_count / (double)(genotype_get_total_base_count(pos->norm_fwd_cvg)+genotype_get_total_base_count(pos->norm_rev_cvg));
	}
	output_set_genotype_representations_for_genotype_t(pos->top_geno->norm_geno,pos->ref_base[0], &norm_geno_rep, &tum_geno_rep);
	write = gzprintf(out,"%d|%d:%d:%d:%d:%d:%d:%d:%d:%d:%5.1e\t",
						norm_geno_rep,tum_geno_rep,
						pos->norm_fwd_cvg->a_count,pos->norm_fwd_cvg->c_count,pos->norm_fwd_cvg->g_count,pos->norm_fwd_cvg->t_count,
						pos->norm_rev_cvg->a_count,pos->norm_rev_cvg->c_count,pos->norm_rev_cvg->g_count,pos->norm_rev_cvg->t_count,
						mut_prop);
	check(write>0,"Error writing normal format for VCF variant line.");
	//TUMOUR FORMAT
	norm_geno_rep = 0;
	tum_geno_rep = 0;
	mut_count = 0;
	mut_count = genotype_get_base_count(pos->tum_fwd_cvg,pos->top_geno->tum_geno->var_base)+genotype_get_base_count(pos->tum_rev_cvg,pos->top_geno->tum_geno->var_base);
	mut_prop = 0.0;
	if(mut_count > 0){
		mut_prop = (double) mut_count / (double)(genotype_get_total_base_count(pos->tum_fwd_cvg)+genotype_get_total_base_count(pos->tum_rev_cvg));
	}
	output_set_genotype_representations_for_genotype_t(pos->top_geno->tum_geno,pos->ref_base[0], &norm_geno_rep, &tum_geno_rep);
	write = gzprintf(out,"%d|%d:%d:%d:%d:%d:%d:%d:%d:%d:%5.1e\n",
						norm_geno_rep,tum_geno_rep,
						pos->tum_fwd_cvg->a_count,pos->tum_fwd_cvg->c_count,pos->tum_fwd_cvg->g_count,pos->tum_fwd_cvg->t_count,
						pos->tum_rev_cvg->a_count,pos->tum_rev_cvg->c_count,pos->tum_rev_cvg->g_count,pos->tum_rev_cvg->t_count,
						mut_prop);
	check(write>0,"Error writing tumour format for VCF variant line.");
	free(top_ref);
	free(top_tum);
	free(sec_ref);
	free(sec_tum);
	return 0;
error:
	if(top_ref) free(top_ref);
	if(top_tum) free(top_tum);
	if(sec_ref) free(sec_ref);
	if(sec_tum) free(sec_tum);
	return -1;
}

char *output_generate_CaVEMan_process_log(char *cave_version){
	char *process = malloc(sizeof(char) * 1000);
	sprintf(process,"##%s=<%s=<.>,%s=<%s>,%s=<\"%s\">,%s=<NORMAL_CONTAMINATION=%g,REF_BIAS=%g,PRIOR_MUT_RATE=%g,PRIOR_SNP_RATE=%g,SNP_CUTOFF=%g,MUT_CUTOFF=%g>>\n##cavemanVersion=%s\n",
						VCF_PROCESSLOG_KEY,VCF_INPUTVCF_KEY,VCF_INPUTSOURCE_KEY,VCF_INPUT_SOURCE_VALUE,
						VCF_INPUTVERSION_KEY,cave_version,
						VCF_INPUTPARAMS_KEY,get_norm_contam(),get_ref_bias(),get_prior_mut_prob(),get_prior_snp_prob(),get_min_snp_prob(),get_min_mut_prob(),
						cave_version);
	return process;
}

char *output_generate_reference_contig_lines(char *bam_file, char *assembly, char *species){
	assert(bam_file != NULL);
	List *contig_list = bam_access_get_contigs_from_bam(bam_file, assembly, species);
	check(contig_list != NULL,"Error fetching contigs from bam file.");
	char *contigs = malloc(sizeof(char) * 1000 * List_count(contig_list));
	strcpy(contigs,"");
	LIST_FOREACH(contig_list, first,next,cur){
		ref_seq_t *ref = (ref_seq_t *)cur->value;
		char contig_str[1000];
		sprintf(contig_str,"##contig=<ID=%s,length=%d,assembly=%s,species=%s>\n",ref->name,ref->length,ref->ass,ref->spp);
		strcat(contigs,contig_str);
	}
	List_clear_destroy(contig_list);
	return contigs;
error:
	if(contig_list) List_clear_destroy(contig_list);
	return NULL;
}


char *output_generate_info_lines(){
	char *info = malloc(sizeof(char) * 5000);
	sprintf(info,"##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n",
		VCF_INFO_KEY,VCF_INFO_TOTAL_DEPTH,
		VCF_INFO_KEY,VCF_INFO_MUT_PROB_SUM,
		VCF_INFO_KEY,VCF_INFO_SNP_PROB_SUM,
        VCF_INFO_KEY,VCF_INFO_VAR_ALL_PROB_SUM,
		VCF_INFO_KEY,VCF_INFO_TOP_GENO,
		VCF_INFO_KEY,VCF_INFO_TOP_GENO_PROB,
		VCF_INFO_KEY,VCF_INFO_SEC_GENO,
		VCF_INFO_KEY,VCF_INFO_SEC_GENO_PROB,
		VCF_INFO_KEY,VCF_INFO_DBNSP_ID);
	return info;
}

char *output_generate_format_lines(){
	char *format = malloc(sizeof(char) * 5000);
	sprintf(format,"##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n##%s=<%s>\n",
				VCF_FORMAT_KEY,VCF_FORMAT_GENOTYPE,
				VCF_FORMAT_KEY,VCF_FORMAT_A_FWD_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_C_FWD_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_G_FWD_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_T_FWD_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_A_REV_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_C_REV_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_G_REV_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_T_REV_COUNT,
				VCF_FORMAT_KEY,VCF_FORMAT_MUT_ALL_PROP);
	return format;
}

int output_vcf_header(gzFile out, char *tum_bam, char *norm_bam, char *ref_seq_loc,
													char *assembly, char *species, char *norm_prot, char *tum_prot,
													char *normal_platform, char *tumour_platform){

	assert(out != NULL);
	char *contigs = NULL;
	char *info = NULL;
	char *process = NULL;
	char *format = NULL;
	char *tumour_name = malloc(sizeof(char) * 150);
	char *normal_name = malloc(sizeof(char) * 150);
	check_mem(tumour_name);
	check_mem(normal_name);

	if(normal_platform == NULL){
		normal_platform = malloc(sizeof(char) * 150);
		check_mem(normal_platform);
	}

	char *res = bam_access_sample_name_platform_from_header(norm_bam,normal_name, normal_platform);
	check(res != NULL,"Error fetching normal sample and platform from bam header.");

	if(tumour_platform == NULL){
		tumour_platform = malloc(sizeof(char) * 150);
		check_mem(tumour_platform);
	}
	res=NULL;
	res = bam_access_sample_name_platform_from_header(tum_bam,tumour_name, tumour_platform);
	check(res != NULL,"Error fetching tumour sample and platform from bam header.");

	check(strcmp(tumour_platform,normal_platform)==0,"Normal and tumour platforms don't match: '%s' ne '%s'",normal_platform,tumour_platform);
	//VCF version (fileformat)
	int write = gzprintf(out,"##%s=%s\n",VCF_VERSION_KEY,VCF_VERSION_VALUE);
	check(write>0,"Error writing version.");
	//fileDate=20120104
	char date[50];
	time_t t = time(NULL);
	strftime(date,sizeof(date),"%Y%m%d",localtime(&t));
	write = gzprintf(out,"##%s=%s\n",VCF_FILEDATE_KEY,date);
	check(write>0,"Error writing filedate.");
	//Reference sequence
	write = gzprintf(out,"##%s=%s\n",VCF_REF_SEQ_KEY,ref_seq_loc);
	check(write>0,"Error writing reference.");
	//vcfProcessLog for CaVEMan
	char *cave_version = CAVEMAN_VERSION;
	process = output_generate_CaVEMan_process_log(cave_version);
	write = gzprintf(out,"%s",process);
	check(write>0,"Error writing CaVEMan version.");

	//Add reference sequence headers
	contigs = output_generate_reference_contig_lines(tum_bam, assembly, species);
	check(contigs != NULL,"Error fetching contigs from bam file.");
	write = gzputs(out,contigs);
	check(write==sizeof(char)*strlen(contigs),"Error (%d) writing contigs.",write);

	//INFO lines
	info = output_generate_info_lines();
	write = gzprintf(out,"%s",info);
	check(write==sizeof(char)*strlen(info),"Error (%d) writing INFO.",write);

	//FORMAT
	format = output_generate_format_lines();
	write = gzprintf(out,"%s",format);
	check(write==sizeof(char)*strlen(format),"Error (%d) writing FORMAT.",write);

	//SAMPLES
	//Normal
	write = gzprintf(out,"##%s=<ID=%s,%s=\"Normal\",%s=.,%s=%s,%s=%s,%s=%s,%s=.>\n",
			VCF_SAMPLE_KEY,VCF_NORMAL_NAME,VCF_DESCRIPTION_KEY,
			VCF_ACCESSION_KEY,VCF_PLATFORM_KEY,normal_platform,
			VCF_PROTOCOL_KEY,norm_prot,
			VCF_INDIVIDUAL_KEY,normal_name,
			VCF_SOURCE_KEY);
	check(write>0,"Error (%lu) writing normal sample.",write);
	//Tumour
	write = gzprintf(out,"##%s=<ID=%s,%s=\"Tumour\",%s=.,%s=%s,%s=%s,%s=%s,%s=.>\n",
			VCF_SAMPLE_KEY,VCF_TUMOUR_NAME,VCF_DESCRIPTION_KEY,
			VCF_ACCESSION_KEY,VCF_PLATFORM_KEY,tumour_platform,
			VCF_PROTOCOL_KEY,tum_prot,
			VCF_INDIVIDUAL_KEY,tumour_name,
			VCF_SOURCE_KEY);
	check(write>0,"Error (%d) writing tumour sample.",write);

	///Finally the line above output.
	write = gzprintf(out,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\t%s\n",VCF_NORMAL_NAME,VCF_TUMOUR_NAME);
	check(write>0,"Error writing column header line.");
	free(tumour_name);
	free(normal_name);
	free(contigs);
	free(info);
	free(format);
	free(process);
	return 0;
error:
	if(contigs) free(contigs);
	if(info) free(info);
	if(format) free(format);
	if(process) free(process);
	if(tumour_name) free(tumour_name);
	if(normal_name) free(normal_name);
	return -1;
}

int output_append_position_to_no_analysis(char *chr_name, int start_one_base, int stop){
	if(List_count(no_analysis_sects) > 0 && (((seq_region_t *)List_last(no_analysis_sects))->end + 1) == start_one_base){
		((seq_region_t *)List_last(no_analysis_sects))->end = stop;
	}else{
		struct seq_region_t *reg = malloc(sizeof(struct seq_region_t));
		check_mem(reg);
		reg->beg = start_one_base;
		reg->end = stop;
		List_push(no_analysis_sects,reg);
	}
	if(List_count(no_analysis_sects) > no_analysis_cache_size){
		int chk = output_flush_no_analysis(chr_name);
		check(chk==0,"Error writing to no analysis output file.");
		no_analysis_sects = List_create();
		check(no_analysis_sects!=NULL,"Error creating new list for no_analysis.");
	}
	return 0;
error:
	return -1;
}

int output_flush_no_analysis(char *chr_name){
  assert(no_analysis_sects!=NULL);
	//output all the no_analysis regions to file....
	LIST_FOREACH(no_analysis_sects, first, next, cur){
		seq_region_t *reg = (seq_region_t *) cur->value;
		int chk = fprintf(no_analysis,"%s\t%d\t%d\n",chr_name,(reg->beg-1),reg->end);
		check(chk>=0,"Error writing to no analysis bed file.");
	}
	List_clear_destroy(no_analysis_sects);
	return 0;
error:
	return -1;
}

void output_set_no_analysis_file(FILE *file){
	no_analysis = file;
}

void output_set_no_analysis_section_list(List *sections){
	no_analysis_sects = sections;
}
