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

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <List.h>
#include <dbg.h>
#include <file_tests.h>
#include <covs_access.h>
#include <cn_access.h>
#include <alg_bean.h>
#include <bam_access.h>
#include <split_access.h>
#include <ignore_reg_access.h>
#include <output.h>
#include <fai_access.h>
#include <algos.h>
#include <config_file_access.h>

static char tum_bam_file[512];// = NULL;
static char norm_bam_file[512];// = NULL;
static char ignore_regions_file[512];// = NULL;
static char *config_file = "caveman.cfg.ini";
static char results[512];// = "results";
static char ref_idx[512];// = "";
static char list_loc[512];// = "splitList";
static char alg_bean_loc[512];// = "alg_bean";
static char version[50];
static char *covariate_file = "covs_arr";
static char *probs_file = "probs_arr";
static char norm_cn_loc[512];
static char tum_cn_loc[512];
static char *norm_plat = NULL;
static char *tum_plat = NULL;
static float min_mut_prob = 0.8;
static float min_snp_prob = 0.95;
static float norm_contam = 0.1;
static float ref_bias = 0.95;
static float prior_mut_prob = 0.000006;
static float prior_snp_prob = 0.0001;
static int min_tum_cvg = 1;
static int min_norm_cvg = 1;
static int estep_max_tumour_coverage = 25000;
static int normal_copy_number = 2;
static int tumour_copy_number = 2;
static int includeSW = 0;
static int includeSingleEnd = 0;
static int includeDups = 0;
static int min_bq = 11;
static int cn = 0;
static int idx;
static int split_size = 50000;
static int debug=0;
static char *assembly = NULL;
static char *species = NULL;
static char *norm_prot = "WGS";
static char *tum_prot = "WGS";
static int max_copy_number = 10;
char *valid_protocols[3] = {"WGS","WXS","RNA"};

void estep_print_usage (int exit_code){
	printf ("Usage: caveman estep -i jobindex [-f file] [-m int] [-k float] [-b float] [-p float] [-q float] [-x int] [-y int] [-c float] [-d float] [-a int]\n\n");
	printf("-i  --index [int]                                Job index (e.g. from $LSB_JOBINDEX)\n\n");
	printf("Optional\n");
	printf("-f  --config-file [file]                         Path to the config file produced by setup. [default:'%s']\n",config_file);
	printf("-m  --min-base-qual [int]                        Minimum base quality for inclusion of a read position [default:%d]\n",min_bq);
	printf("-c  --prior-mut-probability [float]              Prior somatic probability [default:%f]\n",prior_mut_prob);
	printf("-d  --prior-snp-probability [float]              Prior germline mutant probability [default:%f]\n",prior_snp_prob);
	printf("-k  --normal-contamination [float]               Normal contamination of tumour [default:%f]\n",norm_contam);
	printf("-b  --reference-bias [float]                     Reference bias [default:%f]\n",ref_bias);
	printf("-p  --mut-probability-cutoff [float]             Minimum probability call for a somatic mutant position to be output [default:%f]\n",min_mut_prob);
	printf("-q  --snp-probability-cutoff [float]             Minimum probability call for a germline mutant position to be output [default:%f]\n",min_snp_prob);
	printf("-x  --min-tum-coverage [int]                     Minimum tumour coverage for analysis of a position [default:%d]\n",min_tum_cvg);
	printf("-y  --min-norm-coverage [int]                    Minimum normal coverage for analysis of a position [default:%d]\n",min_norm_cvg);
	printf("-a  --split-size [int]                           Size of section to retrieve at a time from bam file. Allows memory footprint tuning [default:%d].\n",split_size);
	printf("-s  --debug                                      Adds an extra output to a debug file. Every base analysed has an output\n");
	printf("-g  --cov-file [file]                            File location of the covariate array. [default:'%s']\n",covariate_file);
	printf("-o  --prob-file [file]                           File location of the prob array. [default:'%s']\n",probs_file);
	printf("-v  --species-assembly [string]                  Species assembly (eg 37/GRCh37), required if bam header SQ lines do not contain AS and SP information.\n");
	printf("-w  --species [string]                           Species name (eg Human), required if bam header SQ lines do not contain AS and SP information.\n");
	printf("-n  --normal-copy-number [int]                   Copy number to use when filling gaps in the normal copy number file [default:%d].\n",normal_copy_number);
	printf("-t  --tumour-copy-number [int]                   Copy number to use when filling gaps in the tumour copy number file [default:%d].\n",tumour_copy_number);
	printf("-l  --normal-protocol [string]                   Normal protocol. Ideally this should match -r but not checked (WGS|WXS|RNA) [default:%s].\n",norm_prot);
	printf("-r  --tumour-protocol [string]                   Tumour protocol. Ideally this should match -l but not checked (WGS|WXS|RNA) [default:%s].\n",tum_prot);
	printf("-P  --normal-platform [string]                   Normal platform. Overrides the values retrieved from bam header.\n");
	printf("-T  --tumour-platform [string]                   Tumour platform. Overrides the values retrieved from bam header.\n");
	printf("-M  --max-copy-number [int]                      Maximum copy number permitted. If exceeded the copy number for the offending region will be set to this value. [default:%d].\n",max_copy_number);
  printf("-h	help                                         Display this usage information.\n");

  exit(exit_code);
}

void estep_setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
             	{"config-file", required_argument, 0, 'f'},
             	{"cov-file",required_argument,0,'g'},
             	{"prob-file",required_argument,0,'o'},
             	{"min-base-qual",required_argument , 0, 'm'},
             	{"index", required_argument, 0, 'i'},
             	{"normal-contamination", required_argument, 0, 'k'},
             	{"prior-snp-probability", required_argument, 0, 'd'},
             	{"prior-mut-probability", required_argument, 0, 'c'},
             	{"reference-bias", required_argument, 0, 'b'},
             	{"mut-probability-cutoff", required_argument, 0, 'p'},
             	{"snp-probability-cutoff", required_argument, 0, 'q'},
             	{"min-tum-coverage", required_argument, 0, 'x'},
             	{"min-norm-coverage", required_argument, 0, 'y'},
             	{"split-size", required_argument, 0, 'a'},
             	{"species-assembly ", required_argument, 0, 'v'},
             	{"species", required_argument, 0, 'w'},
             	{"normal-copy-number", required_argument, 0, 'n'},
             	{"tumour-copy-number", required_argument, 0, 't'},
             	{"normal-protocol", required_argument, 0, 'l'},
             	{"tumour-protocol", required_argument, 0, 'r'},
             	{"max-copy-number", required_argument, 0, 'M'},
             	{"normal-platform", required_argument, 0, 'P'},
             	{"tumour-platform", required_argument, 0, 'T'},
             	{"snp-warnings", no_argument, 0, 'S'},
             	{"help", no_argument, 0, 'h'},
             	{"debug", no_argument, 0, 's'},

             	{ NULL, 0, NULL, 0}
   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "x:y:c:d:p:q:b:k:a:f:i:o:g:m:n:t:v:w:l:M:P:T:r:sh",
                            								long_opts, &index)) != -1){
   	switch(iarg){
   		case 'l':
				norm_prot = optarg;
   			break;

   		case 'r':
				tum_prot = optarg;
   			break;

   		case 'v':
   			assembly = optarg;
   			break;

   		case 'w':
   			species = optarg;
   			break;

   		case 'o':
   			probs_file = optarg;
   			break;

   		case 'g':
   			covariate_file = optarg;
   			break;

   		case 's':
   			debug = 1;
   			break;

   		case 'h':
        estep_print_usage(0);
        break;

      case 'f':
        config_file = optarg;
        break;

      case 'M':
      	if(sscanf(optarg, "%i", &cn) != 1){
      		sentinel("Error parsing -M argument '%s'. Should be an integer > 0",optarg);
      	}
        cn_access_set_max_cn(cn);
        break;

      case 'P':
				norm_plat = optarg;
        break;

      case 'T':
      	tum_plat = optarg;
        break;

      case 'n':
      	if(sscanf(optarg, "%i", &normal_copy_number) != 1){
      		sentinel("Error parsing -n argument '%s'. Should be an integer > 0",optarg);
      	}
        break;

      case 't':
      	if(sscanf(optarg, "%i", &tumour_copy_number) != 1){
      		sentinel("Error parsing -t argument '%s'. Should be an integer > 0",optarg);
      	}
        break;

      case 'i':
      	if(sscanf(optarg, "%i", &idx) != 1){
      		sentinel("Error parsing -i argument '%s'. Should be an integer > 0",optarg);
      	}
        break;

      case 'm':
      	if(sscanf(optarg, "%i", &min_bq) != 1){
      		sentinel("Error parsing -m argument '%s'. Should be an integer >= 0",optarg);
      	}
        break;

			case 'k':
				if(sscanf(optarg, "%f", &norm_contam) != 1){
      		sentinel("Error parsing -k argument '%s'. Should be a float >= 0.0.",optarg);
      	}
				break;

			case 'd':
				if(sscanf(optarg, "%f", &prior_snp_prob) != 1){
      		sentinel("Error parsing -d argument '%s'. Should be a float > 0.0.",optarg);
      	}
				break;

			case 'c':
				if(sscanf(optarg, "%f", &prior_mut_prob) != 1){
      		sentinel("Error parsing -c argument '%s'. Should be a float > 0.0.",optarg);
      	}
				break;

			case 'b':
				if(sscanf(optarg, "%f", &ref_bias) != 1){
      		sentinel("Error parsing -b argument '%s'. Should be a float >= 0.0.",optarg);
      	}
				break;

			case 'p':
				if(sscanf(optarg, "%f", &min_mut_prob) != 1){
      		sentinel("Error parsing -p argument '%s'. Should be a float >= 0.0.",optarg);
      	}
				break;

			case 'q':
				if(sscanf(optarg, "%f", &min_snp_prob) != 1){
      		sentinel("Error parsing -q argument '%s'. Should be a float >= 0.0.",optarg);
      	}
				break;

			case 'x':
				if(sscanf(optarg, "%i", &min_tum_cvg) != 1){
      		sentinel("Error parsing -x argument '%s'. Should be an integer > 0",optarg);
      	}
				break;

			case 'y':
				if(sscanf(optarg, "%i", &min_norm_cvg) != 1){
      		sentinel("Error parsing -y argument '%s'. Should be an integer > 0",optarg);
      	}
				break;

			case 'a':
				if(sscanf(optarg, "%i", &split_size) != 1){
      		sentinel("Error parsing -a argument '%s'. Should be an integer >= 0",optarg);
      	}
				break;

			case 'S':
				set_snp_warnings();
				break;

			case '?':
        estep_print_usage (1);
        break;

      default:
      	estep_print_usage (1);

   	}; // End of args switch statement

   }//End of iteration through options

   //Do some checking to ensure required arguments were passed
   if(idx == 0){
   	estep_print_usage(1);
   }

   //Do some checking to ensure required arguments were passed and are accessible files
   if(check_exist(config_file) != 1){
   	printf("Config file %s does not appear to exist. Have you run previous caveman steps?\n",config_file);
   	estep_print_usage(1);
   }

	 int i=0;
	 int norm_prot_check=0;
	 int tum_prot_check=0;
	 for(i=0;i<3;i++){
		if(strcmp(norm_prot,valid_protocols[i])==0){
			norm_prot_check=1;
		}
		if(strcmp(tum_prot,valid_protocols[i])==0){
			tum_prot_check=1;
		}
	 }

   if(norm_prot_check==0){
		printf("Normal protocol '%s' is invalid should be one of (WGS|WXS|RNA).",norm_prot);
		estep_print_usage(1);
   }

   if(tum_prot_check==0){
		printf("Tumour protocol '%s' is invalid should be one of (WGS|WXS|RNA).",tum_prot);
		estep_print_usage(1);
   }

   set_max_tum_cvg(estep_max_tumour_coverage);

   return;

error:
	estep_print_usage (1);
	return;
}

int estep_main(int argc, char *argv[]){
	estep_setup_options(argc, argv);

	long double ********prob_arr = NULL;
	alg_bean_t *alg = NULL;
	char *fa_file = NULL;
	int ignore_reg_count = 0;
	List *these_regions = NULL;
	struct seq_region_t **ignore_regs = NULL;
	char no_analysis_file_loc[500];
	FILE *no_analysis_file = NULL;
	List *no_analysis_list = NULL;
	char *ref_seq = NULL;
	gzFile debug_file = NULL;
	gzFile snp_file = NULL;
	gzFile mut_file = NULL;

	//Open the config file and do relevant things
	FILE *config = fopen(config_file,"r");
	check(config != NULL,"Failed to open config file for reading. Have you run caveman-setup?");
	int cfg = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,
														ref_idx,ignore_regions_file,alg_bean_loc,results,list_loc,
																					&includeSW,&includeSingleEnd,&includeDups,version,norm_cn_loc,tum_cn_loc);

	check(strcmp(version,CAVEMAN_VERSION)==0,"Stored version in %s %s and current code version %s did not match.",config_file,version,CAVEMAN_VERSION);

	check(cfg==0,"Error parsing config file.");
  bam_access_include_sw(includeSW);
  bam_access_include_se(includeSingleEnd);
  bam_access_include_dup(includeDups);

  set_normal_cn(normal_copy_number);
  set_tumour_cn(tumour_copy_number);

  if(norm_plat==NULL){
  	norm_plat = malloc(sizeof(char) * 50);
		check_mem(norm_plat);
		strcpy(norm_plat,".");
  }
  if(tum_plat==NULL){
		tum_plat = malloc(sizeof(char) * 50);
		check_mem(tum_plat);
		strcpy(tum_plat,".");
  }

	//Load in alg bean
	FILE *alg_bean_file = fopen(alg_bean_loc,"r");
	check(alg_bean_file != 0 ,"Error trying to open alg_bean file: %s.",alg_bean_loc);
	alg = alg_bean_read_file(alg_bean_file);
	check(alg != NULL,"Error reading alg_bean from file.");
	check(fclose(alg_bean_file)==0,"Error closing alg_bean_file.");

	//Load in probability array
	prob_arr = covs_access_read_probs_from_file(probs_file,
														List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
														List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),
																												List_count(alg->ref_base),List_count(alg->call_base));

	//Set the algorithm modifiers and open the bam files
	//Set the min base qual in case it's been changed.
	bam_access_min_base_qual(min_bq);

	//Open the bam files.
	bam_access_openbams(norm_bam_file,tum_bam_file,ref_idx);

	//Estep algorithm modifiers.
	set_min_mut_prob(min_mut_prob);
	set_min_snp_prob(min_snp_prob);
	set_norm_contam(norm_contam);
	set_ref_bias(ref_bias);
	set_prior_mut_prob(prior_mut_prob);
 	set_prior_snp_prob(prior_snp_prob);
 	set_min_tum_cvg(min_tum_cvg);
 	set_min_norm_cvg(min_norm_cvg);

	//Get split section given index
	char chr_name[50];
	int start_zero_based = 0;
	int stop = 0;
	split_access_get_section_from_index(list_loc,chr_name,&start_zero_based,&stop,idx);
	check(stop > 0,"Error fetching region from split file.");
	check(chr_name != NULL, "Error fetching region from split file.");

	//Get the chunk of ref sequence using this split section.
	//Strip the .fai from the fasta file.
	fa_file = malloc(sizeof(char) * (strlen(ref_idx)-3));
	check_mem(fa_file);
	strncpy(fa_file,ref_idx,(strlen(ref_idx)-4));
	fa_file[strlen(ref_idx)-4] = '\0';
	check(fa_file != NULL, "Error decoding FASTA file name");

  //Open no analsyis file here so we can write ignored regions.
	int chk_no = sprintf(no_analysis_file_loc,"%s/%s/%d_%d.no_analysis.bed",results,chr_name,start_zero_based+1,stop);
	check(chk_no>0,"Error generating no analysis file location.");
	no_analysis_file = fopen(no_analysis_file_loc,"w");
	check(no_analysis_file != 0, "Error trying to open no analysis file for output: %s.",no_analysis_file_loc);
	output_set_no_analysis_file(no_analysis_file);

	//Get ignored regions contained in split section
	//Get ignored regions for section and calculate a list of sections to analyse.
	ignore_reg_count = ignore_reg_access_get_ign_reg_count_for_chr(ignore_regions_file,chr_name);
   check(ignore_reg_count >= 0,"Error trying to check the number of ignored regions for this chromosome.");

   //A list structure to store the ignored regions and skipped bases.
   no_analysis_list = List_create();
   output_set_no_analysis_section_list(no_analysis_list);

   //Now create a store for said regions.
   ignore_regs = malloc(sizeof(struct seq_region_t *) *  ignore_reg_count);
   check_mem(ignore_regs);
   check(ignore_reg_access_get_ign_reg_for_chr(ignore_regions_file,chr_name,ignore_reg_count,ignore_regs)==0,"Error fetching ignored regions from file.");

	//Create a list of sections to analyse.
	//Check the contained ignored regions
	//Resolve the ignored regions and start/stop into sections for analysis.
	these_regions = ignore_reg_access_resolve_ignores_to_analysis_sections(start_zero_based+1,stop,ignore_regs,ignore_reg_count);

	//Add each of these ignored regions to the no analysis list
	int i=0;
	for(i=0; i<ignore_reg_count; i++){
		output_append_position_to_no_analysis(chr_name, ignore_regs[i]->beg, ignore_regs[i]->end);
	}

	// If there's only one split it in two.
	if(List_count(these_regions) == 1){
		seq_region_t *old = List_pop(these_regions);
		if(old->beg == old->end){
			seq_region_t *range_1 = malloc(sizeof(seq_region_t));
			range_1->beg = old->beg;
			range_1->end = old->end;
			List_push(these_regions,range_1);
			free(old);
		}else{
			seq_region_t *range_1 = malloc(sizeof(seq_region_t));
			seq_region_t *range_2 = malloc(sizeof(seq_region_t));
			range_1->beg = old->beg;
			int halfway = ((old->end - old->beg)/2);
			range_1->end = old->beg+halfway;
			range_2->beg = old->beg+halfway+1;
			range_2->end = old->end;
			free(old);
			List_push(these_regions,range_1);
			List_push(these_regions,range_2);
		}
	}

	//Generate snp and mut out filenames for this section,
	char snp_out[500];
	char mut_out[500];
	char debug_out[500];

	int chk = sprintf(snp_out,"%s/%s/%d_%d.snps.vcf.gz",results,chr_name,start_zero_based+1,stop);
	check(chk>0,"Error generating snp output file location.");

	chk = sprintf(mut_out,"%s/%s/%d_%d.muts.vcf.gz",results,chr_name,start_zero_based+1,stop);
	check(chk>0,"Error generating mut output file location.");

	chk = sprintf(debug_out,"%s/%s/%d_%d.dbg.vcf.gz",results,chr_name,start_zero_based+1,stop);
	check(chk>0,"Error generating debug file location.");

	//Open files for output
	mut_file = gzopen(mut_out,"wb1");
	check(mut_file != 0, "Error trying to open mut file for output: %s.",mut_out);
	int chk_write = output_vcf_header(mut_file, tum_bam_file, norm_bam_file, fa_file,
																									assembly, species, norm_prot, tum_prot,
																									norm_plat, tum_plat);
	check(chk_write==0,"Error writing header to muts file.");

	snp_file = gzopen(snp_out,"wb1");
	check(snp_file != 0, "Error trying to open snp file for output: %s.",snp_out);
	chk_write = output_vcf_header(snp_file, tum_bam_file, norm_bam_file, fa_file,
																									assembly, species, norm_prot, tum_prot,
																									norm_plat, tum_plat);
	check(chk_write==0,"Error writing header to SNP file.");


	if(debug == 1){
		debug_file = gzopen(debug_out,"wb1");
		check(debug_file != 0, "Error trying to open snp file for output: %s.",debug_out);
		chk_write = output_vcf_header(debug_file, tum_bam_file, norm_bam_file, fa_file,
																									assembly, species, norm_prot, tum_prot,
																									norm_plat, tum_plat);
		check(chk_write==0,"Error writing header to dbg file.");
	}


	//Iterate through analysis sections
	//Iterate through sections.
	LIST_FOREACH(these_regions, first, next, cur){
		printf("Estep section %s:%d-%d\n",chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);
		//Get the reference sequence for this section
		ref_seq = fai_access_get_ref_seqeuence_for_pos(fa_file,chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);

		printf("fetched a reference seq of length %lu for this section.\n",strlen(ref_seq));
		check(ref_seq != NULL,"Error retrieving reference sequence for section %s:%d-%d.",chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);
		//Get all reads or pos pileups for section.
		//Iterate through positions in section and mstep
		int chk = algos_estep_read_position(alg,prob_arr,chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end,ref_seq,
									norm_cn_loc, tum_cn_loc, snp_file, mut_file, debug_file, split_size);
		check(chk==0,"Error running estep for region %s:%d-%d.",chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);
		free(ref_seq);
	}
	int write_check = output_flush_no_analysis(chr_name);
	check(write_check==0,"Error writing no analysis regions to bed file.");

	//Flush and Close output files.
	check(gzclose(snp_file)==0,"Error closing snp file '%s'.",snp_out);
	check(gzclose(mut_file)==0,"Error closing mut file '%s'.",mut_out);
	if(debug_file){
	 	check(gzclose(debug_file)==0,"Error closing debug file.");
	}
	int fcheck = fflush(no_analysis_file);
	check(fcheck==0,"Error flushing no_analysis_file.");
	check(fclose(no_analysis_file)==0,"Error closing no analysis file.");
	//cleanup
	List_clear_destroy(these_regions);
	free(fa_file);
	bam_access_closebams();
	ignore_reg_access_destroy_seq_region_t_arr(ignore_reg_count,ignore_regs);
	covs_access_free_prob_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
						List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),prob_arr);
	alg_bean_destroy(alg);
	return 0;
error:
	if(these_regions) List_clear_destroy(these_regions);
	if(no_analysis_list) List_clear_destroy(no_analysis_list);
	if(fa_file) free(fa_file);
	if(ref_seq) free(ref_seq);
	if(mut_file) gzclose(mut_file);
	if(snp_file) gzclose(snp_file);
	if(debug_file) gzclose(debug_file);
	if(no_analysis_file) fclose(no_analysis_file);
	bam_access_closebams();
	ignore_reg_access_destroy_seq_region_t_arr(ignore_reg_count,ignore_regs);
	if(prob_arr) covs_access_free_prob_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
						List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),prob_arr);
	if(alg) alg_bean_destroy(alg);
	return -1;
}
