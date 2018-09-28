/**   LICENSE
* Copyright (c) 2014-2018 Genome Research Ltd.
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
#include <dbg.h>
#include <mstep.h>
#include <alg_bean.h>
#include <fai_access.h>
#include <bam_access.h>
#include <file_tests.h>
#include <ignore_reg_access.h>
#include <split_access.h>
#include <covs_access.h>
#include <algos.h>
#include <config_file_access.h>
#include <sys/stat.h>

static int includeSW = 0;
static int includeSingleEnd = 0;
static int includeDups = 0;
static int min_bq = 11;
static char tum_bam_file[512];// = NULL;
static char norm_bam_file[512];// = NULL;
static char *config_file = "caveman.cfg.ini";
static char results[512];// = "results";
static char ref_idx[512];// = "";
static char list_loc[512];// = "splitList";
static char alg_bean[512];// = "alg_bean";
static char version[50];// = "alg_bean";
static char ignore_regions_file[512];// = NULL;
static char norm_cn_loc[512];
static char tum_cn_loc[512];
static int split_size = 50000;
static int idx = 0;

void mstep_print_usage (int exit_code){
	printf ("Usage: caveman mstep -i jobindex [-f file] [-m int] [-a int]\n\n");

	printf("-i  --index [int]                 Job index.\n\n");
	printf("Optional\n");
	printf("-f  --config-file [file]          Path to the config file produced by setup [default: '%s'].\n",config_file);
	printf("-a  --split_size [int]            Size of section to retrieve at a time from bam file. Allows memory footprint tuning [default:%d].\n",split_size);
	printf("-m  --min-base-qual [int]         Minimum base quality for inclusion of a read position in analysis [default:%d]\n",min_bq);
	printf("-h	help                          Display this usage information.\n");
  exit(exit_code);
}

int mstep_setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
		         {"config-file", required_argument, 0, 'f'},
             	{"min-base-qual",required_argument , 0, 'm'},
             	{"index", required_argument, 0, 'i'},
             	{"split_size", required_argument, 0, 'a'},
             	{"help", no_argument, 0, 'h'},
             	{ NULL, 0, NULL, 0}
   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "f:i:m:a:h",
                            								long_opts, &index)) != -1){
   	switch(iarg){
   		case 'h':
         	mstep_print_usage(0);
         	break;

      	case 'f':
      		config_file = optarg;
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

      	case 'a':
      		if(sscanf(optarg, "%i", &split_size) != 1){
      			sentinel("Error parsing -a argument '%s'. Should be an integer > 0",optarg);
      		}
					break;

			case '?':
            mstep_print_usage (1);
            break;

      	default:
      		mstep_print_usage (1);

   	}; // End of args switch statement

   }//End of iteration through options

   //Do some checking to ensure required arguments were passed
   if(idx == 0){
   	mstep_print_usage(1);
   }

   //Do some checking to ensure required arguments were passed and are accessible files
   if(check_exist(config_file) != 1){
   	printf("Config file %s does not appear to exist. Have you run previous caveman steps?\n",config_file);
   	mstep_print_usage(1);
   }

   return 0;

error:
	mstep_print_usage (1);
	return -1;
}

int mstep_main(int argc, char *argv[]){
	char *fa_file = NULL;
	uint64_t ********arr_check = NULL;
	char *ref_seq = NULL;
	List *these_regions = NULL;
	struct seq_region_t **ignore_regs = NULL;
	uint64_t ********covs = NULL;
	FILE *alg_bean_file = NULL;
	alg_bean_t *alg = NULL;
    FILE *config = NULL;
    int ignore_reg_count = 0;
    
    int is_err = mstep_setup_options(argc,argv);
	check(is_err==0, "Error parsing options.");

	//Open the config file and do relevant things
	config = fopen(config_file,"r");
	check(config != NULL,"Failed to open config file for reading. Have you run caveman-setup?");

	int cfg = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean,
								results,list_loc,&includeSW,&includeSingleEnd,&includeDups,version,norm_cn_loc,tum_cn_loc);

	check(strcmp(version,CAVEMAN_VERSION)==0,"Stored version in %s %s and current code version %s did not match.",config_file,version,CAVEMAN_VERSION);

	check(cfg==0,"Error parsing config file.");
  bam_access_include_sw(includeSW);
  bam_access_include_se(includeSingleEnd);
  bam_access_include_dup(includeDups);

	//Read in the alg bean
	alg_bean_file = fopen(alg_bean,"r");
	check(alg_bean_file != 0 ,"Error trying to open alg_bean file: %s.",alg_bean);
	alg = alg_bean_read_file(alg_bean_file);
	check(alg != NULL,"Error reading alg_bean from file.");
	check(fclose(alg_bean_file)==0,"Error closing alg bean file.");

	//Get split section from file given the index.
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
	//Iterate through each base in the section and create the genome profile.
	//an array to store the data full of zeroes.

	//Read order (1st/2nd)
	//Strand +/-
	//Lane
	//Read pos
	//Map qual
	//Base qual
	//Ref base
	//Called base
	covs = covs_access_generate_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	check(covs != NULL,"Error generating initial covs array");

	//Set the min base qual in case it's been changed.
	bam_access_min_base_qual(min_bq);

	//Open the bam files.
	bam_access_openbams(norm_bam_file,tum_bam_file,ref_idx);

	printf("Looking at section %s:%d-%d for mstep\n",chr_name,start_zero_based+1,stop);

	//Get ignored regions for section and calculate a list of sections to analyse.
	ignore_reg_count = ignore_reg_access_get_ign_reg_count_for_chr(ignore_regions_file,chr_name);
   check(ignore_reg_count >= 0,"Error trying to check the number of ignored regions for this chromosome.");

   //Now create a store for said regions.
   ignore_regs = malloc(sizeof(struct seq_region_t *) *  ignore_reg_count);
   check_mem(ignore_regs);
   check(ignore_reg_access_get_ign_reg_for_chr(ignore_regions_file,chr_name,ignore_reg_count,ignore_regs)==0,"Error fetching ignored regions from file.");

	//Check the contained ignored regions
	//Resolve the ignored regions and start/stop into sections for analysis.
	these_regions = ignore_reg_access_resolve_ignores_to_analysis_sections(start_zero_based+1,stop,ignore_regs,ignore_reg_count);

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

	//Iterate through sections.
	LIST_FOREACH(these_regions, first, next, cur){
		printf("M-stepping section %s:%d-%d for mstep\n",chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);
		//Get the reference sequence for this section
		ref_seq = fai_access_get_ref_seqeuence_for_pos(fa_file,chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);

		printf("fetched a reference seq of length %lu for this section.\n",strlen(ref_seq));
		check(ref_seq != NULL,"Error retrieving reference sequence for section %s:%d-%d.",chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);
		//Get all reads or pos pileups for section.
		//Iterate through positions in section and mstep
		int chk = algos_mstep_read_position(alg,covs,chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end,ref_seq, split_size);
		check(chk==0,"Error running mstep for region %s:%d-%d.",chr_name,((seq_region_t *)cur->value)->beg,((seq_region_t *)cur->value)->end);
		free(ref_seq);
	}

	//Write this split section array to file.
	char *cov_file = malloc(sizeof(char) * (strlen(results) +1 + strlen(chr_name) + 50));
	strcpy(cov_file,results);

	//Check results directory exists... if not create it.
	int dir_check = mkdir(cov_file,S_IRWXU);
	if(dir_check==0){
		printf("Created directory %s.\n",cov_file);
	}
	sprintf(cov_file,"%s/%s",cov_file,chr_name);
	//check chromosome directory exists - if not create it.
	dir_check = mkdir(cov_file,S_IRWXU);
	if(dir_check==0){
		printf("Created directory %s.\n",cov_file);
	}

	sprintf(cov_file,"%s/%d_%d.covs",cov_file,start_zero_based+1,stop);

	int res = covs_access_write_covs_to_file(cov_file,covs,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
						List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	check(res==0,"Error writing covariate array to file.");
	//Just to be paranoid, read it in again and check we've got the right results.
	arr_check = covs_access_read_covs_from_file(cov_file,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
						List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));

	int check_arrays = cov_access_compare_two_cov_arrays(covs,arr_check,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
						List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	check(check_arrays==0,"Error writing covariate array, when read from file the original and this don't match.");
	covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
			List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base), covs);
	covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
			List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base), arr_check);
	bam_access_closebams();

	free(fa_file);
	free(cov_file);
	alg_bean_destroy(alg);
	ignore_reg_access_destroy_seq_region_t_arr(ignore_reg_count,ignore_regs);
	List_clear_destroy(these_regions);
	return 0;

error:
	bam_access_closebams();
	if(arr_check) covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
			List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base), arr_check);
	if(covs) covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
			List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base), covs);
	if(fa_file) free(fa_file);
	if(ref_seq) free(ref_seq);
	if(alg_bean_file) fclose(alg_bean_file);
	if(alg) alg_bean_destroy(alg);
	ignore_reg_access_destroy_seq_region_t_arr(ignore_reg_count,ignore_regs);
	List_clear_destroy(these_regions);
	return -1;
}
