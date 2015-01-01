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

#include <stdio.h>
#include <assert.h>
#include <dbg.h>
#include <alg_bean.h>
#include <getopt.h>
#include <fai_access.h>
#include <math.h>
#include <file_tests.h>
#include <ignore_reg_access.h>
#include <split_access.h>
#include <covs_access.h>
#include <config_file_access.h>

static char tum_bam_file[512];// = NULL;
static char norm_bam_file[512];// = NULL;
static char results[512];// = "results";
static char *config_file = "caveman.cfg.ini";
static char ref_idx[512];// = "";
static char list_loc[512];// = "splitList";
static char alg_bean_loc[512];// = "alg_bean";
static char version[50];
static char *covariate_file = "covs_arr";
static char *probs_file = "probs_arr";
static long double maxNonRefProbability = 0.25;
static char norm_cn_loc[512];
static char tum_cn_loc[512];

void merge_print_usage (int exit_code){
	printf ("Usage: caveman merge [-f file] [-c file] [-p file] \n\n");
	printf("Optional\n");
	printf("-f  --config-file file               Path to the config file produced by setup. [default: '%s']\n\n",config_file);
	printf("-c  --covariate-file filename        Location to write merged covariate array [default: %s]\n",covariate_file);
	printf("-p  --probabilities-file filename    Location to write probability array [default: %s]\n",probs_file);
	printf("-h	help                             Display this usage information.\n");
  exit(exit_code);
}

void merge_setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
					{"config-file", required_argument, 0, 'f'},
             	{"covariate-file", required_argument, 0, 'c'},
             	{"probabilities-file", required_argument, 0, 'p'},
             	{"help", no_argument, 0, 'h'},
             	{ NULL, 0, NULL, 0}
   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "p:c:f:hb",long_opts, &index)) != -1){
   	switch(iarg){
   		case 'h':
         	merge_print_usage(0);
         	break;

         case 'c':
         	covariate_file = optarg;
         	break;

         case 'p':
         	probs_file = optarg;
         	break;

      	case 'f':
      		config_file = optarg;
      		break;

			case '?':
            merge_print_usage(1);
            break;

      	default:
      		merge_print_usage(1);

   	}; // End of args switch statement

   }//End of iteration through options

   //Do some checking to ensure required arguments were passed and are accessible files
   if(check_exist(config_file) != 1){
   	printf("Config file %s does not appear to exist. Have you run caveman setup?\n",config_file);
   	merge_print_usage(1);
   }

   return;
}

int merge_main(int argc, char *argv[]){

	merge_setup_options(argc,argv);

	alg_bean_t *alg = NULL;
	uint64_t ********covs = NULL;
	long double ********prob_arr = NULL;
	int tmp = 0;
	//Open the config file and do relevant things
	FILE *config = fopen(config_file,"r");
	check(config != NULL,"Failed to open config file for reading. Have you run caveman setup?");

	char ignore_regions_file[512];

	int cfg = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
								results,list_loc,&tmp,&tmp,&tmp,version,norm_cn_loc,tum_cn_loc);

	check(strcmp(version,CAVEMAN_VERSION)==0,"Stored version in %s %s and current code version %s did not match.",config_file,version,CAVEMAN_VERSION);

	check(cfg==0,"Error parsing config file.");

	//Generate an alg bean so we know the sizes...
	FILE *alg_file = fopen(alg_bean_loc,"r");
	alg = alg_bean_read_file(alg_file);
	//Create an empty covariate array.
	covs = covs_access_generate_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	//Get all split sections from split file.
	List *split_sects = split_access_get_all_split_sections(list_loc);
	//Iterate through split sections and retrieve each cov array file
	LIST_FOR_EACH_ELEMENT(split_sects, first, next, cur) {
		char *chr = ((seq_region_t *)cur)->chr_name;
		int start = ((seq_region_t *)cur)->beg;
		int stop = ((seq_region_t *)cur)->end;
		char cov_loc[500] = "";
		int chck = sprintf(cov_loc,"%s/%s/%d_%d.covs",results,chr,start,stop);
		check(chck>0,"Error creating path to covs file.");
		uint64_t ********to_append = covs_access_read_covs_from_file(cov_loc,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
		//append retrieved array to the new one
		covs_access_merge_count_arrays(covs,to_append,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
		covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),to_append);
		free(((seq_region_t *)cur)->chr_name);
	}
	//Destroy the list of split sections.
	List_clear_destroy(split_sects);

	//Write the merged cov array to file.
	int chk_covs = covs_access_write_covs_to_file(covariate_file,covs,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	check(chk_covs==0,"Error writing covariate array to file.");
	//read in the written file and check we've still got the right results.
	uint64_t ********arr_check = covs_access_read_covs_from_file(covariate_file,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));

	int check_arrays = cov_access_compare_two_cov_arrays(covs,arr_check,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	check(check_arrays==0,"Error writing covariate array, when read back in they don't match.");


	//Use the merged cov array to create the probabilities array.
	prob_arr = covs_access_generate_probability_array(covs,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));


	int i,j,k,m,n,p,r,s = 0;
	for(i=0;i<List_count(alg->read_order);i++){
		for(j=0;j<List_count(alg->strand);j++){
			for(k=0;k<List_count(alg->lane);k++){
				for(m=0;m<List_count(alg->rd_pos);m++){
					for(n=0;n<List_count(alg->map_qual);n++){
						for(p=0;p<List_count(alg->base_qual);p++){
							for(r=0;r<List_count(alg->ref_base);r++){
								for(s=0;s<List_count(alg->call_base);s++){
									if(r!=s && expl(prob_arr[i][j][k][m][n][p][r][s]) > maxNonRefProbability){
										fprintf(stderr,"WARNING: At a covariate location where ref base != called base a probability %5.1Le higher than the default %5.1Le was detected.\n",expl(prob_arr[i][j][k][m][n][p][r][s]),maxNonRefProbability);
										fprintf(stderr,"This is likely to cause erroneous genotype calls and may be a symptom of some form of incorrect quality manipulation.\n");
									}
								}
							}
						}
					}
				}
			}
		}
	}

	//Write the probability array to file.
	int check_probs = covs_access_write_probs_to_file(probs_file,prob_arr,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	check(check_probs==0,"Error writing probability array to file.");
	//Free up the arrays & alg bean.
	covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),covs);
	covs_access_free_prob_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
					List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),prob_arr);
	alg_bean_destroy(alg);
	//Done
	return 0;
error:
	if(covs) covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),covs);
	if(prob_arr) covs_access_free_prob_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
						List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),prob_arr);
	if(alg) alg_bean_destroy(alg);
	return -1;
}
