/**   LICENSE
* Copyright (c) 2014 Genome Research Ltd. 
* 
* Author: Cancer Genome Project cgpit@sanger.ac.uk 
* 
* This file is part of caveman_c. 
* 
* caveman_c is free software: you can redistribute it and/or modify it under 
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
#include <config_file_access.h>
#include <string.h>

char *cfg_file = "tests/test.cfg";
char *cfg_file_write = "tests/test.w.cfg";
char *tum_bam_file_exp = "tum_test.bam";
char *norm_bam_file_exp = "norm_test.bam";
char *ref_idx_exp = "index_test.fa.fai";
char *ignore_regions_file_exp = "ignore_test.tab";
char *alg_bean_loc_exp = "alg_bean_test";
char *results_exp = "results_test";
char *list_loc_exp = "splitList_test";
int includeSW_exp = 0;
int includeSingleEnd_exp = 0;
int includeDups_exp = 0;


char *test_config_file_access_read_config_file(){
	char *tum_bam_file = malloc(sizeof(char) * 50);
	char *norm_bam_file = malloc(sizeof(char) * 50);
	char *ref_idx = malloc(sizeof(char) * 50);
	char *ignore_regions_file = malloc(sizeof(char) * 50);
	char *alg_bean_loc = malloc(sizeof(char) * 50);
	char *results = malloc(sizeof(char) * 50);
	char *list_loc = malloc(sizeof(char) * 50);
	int includeSW = 1;
	int includeSingleEnd = 1;
	int includeDups = 1;
	
	FILE *config = fopen(cfg_file,"r");
	mu_assert(config!= NULL,"Problem opening config file for reading.");
	int run = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
										results,list_loc,&includeSW,&includeSingleEnd,&includeDups);
										
	mu_assert(run==0,"Error reading config file.");
	fclose(config);
	mu_assert(strcmp(tum_bam_file,tum_bam_file_exp)==0,"Error with mut bam path");
	mu_assert(strcmp(norm_bam_file,norm_bam_file_exp)==0,"Error with norm bam path");
	mu_assert(strcmp(ref_idx,"index_test.fa.fai")==0,"Error with ref index path");
	mu_assert(strcmp(ignore_regions_file,"ignore_test.tab")==0,"Error with ignore file path");
	mu_assert(strcmp(alg_bean_loc,"alg_bean_test")==0,"Error with alg bean path");
	mu_assert(strcmp(results,"results_test")==0,"Error with results file path");
	mu_assert(strcmp(list_loc,"splitList_test")==0,"Error with split file path");
	mu_assert(includeSW==0,"Error with SW include");
	mu_assert(includeSingleEnd==0,"Error with SE include");
	mu_assert(includeDups==0,"Error with duplicates include");
	
	return NULL;
}

char *test_config_file_access_write_config_file(){
	char *tum_bam_file = malloc(sizeof(char) * 50);
	char *norm_bam_file = malloc(sizeof(char) * 50);
	char *ref_idx = malloc(sizeof(char) * 50);
	char *ignore_regions_file = malloc(sizeof(char) * 50);
	char *alg_bean_loc = malloc(sizeof(char) * 50);
	char *results = malloc(sizeof(char) * 50);
	char *list_loc = malloc(sizeof(char) * 50);
	int includeSW = 1;
	int includeSingleEnd = 1;
	int includeDups = 1;
	FILE *config_w = fopen(cfg_file_write,"w");
	int res = config_file_access_write_config_file(config_w,tum_bam_file_exp,norm_bam_file_exp,ref_idx_exp,ignore_regions_file_exp,
						alg_bean_loc_exp,results_exp,list_loc_exp,includeSW_exp,includeSingleEnd_exp,includeDups_exp);
	mu_assert(res==0,"Error writing config file.");
	fclose(config_w);
	
	FILE *config = fopen(cfg_file_write,"r");
	mu_assert(config!= NULL,"Problem opening config file for reading.");
	int run = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
										results,list_loc,&includeSW,&includeSingleEnd,&includeDups);
										
	mu_assert(run==0,"Error reading config file.");
	fclose(config);
	
	mu_assert(strcmp(tum_bam_file,"tum_test.bam")==0,"Error with mut bam path");
	mu_assert(strcmp(norm_bam_file,"norm_test.bam")==0,"Error with norm bam path");
	mu_assert(strcmp(ref_idx,"index_test.fa.fai")==0,"Error with ref index path");
	mu_assert(strcmp(ignore_regions_file,"ignore_test.tab")==0,"Error with ignore file path");
	mu_assert(strcmp(alg_bean_loc,"alg_bean_test")==0,"Error with alg bean path");
	mu_assert(strcmp(results,"results_test")==0,"Error with results file path");
	mu_assert(strcmp(list_loc,"splitList_test")==0,"Error with split file path");
	mu_assert(includeSW==0,"Error with SW include");
	mu_assert(includeSingleEnd==0,"Error with SE include");
	mu_assert(includeDups==0,"Error with duplicates include");
	
	int del = remove(cfg_file_write);
	mu_assert(del==0,"Problem removing config file written for tests.");
	
	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_config_file_access_read_config_file);
   mu_run_test(test_config_file_access_write_config_file);
   return NULL;
}

RUN_TESTS(all_tests);