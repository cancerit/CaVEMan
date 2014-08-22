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
#include <config_file_access.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>

char *cfg_file = "./testData/test.cfg";
char *cfg_file_write = "./testData/test.w.cfg";
char *tum_bam_file_exp = "./testData/mstep_test_mt.bam";
char *norm_bam_file_exp = "./testData/mstep_test_wt.bam";
char *ref_idx_exp = "./testData/ref.fai";
char *ignore_regions_file_exp = "./testData/ign.test";
char *alg_bean_loc_exp = "./alg_bean_test_new";
char *results_exp = "./results_test";
char *list_loc_exp = "./splitList_test";
char *norm_cn_exp = "./testData/wc.cave.cn";
char *tum_cn_exp = "./testData/mc.cave.cn";
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
	char *version = malloc(sizeof(char) * 50);
	char *norm_cn = malloc(sizeof(char) * 50);
	char *tum_cn = malloc(sizeof(char) * 50);
	int includeSW = 1;
	int includeSingleEnd = 1;
	int includeDups = 1;

	FILE *config = fopen(cfg_file,"r");
	mu_assert(config!= NULL,"Problem opening config file for reading.");
	int run = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
										results,list_loc,&includeSW,&includeSingleEnd,&includeDups,version,norm_cn,tum_cn);

	mu_assert(run==0,"Error reading config file.");
	fclose(config);
	mu_assert(strcmp(tum_bam_file,tum_bam_file_exp)==0,"Error with mut bam path");
	mu_assert(strcmp(norm_bam_file,norm_bam_file_exp)==0,"Error with norm bam path");
	mu_assert(strcmp(ref_idx,ref_idx_exp)==0,"Error with ref index path");
	mu_assert(strcmp(ignore_regions_file,ignore_regions_file_exp)==0,"Error with ignore file path");
	mu_assert(strcmp(alg_bean_loc,alg_bean_loc_exp)==0,"Error with alg bean path");
	mu_assert(strcmp(results,results_exp)==0,"Error with results file path");
	mu_assert(strcmp(list_loc,list_loc_exp)==0,"Error with split file path");
	mu_assert(includeSW==0,"Error with SW include");
	mu_assert(includeSingleEnd==0,"Error with SE include");
	mu_assert(includeDups==0,"Error with duplicates include");
	mu_assert(strcmp(version,"TEST_VERSION")==0,"Error with version");
	mu_assert(strcmp(norm_cn,norm_cn_exp)==0,"Error with version");
	mu_assert(strcmp(tum_cn,tum_cn_exp)==0,"Error with version");

	return NULL;
}

char *test_config_file_access_write_config_file(){
	char *tum_bam_file = malloc(sizeof(char) * PATH_MAX);
	char *norm_bam_file = malloc(sizeof(char) * PATH_MAX);
	char *ref_idx = malloc(sizeof(char) * PATH_MAX);
	char *ignore_regions_file = malloc(sizeof(char) * PATH_MAX);
	char *alg_bean_loc = malloc(sizeof(char) * PATH_MAX);
	char *results = malloc(sizeof(char) * PATH_MAX);
	char *list_loc = malloc(sizeof(char) * PATH_MAX);
	char *version = malloc(sizeof(char) * PATH_MAX);
	char *norm_cn = malloc(sizeof(char) * PATH_MAX);
	char *tum_cn = malloc(sizeof(char) * PATH_MAX);
	int includeSW = 1;
	int includeSingleEnd = 1;
	int includeDups = 1;
	FILE *config_w = fopen(cfg_file_write,"w");
	int res = config_file_access_write_config_file(config_w,tum_bam_file_exp,norm_bam_file_exp,ref_idx_exp,ignore_regions_file_exp,
						alg_bean_loc_exp,results_exp,list_loc_exp,includeSW_exp,includeSingleEnd_exp,includeDups_exp,norm_cn_exp,tum_cn_exp);
	mu_assert(res==0,"Error writing config file.");
	fclose(config_w);

	FILE *config = fopen(cfg_file_write,"r");
	mu_assert(config!= NULL,"Problem opening config file for reading.");
	int run = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
										results,list_loc,&includeSW,&includeSingleEnd,&includeDups,version,norm_cn,tum_cn);

	mu_assert(run==0,"Error reading config file.");
	fclose(config);
	char *tmp = malloc(sizeof(char) * PATH_MAX);
	char *chk = realpath(tum_bam_file_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(tum_bam_file,tmp)==0,"Error with mut bam path");
	chk = realpath(norm_bam_file_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(norm_bam_file,tmp)==0,"Error with norm bam path");
	chk = realpath(ref_idx_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(ref_idx,tmp)==0,"Error with ref index path");
	chk = realpath(ignore_regions_file_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(ignore_regions_file,tmp)==0,"Error with ignore file path");
	chk = realpath(alg_bean_loc_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(alg_bean_loc,tmp)==0,"Error with alg bean path");
	chk = realpath(results_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(results,tmp)==0,"Error with results file path");
	chk = realpath(list_loc_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(list_loc,tmp)==0,"Error with split file path");
	mu_assert(strcmp(version,CAVEMAN_VERSION)==0,"Error with version");
	mu_assert(includeSW==0,"Error with SW include");
	mu_assert(includeSingleEnd==0,"Error with SE include");
	mu_assert(includeDups==0,"Error with duplicates include");
	chk = realpath(norm_cn_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(norm_cn,tmp)==0,"Error with norm CN");
	chk = realpath(tum_cn_exp,tmp);
	mu_assert(chk!=NULL,"realpath allocation");
	mu_assert(strcmp(tum_cn,tmp)==0,"Error with tum CN");

	int del = remove(cfg_file_write);
	mu_assert(del==0,"Problem removing config file written for tests.");
	del = remove(alg_bean_loc_exp);
	mu_assert(del==0,"Problem removing alg_bean made by tests.");
	del = remove(list_loc_exp);
	mu_assert(del==0,"Problem removing splitList made by tests.");
	del = rmdir(results_exp);
	mu_assert(del==0,"Problem removing results folder made by tests.");

	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_config_file_access_read_config_file);
   mu_run_test(test_config_file_access_write_config_file);
   return NULL;
}

RUN_TESTS(all_tests);
