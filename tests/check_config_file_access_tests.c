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

#include <stdlib.h>
#include <check.h>
#include <config_file_access.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <stdlib.h>

#include "check_config_file_access_tests.h"

char *cfg_file = "../testData/test.cfg";
char *cfg_file_write = "../testData/test.w.cfg";
char *tum_bam_file_exp = "../testData/mstep_test_mt.bam";
char *norm_bam_file_exp = "../testData/mstep_test_wt.bam";
char *ref_idx_exp = "../testData/ref.fai";
char *ignore_regions_file_exp = "../testData/ign.test";
char *alg_bean_loc_exp = "../alg_bean_test_new";
char *results_exp = "../results_test";
char *list_loc_exp = "../splitList_test";
char *norm_cn_exp = "../testData/wc.cave.cn";
char *tum_cn_exp = "../testData/mc.cave.cn";
int includeSW_exp = 0;
int includeSingleEnd_exp = 0;
int includeDups_exp = 0;


START_TEST(test_config_file_access_read_config_file){
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
	ck_assert_msg(config!= NULL,"Problem opening config file for reading.");
	int run = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,
                                                            alg_bean_loc,results,list_loc,&includeSW,
                                                            &includeSingleEnd,&includeDups,version,norm_cn,tum_cn);

	ck_assert_msg(run==0,"Error reading config file.");
	fclose(config);
	ck_assert_msg(strcmp(tum_bam_file,tum_bam_file_exp)==0,"Error with mut bam path %s\t%s",tum_bam_file,tum_bam_file_exp);
	ck_assert_msg(strcmp(norm_bam_file,norm_bam_file_exp)==0,"Error with norm bam path");
	ck_assert_msg(strcmp(ref_idx,ref_idx_exp)==0,"Error with ref index path");
	ck_assert_msg(strcmp(ignore_regions_file,ignore_regions_file_exp)==0,"Error with ignore file path");
	ck_assert_msg(strcmp(alg_bean_loc,alg_bean_loc_exp)==0,"Error with alg bean path");
	ck_assert_msg(strcmp(results,results_exp)==0,"Error with results file path");
	ck_assert_msg(strcmp(list_loc,list_loc_exp)==0,"Error with split file path");
	ck_assert_msg(includeSW==0,"Error with SW include");
	ck_assert_msg(includeSingleEnd==0,"Error with SE include");
	ck_assert_msg(includeDups==0,"Error with duplicates include");
	ck_assert_msg(strcmp(version,"TEST_VERSION")==0,"Error with version");
	ck_assert_msg(strcmp(norm_cn,norm_cn_exp)==0,"Error with version");
	ck_assert_msg(strcmp(tum_cn,tum_cn_exp)==0,"Error with version");
    free(tum_bam_file);
    free(norm_bam_file);
    free(ref_idx);
    free(ignore_regions_file);
    free(alg_bean_loc);
    free(results);
    free(list_loc);
    free(version);
    free(norm_cn);
    free(tum_cn);
}
END_TEST

START_TEST(test_config_file_access_write_config_file){
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
	ck_assert_msg(res==0,"Error writing config file.");
	fclose(config_w);

	FILE *config = fopen(cfg_file_write,"r");
	ck_assert_msg(config!= NULL,"Problem opening config file for reading.");
	int run = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
										results,list_loc,&includeSW,&includeSingleEnd,&includeDups,version,norm_cn,tum_cn);

	ck_assert_msg(run==0,"Error reading config file.");
	fclose(config);
	ck_assert_msg(strcmp(tum_bam_file,tum_bam_file_exp)==0,"Error with mut bam path");
	ck_assert_msg(strcmp(norm_bam_file,norm_bam_file_exp)==0,"Error with norm bam path");
	ck_assert_msg(strcmp(ref_idx,ref_idx_exp)==0,"Error with ref index path");
	ck_assert_msg(strcmp(ignore_regions_file,ignore_regions_file_exp)==0,"Error with ignore file path");
	ck_assert_msg(strcmp(alg_bean_loc,alg_bean_loc_exp)==0,"Error with alg bean path");
	ck_assert_msg(strcmp(results,results_exp)==0,"Error with results file path");
	ck_assert_msg(strcmp(list_loc,list_loc_exp)==0,"Error with split file path");
	ck_assert_msg(strcmp(version,CAVEMAN_VERSION)==0,"Error with version");
	ck_assert_msg(includeSW==0,"Error with SW include");
	ck_assert_msg(includeSingleEnd==0,"Error with SE include");
	ck_assert_msg(includeDups==0,"Error with duplicates include");
	ck_assert_msg(strcmp(norm_cn,norm_cn_exp)==0,"Error with norm CN");
	ck_assert_msg(strcmp(tum_cn,tum_cn_exp)==0,"Error with tum CN");
	int del = remove(cfg_file_write);
	ck_assert_msg(del==0,"Problem removing config file written for tests.");
	del = rmdir(results_exp);
	ck_assert_msg(del==0,"Problem removing results folder made by tests.");
    free(tum_bam_file);
    free(norm_bam_file);
    free(ref_idx);
    free(ignore_regions_file);
    free(alg_bean_loc);
    free(results);
    free(list_loc);
    free(version);
    free(norm_cn);
    free(tum_cn);
}
END_TEST

Suite * check_config_file_access_tests_suite(void){
    Suite *s;
    TCase *tc_config_file_access;

    s = suite_create("config_file_access_tests");

    /* Core test case */
    tc_config_file_access = tcase_create("config file access testing");
    tcase_add_test(tc_config_file_access, test_config_file_access_read_config_file);
    tcase_add_test(tc_config_file_access, test_config_file_access_write_config_file);

    suite_add_tcase (s, tc_config_file_access);

    return s;
}
