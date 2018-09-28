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
#include <math.h>
#include <stdio.h>
#include <covs_access.h>
#include "check_covs_access_tests.h"

char *test_cov_out_loc = "../testData/covs.test.out";
char *test_cov_loc = "../testData/covs.test";
char *prob_arr_loc = "../testData/probs_arr";
//covs test is a small array in the form...
//arr[0][0][0][0][0][0][0][0] = 7;
//arr[0][0][0][0][0][0][0][1] = 12;

int dim1 = 1;
int dim2 = 1;
int dim3 = 1;
int dim4 = 1;
int dim5 = 1;
int dim6 = 1;
int dim7 = 1;
int dim8 = 2;

START_TEST(test_covs_access_generate_cov_array_given_dimensions){
	uint64_t ********arr = covs_access_generate_cov_array_given_dimensions( dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	ck_assert_msg(arr != NULL,"Array not generated.\n");
	ck_assert_msg(arr[0][0][0][0][0][0][0][0] == 0, "Incorrectly formed array created.\n");
	ck_assert_msg(arr[0][0][0][0][0][0][0][1] == 0, "Incorrectly formed array created.\n");
	covs_access_free_cov_array_given_dimensions(dim1,dim2, dim3, dim4, dim5, dim6, dim7, dim8, arr);
}
END_TEST

START_TEST(test_covs_access_write_covs_to_file){
	uint64_t ********arr = covs_access_generate_cov_array_given_dimensions( dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	ck_assert_msg(arr != NULL,"Array not generated.\n");
	arr[0][0][0][0][0][0][0][0] = 7;
	arr[0][0][0][0][0][0][0][1] = 12;
	covs_access_write_covs_to_file(test_cov_out_loc,arr,dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	//void covs_access_write_covs_to_file(char *file_loc,uint64_t ********arr,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8)
	//Load in the cov file and check it matches.
	uint64_t ********arr_2 = covs_access_read_covs_from_file(test_cov_out_loc, dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	ck_assert_msg(arr_2 != NULL,"Array not generated.\n");
	ck_assert_msg(cov_access_compare_two_cov_arrays(arr,arr_2,dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8)==0,"Arrays after writing and reading in did not match.\n");
	covs_access_free_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, arr);
	covs_access_free_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, arr_2);
}
END_TEST

START_TEST(test_covs_access_read_covs_from_file){
	//uint64_t ********covs_access_read_covs_from_file(char *file_loc,int dim1,int dim2, int dim3, int dim4, int dim5, int dim6, int dim7, int dim8)
	uint64_t ********exp = covs_access_generate_cov_array_given_dimensions( dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	ck_assert_msg(exp != NULL,"Array not generated.\n");
	exp[0][0][0][0][0][0][0][0] = 7;
	exp[0][0][0][0][0][0][0][1] = 12;
	uint64_t ********got = covs_access_read_covs_from_file(test_cov_loc,dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	ck_assert_msg(got != NULL,"Array not generated.\n");
	ck_assert_msg(cov_access_compare_two_cov_arrays(exp,got,dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8)==0,"Arrays after reading in did not match.\n");
	covs_access_free_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, exp);
	covs_access_free_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, got);
}
END_TEST

START_TEST(test_cov_access_compare_two_cov_arrays){
	uint64_t ********arr = covs_access_generate_cov_array_given_dimensions( dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	arr[0][0][0][0][0][0][0][0] = 7;
	arr[0][0][0][0][0][0][0][1] = 12;
	uint64_t ********arr_2 = covs_access_generate_cov_array_given_dimensions( dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8);
	arr_2[0][0][0][0][0][0][0][0] = 7;
	arr_2[0][0][0][0][0][0][0][1] = 12;
	ck_assert_msg(cov_access_compare_two_cov_arrays(arr,arr_2,dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8)==0,"Arrays do not appear to be equal when they should be.\n");
	arr_2[0][0][0][0][0][0][0][1] = 1;
	ck_assert_msg(cov_access_compare_two_cov_arrays(arr,arr_2,dim1, dim2,  dim3,  dim4,  dim5,  dim6,  dim7,  dim8)!=0,"Arrays appear to be equal when they aren't.\n");
	covs_access_free_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, arr);
	covs_access_free_cov_array_given_dimensions(dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, arr_2);
}
END_TEST

START_TEST(test_covs_access_generate_probability_array){
	uint64_t ********arr = covs_access_generate_cov_array_given_dimensions( 1, 1,  1,  1,  1,  1,  4,  4);
	ck_assert_msg(arr != NULL,"Array not generated.\n");
	arr[0][0][0][0][0][0][0][0] = 1;
	arr[0][0][0][0][0][0][0][1] = 4;
	arr[0][0][0][0][0][0][0][2] = 5;
	arr[0][0][0][0][0][0][0][3] = 1;
	arr[0][0][0][0][0][0][1][0] = 2;
	arr[0][0][0][0][0][0][1][1] = 2;
	arr[0][0][0][0][0][0][1][2] = 2;
	arr[0][0][0][0][0][0][1][3] = 2;
	arr[0][0][0][0][0][0][2][0] = 1;
	arr[0][0][0][0][0][0][2][1] = 1;
	arr[0][0][0][0][0][0][2][2] = 4;
	arr[0][0][0][0][0][0][2][3] = 4;
	arr[0][0][0][0][0][0][3][0] = 3;
	arr[0][0][0][0][0][0][3][1] = 5;
	arr[0][0][0][0][0][0][3][2] = 0;
	arr[0][0][0][0][0][0][3][3] = 0;

	long double ********probs = covs_access_generate_probability_array(arr,1, 1,  1,  1,  1,  1,  4,  4);
	ck_assert_msg(probs != NULL,"Array not generated.\n");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][0]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][1]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][2]-logl(0.5))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][3]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][0]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][1]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][2]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][3]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][0]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][1]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][2]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][3]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][0]-logl(0.3))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][1]-logl(0.5))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][2]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][3]-logl(0.1))<0.0001,"Error in prob calculation");

	covs_access_free_cov_array_given_dimensions(1, 1,  1,  1,  1,  1,  4,  4, arr);
	covs_access_free_prob_array_given_dimensions(1, 1,  1,  1,  1,  1,  4,  4, probs);
}
END_TEST

START_TEST(test_covs_access_write_read_probs_to_file){
	uint64_t ********arr = covs_access_generate_cov_array_given_dimensions( 1, 1,  1,  1,  1,  1,  4,  4);
	ck_assert_msg(arr != NULL,"Array not generated.\n");
	arr[0][0][0][0][0][0][0][0] = 1;
	arr[0][0][0][0][0][0][0][1] = 4;
	arr[0][0][0][0][0][0][0][2] = 5;
	arr[0][0][0][0][0][0][0][3] = 1;
	arr[0][0][0][0][0][0][1][0] = 2;
	arr[0][0][0][0][0][0][1][1] = 2;
	arr[0][0][0][0][0][0][1][2] = 2;
	arr[0][0][0][0][0][0][1][3] = 2;
	arr[0][0][0][0][0][0][2][0] = 1;
	arr[0][0][0][0][0][0][2][1] = 1;
	arr[0][0][0][0][0][0][2][2] = 4;
	arr[0][0][0][0][0][0][2][3] = 4;
	arr[0][0][0][0][0][0][3][0] = 3;
	arr[0][0][0][0][0][0][3][1] = 5;
	arr[0][0][0][0][0][0][3][2] = 0;
	arr[0][0][0][0][0][0][3][3] = 0;

	char *test_out_probs = "../testData/probs.test.out";
	long double ********probs = covs_access_generate_probability_array(arr,1, 1,  1,  1,  1,  1,  4,  4);
	covs_access_write_probs_to_file(test_out_probs,probs,1, 1,  1,  1,  1,  1,  4,  4);

	ck_assert_msg(probs != NULL,"Array not generated.\n");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][0]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][1]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][2]-logl(0.5))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][0][3]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][0]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][1]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][2]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][1][3]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][0]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][1]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][2]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][2][3]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][0]-logl(0.3))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][1]-logl(0.5))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][2]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs[0][0][0][0][0][0][3][3]-logl(0.1))<0.0001,"Error in prob calculation");


	long double  ********probs_2 = covs_access_read_probs_from_file(test_out_probs,1, 1,  1,  1,  1,  1,  4,  4);
	ck_assert_msg(probs_2 != NULL, "Error reading written probs array");

	ck_assert_msg(cov_access_compare_two_prob_arrays(probs,probs_2,1, 1,  1,  1,  1,  1,  4,  4)==0,"Prob arrays after reading/writing do not match.");

	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][0][0]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][0][1]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][0][2]-logl(0.5))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][0][3]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][1][0]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][1][1]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][1][2]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][1][3]-logl(0.25))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][2][0]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][2][1]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][2][2]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][2][3]-logl(0.4))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][3][0]-logl(0.3))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][3][1]-logl(0.5))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][3][2]-logl(0.1))<0.0001,"Error in prob calculation");
	ck_assert_msg(abs(probs_2[0][0][0][0][0][0][3][3]-logl(0.1))<0.0001,"Error in prob calculation");


	covs_access_free_cov_array_given_dimensions(1, 1,  1,  1,  1,  1,  4,  4, arr);
	covs_access_free_prob_array_given_dimensions(1, 1,  1,  1,  1,  1,  4,  4, probs);
	covs_access_free_prob_array_given_dimensions(1, 1,  1,  1,  1,  1,  4,  4, probs_2);
	remove(test_out_probs);
}
END_TEST

Suite * check_covs_access_tests_suite(void){
    Suite *s;
    TCase *tc_covs_access;

    s = suite_create("covs_access_tests");

    /* Core test case */
    tc_covs_access = tcase_create("covs file access testing");

    tcase_add_test(tc_covs_access, test_covs_access_generate_cov_array_given_dimensions);
    tcase_add_test(tc_covs_access, test_covs_access_generate_probability_array);
    tcase_add_test(tc_covs_access, test_cov_access_compare_two_cov_arrays);
    tcase_add_test(tc_covs_access, test_covs_access_read_covs_from_file);
    tcase_add_test(tc_covs_access, test_covs_access_write_covs_to_file);
    tcase_add_test(tc_covs_access, test_covs_access_write_read_probs_to_file);

    suite_add_tcase (s, tc_covs_access);

    return s;
}
