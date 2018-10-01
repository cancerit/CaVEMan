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
#include <cn_access.h>

#include "check_cn_access_tests.h"

char *norm_cn_file = "../testData/wc.cave.cn";
char *tum_cn_file = "../testData/mc.cave.cn";
char *zeroes_cn_file = "../testData/zeroes.cave.cn";
char *zeroes_cn_bed = "../testData/zeroes.cave.cn.bed";

START_TEST(test_cn_access_get_copy_number_for_location){
  int is_normal = 1;
  char *chr = "10";
  int pos = 17335160;
  int exp_norm = 2;
  int exp_tum = 4;
  int res = cn_access_get_copy_number_for_location(norm_cn_file,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem opening normal cn file file for reading.");
  ck_assert_msg(res==exp_norm,"Wrong copy number found using normal file.");
  is_normal = 0;
  clear_copy_number_store();

  res = cn_access_get_copy_number_for_location(tum_cn_file,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem opening tumour cn file file for reading.");
  ck_assert_msg(res==exp_tum,"Wrong copy number found using tumour file.");
  clear_copy_number_store();

  //Test outside the boundary of the file to assert 0 is returned.
  pos = 999999999;
  res = cn_access_get_copy_number_for_location(tum_cn_file,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem opening normal cn file file for reading.");
  ck_assert_msg(res==0,"Wrong copy number found using boundary outside the file.");

  clear_copy_number_store();
  //Try with no cn file set it should still work
  char* null_file = NULL;
  res = cn_access_get_copy_number_for_location(null_file,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem providing no cn file file for reading.");
  ck_assert_msg(res==0,"Wrong copy number found using non existant file.");

  clear_copy_number_store();
  is_normal = 1;
  res = cn_access_get_copy_number_for_location(null_file,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem providing no cn file file for reading.");
  ck_assert_msg(res==0,"Wrong copy number found using non existant file.");
}
END_TEST

START_TEST(test_cn_access_get_copy_number_for_location_zeroes){
    int is_normal = 1;
    char *chr = "1";
    int pos = 61735;
    int exp_cn = 0;
    int res = cn_access_get_copy_number_for_location(zeroes_cn_file,chr,pos,is_normal);
    ck_assert_msg(res!=-1,"Problem reading normal cn file file for reading.");
    ck_assert_msg(res==exp_cn,"Wrong copy number found using zeroes file.");
    clear_copy_number_store();
}
END_TEST

START_TEST(test_cn_access_get_copy_number_for_location_bed){
  int is_normal = 1;
  char *chr = "1";
  int pos = 61734;
  int exp_cn = 4;
  int res = cn_access_get_copy_number_for_location(zeroes_cn_file,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem reading normal cn file file for reading.");
  ck_assert_msg(res==exp_cn,"Wrong copy number found using zeroes file.");
	clear_copy_number_store();

	exp_cn = 0;
	res = cn_access_get_copy_number_for_location(zeroes_cn_bed,chr,pos,is_normal);
  ck_assert_msg(res!=-1,"Problem reading normal cn file file for reading.");
  ck_assert_msg(res==exp_cn,"Wrong copy number found using zeroes file.");
	clear_copy_number_store();
}
END_TEST

START_TEST(test_cn_access_get_mean_cn_for_range){
	int exp_mean = 3;
	int start = 114629220;
	int stop = 152555165;
	char *chr = "1";
	int got = cn_access_get_mean_cn_for_range(tum_cn_file,chr,start,stop,1);
	ck_assert_msg(got==exp_mean,"Test tumour mean cn 3");
	got=0;
	got = cn_access_get_mean_cn_for_range(tum_cn_file,chr,start,stop,0);
	ck_assert_msg(got==exp_mean,"Test normal mean cn 3");

	chr = "4";
	start = 34779055;
	stop = 34824724;
	got = cn_access_get_mean_cn_for_range(tum_cn_file,chr,start,stop,0);
	ck_assert_msg(got==exp_mean, "Test normal mean cn 2");
	exp_mean = 2;
}
END_TEST

Suite * check_cn_access_tests_suite(void){
    Suite *s;
    TCase *tc_cn_access;

    s = suite_create("cn_access_tests");

    /* Core test case */
    tc_cn_access = tcase_create("cn access testing");
    tcase_add_test(tc_cn_access, test_cn_access_get_copy_number_for_location);
    tcase_add_test(tc_cn_access, test_cn_access_get_copy_number_for_location_zeroes);
    tcase_add_test(tc_cn_access, test_cn_access_get_mean_cn_for_range);

    suite_add_tcase (s, tc_cn_access);
    return s;
}
