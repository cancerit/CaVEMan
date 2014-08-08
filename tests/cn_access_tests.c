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
#include <cn_access.h>

char *norm_cn_file = "testData/wc.cave.cn";
char *tum_cn_file = "testData/mc.cave.cn";
char *zeroes_cn_file = "testData/zeroes.cave.cn";
char *zeroes_cn_bed = "testData/zeroes.cave.cn.bed";

char *test_cn_access_get_copy_number_for_location(){
  int is_normal = 1;
  char *chr = "10";
  int pos = 17335160;
  int exp_norm = 2;
  int exp_tum = 4;
  int res = cn_access_get_copy_number_for_location(norm_cn_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem opening normal cn file file for reading.");
  mu_assert(res==exp_norm,"Wrong copy number found using normal file.");
  is_normal = 0;
  clear_copy_number_store();

  res = cn_access_get_copy_number_for_location(tum_cn_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem opening tumour cn file file for reading.");
  mu_assert(res==exp_tum,"Wrong copy number found using tumour file.");
  clear_copy_number_store();

  //Test outside the boundary of the file to assert 0 is returned.
  pos = 10000000000;
  res = cn_access_get_copy_number_for_location(tum_cn_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem opening normal cn file file for reading.");
  mu_assert(res==0,"Wrong copy number found using boundary outside the file.");

  clear_copy_number_store();
  //Try with no cn file set it should still work
  char* null_file = NULL;
  res = cn_access_get_copy_number_for_location(null_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem providing no cn file file for reading.");
  mu_assert(res==0,"Wrong copy number found using non existant file.");

  clear_copy_number_store();
  is_normal = 1;
  res = cn_access_get_copy_number_for_location(null_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem providing no cn file file for reading.");
  mu_assert(res==0,"Wrong copy number found using non existant file.");
  return NULL;
}

char *test_cn_access_get_copy_number_for_location_zeroes(){
  int is_normal = 1;
  char *chr = "1";
  int pos = 61735;
  int exp_cn = 0;
  int res = cn_access_get_copy_number_for_location(zeroes_cn_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem reading normal cn file file for reading.");
  mu_assert(res==exp_cn,"Wrong copy number found using zeroes file.");
	clear_copy_number_store();
  return NULL;
}

char *test_cn_access_get_copy_number_for_location_bed(){
  int is_normal = 1;
  char *chr = "1";
  int pos = 61734;
  int exp_cn = 4;
  int res = cn_access_get_copy_number_for_location(zeroes_cn_file,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem reading normal cn file file for reading.");
  mu_assert(res==exp_cn,"Wrong copy number found using zeroes file.");
	clear_copy_number_store();

	exp_cn = 0;
	res = cn_access_get_copy_number_for_location(zeroes_cn_bed,chr,pos,is_normal);
  mu_assert(res!=-1,"Problem reading normal cn file file for reading.");
  mu_assert(res==exp_cn,"Wrong copy number found using zeroes file.");
	clear_copy_number_store();
  return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_cn_access_get_copy_number_for_location);
   mu_run_test(test_cn_access_get_copy_number_for_location_zeroes);
   return NULL;
}

RUN_TESTS(all_tests);
