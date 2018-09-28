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
#include <split_access.h>
#include <string.h>
#include "check_split_access_tests.h"

char *test_split_file = "../testData/split.test";

START_TEST(test_split_access_get_section_from_index){
		//(char *file_loc, char *chr, int *start_zero_based, int *stop, int index){
	int index = 1;
	char *chr = malloc(sizeof(char) * 20);
	int start = 0;
	int stop = 0;
	split_access_get_section_from_index(test_split_file,chr,&start,&stop,index);
	ck_assert_msg(strcmp(chr,"1")==0,"Incorrect chromsome retrieved.\n");
	ck_assert_msg(start==0,"Incorrect start retrieved.\n");
	ck_assert_msg(stop==10000,"Incorrect stop retrieved.\n");
	index = 2;
	free(chr);
	chr = malloc(sizeof(char) * 20);
	start = 0;
	stop = 0;
	split_access_get_section_from_index(test_split_file,chr,&start,&stop,index);
	ck_assert_msg(strcmp(chr,"4")==0,"Incorrect chromsome retrieved 2.\n");
	ck_assert_msg(start==11,"Incorrect start retrieved 2.\n");
	ck_assert_msg(stop==30000,"Incorrect stop retrieved 2.\n");
    free(chr);
}
END_TEST

Suite * check_split_access_tests_suite(void){
    Suite *s;
    TCase *tc_split_access;

    s = suite_create("split_access_tests");

    /* Core test case */
    tc_split_access = tcase_create("split access testing");

    tcase_add_test(tc_split_access, test_split_access_get_section_from_index);

    suite_add_tcase (s, tc_split_access);

    return s;
}
