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
#include <fai_access.h>
#include "check_fai_access_tests.h"

char *fai_test_file = "../testData/genome.fa.fai";
char *fa_test_file = "../testData/genome.fa";

START_TEST(test_fai_access_get_name_from_index){
	int index = 2;
	char *exp_chr = "II";
	int exp_length = 15279345;
	char *chr = malloc(sizeof(char) * 50);
	int length_fetched = 0;
	fai_access_get_name_from_index(index, fai_test_file, chr, &length_fetched);
	ck_assert_msg(strcmp(chr,exp_chr)==0,"Wrong chromosome retrieved from fai file.");
	ck_assert_msg(length_fetched == exp_length,"Wrong length retrieved from fai file.");
	exp_chr = "I";
	index = 1;
	exp_length = 15072423;
	fai_access_get_name_from_index(index, fai_test_file, chr, &length_fetched);
	ck_assert_msg(strcmp(chr,exp_chr)==0,"Wrong chromosome retrieved from fai file.");
	ck_assert_msg(length_fetched == exp_length,"Wrong length retrieved from fai file.");
	free(chr);
}
END_TEST

START_TEST(test_fai_access_get_ref_seqeuence_for_pos){
	int from = 5000;
	int to = 1510;
	char *seq_name = "I";
	char *exp_seq = "AAACTGGTTCA";
	char *seq = fai_access_get_ref_seqeuence_for_pos(fa_test_file,seq_name,from,to);
	ck_assert_msg(seq != NULL,"NULL sequence returned.\n");
	ck_assert_msg(strcmp(seq,exp_seq),"Sequence retrieved doesn't match expected.");
	free(seq);
}
END_TEST

Suite * check_fai_access_tests_suite(void){
    Suite *s;
    TCase *tc_fai_access;

    s = suite_create("fai_access_tests");

    /* Core test case */
    tc_fai_access = tcase_create("fai access testing");
    tcase_add_test(tc_fai_access, test_fai_access_get_name_from_index);
    tcase_add_test(tc_fai_access, test_fai_access_get_ref_seqeuence_for_pos);

    suite_add_tcase (s, tc_fai_access);

    return s;
}

