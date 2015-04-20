/**   LICENSE
* Copyright (c) 2014-2015 Genome Research Ltd.
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
#include <fai_access.h>

char *fai_test_file = "testData/genome.fa.fai";
char *fa_test_file = "testData/genome.fa";

char *test_fai_access_get_name_from_index(){
	int index = 2;
	char *exp_chr = "II";
	int exp_length = 15279345;
	char *chr = malloc(sizeof(char) * 50);
	int length_fetched = 0;
	fai_access_get_name_from_index(index, fai_test_file, chr, &length_fetched);
	mu_assert(strcmp(chr,exp_chr)==0,"Wrong chromosome retrieved from fai file.");
	mu_assert(length_fetched == exp_length,"Wrong length retrieved from fai file.");
	exp_chr = "I";
	index = 1;
	exp_length = 15072423;
	fai_access_get_name_from_index(index, fai_test_file, chr, &length_fetched);
	mu_assert(strcmp(chr,exp_chr)==0,"Wrong chromosome retrieved from fai file.");
	mu_assert(length_fetched == exp_length,"Wrong length retrieved from fai file.");
	free(chr);
	return NULL;
}

char *test_fai_access_get_ref_seqeuence_for_pos(){
	int from = 5000;
	int to = 1510;
	char *seq_name = "I";
	char *exp_seq = "AAACTGGTTCA";
	char *seq = fai_access_get_ref_seqeuence_for_pos(fa_test_file,seq_name,from,to);
	mu_assert(seq != NULL,"NULL sequence returned.\n");
	mu_assert(strcmp(seq,exp_seq),"Sequence retrieved doesn't match expected.");
	free(seq);
	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_fai_access_get_name_from_index);
   mu_run_test(test_fai_access_get_ref_seqeuence_for_pos);
   return NULL;
}

RUN_TESTS(all_tests);
