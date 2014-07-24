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
#include <split_access.h>
#include <string.h>

char *test_split_file = "tests/split.test";

char *test_split_access_get_section_from_index(){
		//(char *file_loc, char *chr, int *start_zero_based, int *stop, int index){
	int index = 1;
	char *chr = malloc(sizeof(char) * 20);
	int start = 0;
	int stop = 0;
	split_access_get_section_from_index(test_split_file,chr,&start,&stop,index);
	mu_assert(strcmp(chr,"1")==0,"Incorrect chromsome retrieved.\n");
	mu_assert(start==0,"Incorrect start retrieved.\n");
	mu_assert(stop==10000,"Incorrect stop retrieved.\n");
	index = 2;
	free(chr);
	chr = malloc(sizeof(char) * 20);
	start = 0;
	stop = 0;
	split_access_get_section_from_index(test_split_file,chr,&start,&stop,index);
	mu_assert(strcmp(chr,"4")==0,"Incorrect chromsome retrieved 2.\n");
	mu_assert(start==11,"Incorrect start retrieved 2.\n");
	mu_assert(stop==30000,"Incorrect stop retrieved 2.\n");
	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_split_access_get_section_from_index);
   return NULL;
}

RUN_TESTS(all_tests);