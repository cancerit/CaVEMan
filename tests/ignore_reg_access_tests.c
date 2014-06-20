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
#include <ignore_reg_access.h>

char *test_ign_file = "tests/ign.test";
//Above file contains:
//1	200	5000
//10	13	5678
//10	5679	5699

char *test_ignore_reg_access_get_ign_reg_count_for_chr(){
	char *chr = "1";
	int res = ignore_reg_access_get_ign_reg_count_for_chr(test_ign_file, chr);
	mu_assert(res==1,"Incorrect count of ignored regions returned.\n");
	chr = "10";
	res = ignore_reg_access_get_ign_reg_count_for_chr(test_ign_file, chr);
	mu_assert(res==2,"Incorrect count of ignored regions returned.\n");
	return NULL;
}

char *test_ignore_reg_access_get_ign_reg_overlap(){
	//seq_region_t *ignore_reg_access_get_ign_reg_overlap(int pos, struct seq_region_t **regions, int entry_count)
	return NULL;
}

char *test_ignore_reg_access_get_ign_reg_for_chr(){
	//void ignore_reg_access_get_ign_reg_for_chr(char *ign_file,char *chr, int entry_count, struct seq_region_t **regions)
	return NULL;
}

char *test_ignore_reg_access_get_ign_reg_contained(){
	//List *ignore_reg_access_get_ign_reg_contained(int from, int to, struct seq_region_t **regions, int entry_count)
	return NULL;
}

char *test_ignore_reg_access_resolve_ignores_to_analysis_sections(){
	//List *ignore_reg_access_resolve_ignores_to_analysis_sections(int start, int end, struct seq_region_t **regions, int entry_count)
	return NULL;
}

char *test_ignore_reg_access_destroy_seq_region_t_arr(){
	//void ignore_reg_access_destroy_seq_region_t_arr(int entry_count, seq_region_t **regions)
	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_ignore_reg_access_get_ign_reg_count_for_chr);
   return NULL;
}

RUN_TESTS(all_tests);