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
#include <alg_bean.h>
#include <List.h>
#include <math.h>

char *norm_file = "testData/mt.bam";
char *tum_file = "testData/wt.bam";
char *test_alg_out_loc = "testData/test_alg_bean.out";
char *test_alg_loc = "testData/test_alg_bean";

static alg_bean_t *test_bean;

int check_base_lists(List *li){
	int i=0;
	LIST_FOREACH(li,first,next,cur){
		if(i==0){
			check(strcmp((char *)cur->value,"A") == 0,"Incorrect default base for entry 1.\n");
		}else if(i==1){
			check(strcmp((char *)cur->value,"C") == 0,"Incorrect default base for entry 2.\n");
		}else if(i==2){
			check(strcmp((char *)cur->value,"G") == 0,"Incorrect default base for entry 3.\n");
		}else if(i==3){
			check(strcmp((char *)cur->value,"T") == 0,"Incorrect default base for entry 4.\n");
		}else{
			sentinel("Too many entries in base List: %d\n",i+1);
		}
		i++;
	}
	return 1;
error:
	return 0;
}

int check_base_quality(List *li){
	int i=0;
	LIST_FOREACH(li,first,next,cur){
		if(i==0){
			check(((alg_bean_intrange *)cur->value)->from == 0,"Incorrect default base qual from for entry 1.\n");
			check(((alg_bean_intrange *)cur->value)->to == 9,"Incorrect default base qual to for entry 1.\n");
		}else if(i==1){
			check(((alg_bean_intrange *)cur->value)->from == 10,"Incorrect default base qual from for entry 2.\n");
			check(((alg_bean_intrange *)cur->value)->to == 19,"Incorrect default base qual to for entry 2.\n");
		}else if(i==2){
			check(((alg_bean_intrange *)cur->value)->from == 20,"Incorrect default base qual from for entry  3.\n");
			check(((alg_bean_intrange *)cur->value)->to == 24,"Incorrect default base qual to for entry 3.\n");
		}else if(i==3){
			check(((alg_bean_intrange *)cur->value)->from == 25,"Incorrect default base qual from for entry  4.\n");
			check(((alg_bean_intrange *)cur->value)->to == 29,"Incorrect default base qual to for entry 4.\n");
		}else if(i==4){
			check(((alg_bean_intrange *)cur->value)->from == 30,"Incorrect default base qual from for entry  5.\n");
			check(((alg_bean_intrange *)cur->value)->to == 34,"Incorrect default base qual to for entry 5.\n");
		}else if(i==5){
			check(((alg_bean_intrange *)cur->value)->from == 35,"Incorrect default base qual from for entry 6.\n");
			check(((alg_bean_intrange *)cur->value)->to == 39,"Incorrect default base qual to for entry 6.\n");
		}else if(i==6){
			check(((alg_bean_intrange *)cur->value)->from == 40,"Incorrect default base qual from for entry 7.\n");
			check(((alg_bean_intrange *)cur->value)->to == 200,"Incorrect default base qual to for entry 7.\n");
		}else{
			sentinel("Too many entries in base qual List: %d\n",i+1);
		}
		i++;
	}
	return 1;
error:
	return 0;
}

int check_map_quality(List *li){
	int i=0;
	LIST_FOREACH(li,first,next,cur){
		if(i==0){
			check(((alg_bean_intrange *)cur->value)->from == 0,"Incorrect default map qual from for entry 1.\n");
			check(((alg_bean_intrange *)cur->value)->to == 60,"Incorrect default map qual to for entry 1.\n");
		}else if(i==1){
			check(((alg_bean_intrange *)cur->value)->from == 255,"Incorrect default map qual from for entry 2.\n");
			check(((alg_bean_intrange *)cur->value)->to == 255,"Incorrect default map qual to for entry 2.\n");
		}else{
			sentinel("Too many entries in map qual List: %d\n",i+1);
		}
		i++;
	}
	return 1;
error:
	return 0;
}

int check_read_position(List *li){
	int i=0;
	LIST_FOREACH(li,first,next,cur){
		if(i==0){
			float exp = 2.63;
			check(fabs(*(float *)cur->value - exp) < 0.00001,
				"Incorrect default read position percentage for entry 1. %f != %f\n",*(float *)cur->value,exp);
		}else if(i==1){
			float exp = 47.37;
			check(fabs(*(float *)cur->value - exp) < 0.00001,
				"Incorrect default read position percentage for entry 2. %f != %f\n",*(float *)cur->value,exp);
		}else if(i==2){
			float exp = 23.68;
			check(fabs(*(float *)cur->value - exp) < 0.00001,
				"Incorrect default read position percentage for entry 3. %f != %f\n",*(float *)cur->value,exp);
		}else if(i==3){
			float exp = 13.16;
			check(fabs(*(float *)cur->value - exp) < 0.00001,
				"Incorrect default read position percentage for entry 4. %f != %f\n",*(float *)cur->value,exp);
		}else if(i==4){
			float exp = 13.16;
			check(fabs(*(float *)cur->value - exp) < 0.00001,
				"Incorrect default read position percentage for entry 5. %f != %f\n",*(float *)cur->value,exp);
		}else{
			sentinel("Too many entries in read pos array: %d\n",i+1);
		}
		i++;
	}
	return 1;
error:
	return 0;
}

int check_alg_bean_defaults(alg_bean_t *test_bean){
	//lane ids:	6514_2	6640_4
	//We actually have 4 lanes despite 2 being the same as they're copied headers!
	check(test_bean->lane != NULL,"Lane list unset.\n");
	LIST_FOREACH(test_bean->lane,first,next,cur){
		char *lane_id = ((char *)cur->value);
		check((strcmp(lane_id,"6514_2_1") == 0 || strcmp(lane_id,"6514_2_0") == 0 || strcmp(lane_id,"6640_4_1") == 0 || strcmp(lane_id,"6640_4_0") == 0),
								"Unexpected lane id recorded.\n");
	}
	check(List_count(test_bean->lane) == 4,"Incorrect lane size recorded in size param.\n");

	check(test_bean->rd_pos != NULL,"rd_pos list unset.\n");
	check(List_count(test_bean->rd_pos) == 5,"Incorrect read_pos array size.\n");
	check(check_read_position(test_bean->rd_pos),"Error in read_position default checks.\n");

	check(test_bean->base_qual != NULL,"base_qual list unset.\n");
	check(List_count(test_bean->base_qual) == 7,"Incorrect base_quality array size.\n");
	check(check_base_quality(test_bean->base_qual),"Error in base_quality default checks.\n");

	check(test_bean->map_qual != NULL,"map_qual list unset.\n");
	check(List_count(test_bean->map_qual) == 2,"Incorrect map_quality array size.\n");
	check(check_map_quality(test_bean->map_qual),"Error in map_quality default checks.\n");

	check(test_bean->read_order != NULL,"read_order list unset.\n");
	check(List_count(test_bean->read_order) == 2,"Incorrect read_order array size.\n");

	check(test_bean->ref_base != NULL,"ref_base list unset.\n");
	check(List_count(test_bean->ref_base) == 4,"Incorrect ref_base array size.\n");
	check(check_base_lists(test_bean->ref_base),"Error in ref_base default checks.\n");

	check(test_bean->call_base != NULL,"call_base list unset.\n");
	check(List_count(test_bean->call_base) == 4,"Incorrect call_base array size.\n");
	check(check_base_lists(test_bean->call_base),"Error in call_base default checks.\n");

	check(test_bean->strand != NULL,"strand list unset.\n");
	check(List_count(test_bean->strand) == 2,"Incorrect strand array size.\n");
	return 0;
error:
	return 1;
}

char *test_alg_bean_read_file(){
	FILE *open = fopen(test_alg_loc,"r");
	mu_assert(open != NULL, "Problem opening test alg_bean file");
	test_bean = alg_bean_read_file(open);
	fclose(open);
	mu_assert(check_alg_bean_defaults(test_bean) == 0, "Error with loaded default bean.\n");
	alg_bean_destroy(test_bean);
	return NULL;
}

char *test_alg_bean_get_index_for_str_arr(){
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(norm_file,tum_file);
	//We know the
	//6514_2_PD4107a
	//6640_4_PD4107a
	//6514_2_PD4107a
	//6640_4_PD4107a
	//Expect the following results:
	mu_assert(alg_bean_get_index_for_str_arr(test_bean->lane,"6514_2_1")==0,"Wrong index returned for lane check 1.\n");
	mu_assert(alg_bean_get_index_for_str_arr(test_bean->lane,"6640_4_1")==1,"Wrong index returned for lane check 2.\n");
	mu_assert(alg_bean_get_index_for_str_arr(test_bean->lane,"6514_2_0")==2,"Wrong index returned for lane check 3.\n");
	mu_assert(alg_bean_get_index_for_str_arr(test_bean->lane,"6640_4_0")==3,"Wrong index returned for lane check 4.\n");
	alg_bean_destroy(test_bean);
	return NULL;
}

char *test_alg_bean_get_index_for_intrange_arr(){
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(norm_file,tum_file);
	//Expect the following results:
	mu_assert(alg_bean_get_index_for_intrange_arr(test_bean->base_qual,8)==0,"Wrong index returned for base qual check 1.\n");
	mu_assert(alg_bean_get_index_for_intrange_arr(test_bean->base_qual,10)==1,"Wrong index returned for base qual check 2.\n");
	mu_assert(alg_bean_get_index_for_intrange_arr(test_bean->base_qual,25)==3,"Wrong index returned for base qual check 3.\n");
	mu_assert(alg_bean_get_index_for_intrange_arr(test_bean->base_qual,35)==5,"Wrong index returned for base qual check 4.\n");
	return NULL;
}

char *test_alg_bean_get_index_for_char_arr(){
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(norm_file,tum_file);
	//Expect the following results:
	mu_assert(alg_bean_get_index_for_char_arr(test_bean->call_base,"A")==0,"Wrong index returned for base check 1.\n");
	mu_assert(alg_bean_get_index_for_char_arr(test_bean->call_base,"C")==1,"Wrong index returned for base check 2.\n");
	mu_assert(alg_bean_get_index_for_char_arr(test_bean->call_base,"G")==2,"Wrong index returned for base check 3.\n");
	mu_assert(alg_bean_get_index_for_char_arr(test_bean->call_base,"T")==3,"Wrong index returned for base check 4.\n");
	return NULL;
}

char *test_alg_bean_get_index_for_read_pos_prop_arr(){
	test_bean = alg_bean_generate_default_alg_bean(norm_file,tum_file);
	//Expect the following results:
	mu_assert(alg_bean_get_index_for_read_pos_prop_arr(test_bean->rd_pos,2,75)==0,"Wrong index returned for RD pos check 1.\n");
	mu_assert(alg_bean_get_index_for_read_pos_prop_arr(test_bean->rd_pos,3,75)==1,"Wrong index returned for RD pos check 2.\n");
	mu_assert(alg_bean_get_index_for_read_pos_prop_arr(test_bean->rd_pos,4,150)==0,"Wrong index returned for RD pos check 3.\n");
	mu_assert(alg_bean_get_index_for_read_pos_prop_arr(test_bean->rd_pos,5,150)==1,"Wrong index returned for RD pos check 4.\n");
	return NULL;
}

char *test_alg_bean_generate_default_alg_bean(){
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(norm_file,tum_file);
	mu_assert(test_bean != NULL, "Couldn't retrieve default alg_bean.\n");

	mu_assert(check_alg_bean_defaults(test_bean) == 0, "Error with loaded default bean.\n");
	return NULL;
}

char *test_alg_bean_hard_copy_char_list(){
	test_bean = alg_bean_generate_default_alg_bean(norm_file,tum_file);
	List *one = test_bean->lane;
	List *joined = List_create();
	alg_bean_hard_copy_char_list(joined,one);
	mu_assert(List_count(one) == List_count(joined),"Wrong number of elements after hard copy");
	List_clear_destroy(one);
	List_clear_destroy(joined);
	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_alg_bean_generate_default_alg_bean);
   mu_run_test(test_alg_bean_read_file);
	mu_run_test(test_alg_bean_get_index_for_str_arr);
	mu_run_test(test_alg_bean_get_index_for_intrange_arr);
	mu_run_test(test_alg_bean_get_index_for_char_arr);
	mu_run_test(test_alg_bean_get_index_for_read_pos_prop_arr);
	mu_run_test(test_alg_bean_hard_copy_char_list);
   return NULL;
}

RUN_TESTS(all_tests);
