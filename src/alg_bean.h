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

#ifndef _alg_bean_h
#define _alg_bean_h

#include <stdio.h>
#include <BasicLists.h>
#include "sam.h"

typedef struct alg_bean_intrange {
	int from;
	int to;
} alg_bean_intrange;

#define ELEMENT_TYPE alg_bean_intrange
#define ELEMENTS_PER_NODE 16
#include <List.h>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

typedef struct alg_bean_t {
	float_List *rd_pos;
	int rd_pos_size;
	alg_bean_intrange_List *base_qual;
	int base_qual_size;
	alg_bean_intrange_List *map_qual;
	int map_qual_size;
	String_List *lane;	
	int lane_size;
	int_List *read_order;
	int read_order_size;
	String_List *ref_base;
	int ref_base_size;
	String_List *call_base;
	int call_base_size;
	int_List *strand; 
	int strand_size;
} alg_bean_t;

int alg_bean_create_default_file(FILE *file, char *norm, char *tum);
int alg_bean_write_file(FILE *file, alg_bean_t *bean);
alg_bean_intrange_List *alg_bean_parse_int_range(char *txt);
String_List *alg_bean_parse_str_list(char *txt);
float_List *alg_bean_parse_float_list(char *txt);
alg_bean_t *alg_bean_read_file(FILE *file);
alg_bean_t *alg_bean_generate_default_alg_bean(char *norm, char *tum);
void alg_bean_destroy(alg_bean_t *bean);
String_List *alg_bean_hard_copy_char_list(String_List *new_list, String_List *old);
int alg_bean_get_index_for_str_arr(String_List *list,char *value);
int alg_bean_get_index_for_intrange_arr(alg_bean_intrange_List *list,int value);
int alg_bean_get_index_for_char_arr(String_List *list,char *value);
int alg_bean_get_index_for_read_pos_prop_arr(float_List *list,int pos,int rd_len);

#define CEIL(a, b) (((a) / (b)) + (((a) % (b)) > 0 ? 1 : 0))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define ABS(a)     (((a) < 0) ? -(a) : (a))

#endif
