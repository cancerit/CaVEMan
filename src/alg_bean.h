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

#ifndef _alg_bean_h
#define _alg_bean_h

#include <stdio.h>
#include "List.h"
#include "khash.h"

//New hash to store unique readlengths
KHASH_MAP_INIT_INT(readlenpos, List *)

typedef struct alg_bean_intrange{
	int from;
	int to;
} alg_bean_intrange;

typedef struct alg_bean_t{
	 List *rd_pos;
	 int rd_pos_size;
	 List *base_qual;
	 int base_qual_size;
	 List *map_qual;
	 int map_qual_size;
	 List *lane;
	 int lane_size;
	 List *read_order;
	 int read_order_size;
	 List *ref_base;
	 int ref_base_size;
	 List *call_base;
	 int call_base_size;
	 List *strand;
	 int strand_size;
     khash_t(readlenpos) *read_len_pos;
} alg_bean_t;

int alg_bean_create_default_file(FILE *file, char *norm, char *tum);
int alg_bean_write_file(FILE *file, alg_bean_t *bean);
List *alg_bean_parse_int_range(char *txt);
List *alg_bean_parse_str_list(char *txt);
List *alg_bean_parse_float_list(char *txt);
alg_bean_t *alg_bean_read_file(FILE *file);
alg_bean_t *alg_bean_generate_default_alg_bean(char *norm, char *tum);
void alg_bean_destroy(alg_bean_t *bean);
List *alg_bean_hard_copy_char_list(List *new_list, List *old);
int alg_bean_get_index_for_str_arr(List *list,char *value);
int alg_bean_get_index_for_intrange_arr(List *list,int value);
int alg_bean_get_index_for_char_arr(List *list,char *value);
int alg_bean_get_index_for_read_pos_prop_arr(void *_hash, int pos,int rd_len);
int alg_bean_add_read_length_arrs(alg_bean_t *bean, char* list_loc, char* contig);
List *alg_bean_get_position_list_from_read_pos_proportion_arr(List *list,int rd_len);

#define CEIL(a, b) (((a) / (b)) + (((a) % (b)) > 0 ? 1 : 0))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define ABS(a)     (((a) < 0) ? -(a) : (a))

#endif
