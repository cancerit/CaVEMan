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

#ifndef _bam_access_h
#define _bam_access_h

#include <stdio.h>
#include <List.h>
#include <alg_bean.h>
#include <List_algos.h>
#include "sam.h"

typedef struct file_holder{
   int beg, end;
	samfile_t *in;
	bam_index_t *idx;
	List *reads;
	alg_bean_t *bean;
	int **base_counts;
	int base_counts_size;
	char *bam_access_bases;
} file_holder;

typedef struct{
	int beg, end;
	samfile_t *in;
} tmpstruct_t;

typedef struct ref_seq_t{
	int length;
	char *ass;
	char name[100];
	char *spp;
} ref_seq_t;

typedef struct read_pos_t{
	char called_base;
	int rd_pos;
	int base_qual;
	int map_qual;
	int lane_i;
	int read_order;
	int strand;
	int ref_pos;
	int rd_len;
	int normal;
	//char *rd_nom;
	long double ref_base_probs[4];
} read_pos_t;

int bam_access_openbams(char *normFile, char *tumFile);

int bam_access_get_count_for_region(char *chr_name, int start, int stop);

void bam_access_closebams();

int bam_access_get_count_with_bam(char *chr_name, int start, int stop, file_holder *fh);

void bam_access_include_sw(int inc);

void bam_access_include_dup(int inc);

void bam_access_include_se(int inc);

void bam_access_min_base_qual(int qual);

List *bam_access_get_lane_list_from_header(char *bam_file_loc, char *isnorm);

List *bam_access_get_reads_at_this_pos(char *chr_name, int start, int stop, int sorted, alg_bean_t *bean);

char *bam_access_sample_name_platform_from_header(char *bam_file,char *sample, char *plat);

void List_insert_sorted(List *list, void *value, List_compare cmp);

List *bam_access_get_contigs_from_bam(char *bam_file, char *assembly, char *species);

file_holder *bam_access_get_by_position_counts(char *normFile, char *chr, int start, int end);

#endif
