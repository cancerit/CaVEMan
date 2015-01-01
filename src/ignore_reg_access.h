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

#ifndef _ignore_reg_access_h
#define _ignore_reg_access_h

#include <stdio.h>

typedef struct seq_region_t{  
  int beg;
  int end;
  char *chr_name;
  int val;
} seq_region_t;

#define ELEMENT_TYPE seq_region_t
#define ELEMENTS_PER_NODE 8
#include <List.h>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

int ignore_reg_access_get_ign_reg_count_for_chr(char *ign_file,char *chr);
int ignore_reg_access_get_ign_reg_for_chr(char *ign_file,char *chr,int entry_count,struct seq_region_t **regions);
seq_region_t *ignore_reg_access_get_ign_reg_overlap(int pos, struct seq_region_t **regions, int entry_count);
void ignore_reg_access_destroy_seq_region_t_arr(int entry_count, struct seq_region_t **regions);
seq_region_t_List *ignore_reg_access_get_ign_reg_contained(int from, int to, struct seq_region_t **regions, int entry_count);
seq_region_t_List *ignore_reg_access_resolve_ignores_to_analysis_sections(int from, int to, struct seq_region_t **regions, int entry_count);

#endif
