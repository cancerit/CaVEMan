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

#ifndef _ignore_reg_access_h
#define _ignore_reg_access_h

#include <stdio.h>
#include <List.h>

typedef struct seq_region_t{
   int beg;
   int end;
	char *chr_name;
	int val;
} seq_region_t;

int ignore_reg_access_get_ign_reg_count_for_chr(char *ign_file,char *chr);
int ignore_reg_access_get_ign_reg_for_chr(char *ign_file,char *chr,int entry_count,struct seq_region_t **regions);
seq_region_t *ignore_reg_access_get_ign_reg_overlap(int pos, struct seq_region_t **regions, int entry_count);
void ignore_reg_access_destroy_seq_region_t_arr(int entry_count, struct seq_region_t **regions);
List *ignore_reg_access_get_ign_reg_contained(int from, int to, struct seq_region_t **regions, int entry_count);
List *ignore_reg_access_resolve_ignores_to_analysis_sections(int from, int to, struct seq_region_t **regions, int entry_count);

#endif
