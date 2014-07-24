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

#ifndef _split_h
#define _split_h

#include <List.h>
#include <ignore_reg_access.h>

int get_read_counts_with_ignore(List *ignore,int start, int stop, char *chr);
int shrink_section_to_size(char *chr_name,int sect_start, int sect_stop, struct seq_region_t **ignore_regs, int ignore_reg_count);
int split_main(int argc, char *argv[]);

#endif
