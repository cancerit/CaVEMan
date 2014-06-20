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

#ifndef _split_access_h
#define _split_access_h

#include <stdio.h>
#include <List.h>

int split_access_print_section(FILE *output, char *chr, int start_zero_based, int stop);
void split_access_get_section_from_index(char *file_loc, char *chr, int *start_zero_based, int *stop, int index);
List *split_access_get_all_split_sections(char *file_loc);

#endif