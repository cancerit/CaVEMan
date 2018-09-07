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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "dbg.h"
#include <string.h>
#include <limits.h>
#include <alg_bean.h>
#include <ignore_reg_access.h>

int ignore_reg_access_get_ign_reg_count_for_chr(char *ign_file, char *chr){
	assert(ign_file != NULL);
	assert(chr != NULL);
	FILE *file = fopen(ign_file,"r");
	check(file != NULL,"Couldn't open ignored region file: %s.",ign_file);
	//Read all lines, only including the ones that are of the correct chromosome in the count
	int entry_count = 0;
	//read and count
	char rd[200];
	while(fgets(rd, 200, file) != NULL){
		check(rd != NULL,"Invalid line read in ignored region file.");
		char *chr_nom = malloc(sizeof(char)*50);
		check_mem(chr_nom);
		int chk = sscanf(rd,"%s",chr_nom);
		check(chk == 1,"Incorrect line read.\n");
		if(strcmp(chr_nom,chr) == 0){
			entry_count++;
		}
		free(chr_nom);
	}
	check(fclose(file)==0,"Error closing ignored region file '%s'.",ign_file);
	return entry_count;
error:
	if(file) fclose(file);
	return -1;
}

seq_region_t *ignore_reg_access_get_ign_reg_overlap(int pos, struct seq_region_t **regions, int entry_count){
	int i=0;
	for(i=0; i<entry_count; i++){
		if(regions[i]->beg <= pos && regions[i]->end >= pos){
			seq_region_t *reg_copy = malloc(sizeof(struct seq_region_t));
			check_mem(reg_copy);
			reg_copy->beg = regions[i]->beg;
			reg_copy->end = regions[i]->end;
			return reg_copy;
		}
	}
error:
	return NULL;
}

int ignore_reg_access_get_ign_reg_for_chr(char *ign_file,char *chr, int entry_count, struct seq_region_t **regions){
	assert(ign_file != NULL);
	assert(chr != NULL);
	assert(entry_count >= 0);
	assert(sizeof(regions)>0);
	int is_bed = 0;
	if(entry_count == 0){
		return 0;
	}
	//assign the right size to the array
	//then reread so we can parse the actual lines.
	//Check for bed extension
	const char *ext = strrchr(ign_file, '.');
	if(ext && ext != ign_file && strcmp(ext+1,"bed")==0){
		is_bed = 1;
	}

	FILE *file = fopen(ign_file,"r");
	check(file != NULL,"Couldn't open ignored region file: %s.",ign_file);
	int found_count = 0;

	char rd[200];
	while(fgets(rd, 200, file) != NULL){
		check(rd != NULL,"Invalid line read in ignored region file.");
		char *chr_nom = malloc(sizeof(char)*50);
		check_mem(chr_nom);
		int beg,end;
		int chk = sscanf(rd,"%s\t%d\t%d",chr_nom,&beg,&end);
		if(chk==3){
			if(strcmp(chr_nom,chr) == 0){
				regions[found_count] = malloc(sizeof(struct seq_region_t));
				check_mem(regions[found_count]);
				regions[found_count]->beg = beg + is_bed;
				regions[found_count]->end = end;
				found_count++;
			}
		}else if(1==sscanf(rd,"%s",chr_nom)){//Check for just a chromosome.
			if(strcmp(chr_nom,chr) == 0){
				regions[found_count] = malloc(sizeof(struct seq_region_t));
				check_mem(regions[found_count]);
				regions[found_count]->beg = 1;
				regions[found_count]->end = INT_MAX;
				found_count++;
			}
		}else{
			free(chr_nom);
			sentinel("Incorrect line read from ignore file %s.",rd);
		}
		free(chr_nom);
	}
	check(entry_count == found_count,"Wrong number of lines found %d for chr: %s. Expected %d.",found_count,chr,entry_count);
	check(fclose(file)==0,"Error closing ignored region file '%s'.",ign_file);
	return 0;

error:
	if(file) fclose(file);
	if(regions) ignore_reg_access_destroy_seq_region_t_arr(entry_count, regions);
	return -1;
}

List *ignore_reg_access_get_ign_reg_contained(int from, int to, struct seq_region_t **regions, int entry_count){
	List *li = List_create();
	int i=0;
	for(i=0; i<entry_count; i++){
		if(regions[i]->beg >= from && regions[i]->end <= to){
			//Make a copy of this region and put in the list
			seq_region_t *reg_copy = malloc(sizeof(struct seq_region_t));
			check_mem(reg_copy);
			reg_copy->beg = regions[i]->beg;
			reg_copy->end = regions[i]->end;
			List_push(li,reg_copy);
		}
	}
	return li;
error:
	return NULL;
}

List *ignore_reg_access_resolve_ignores_to_analysis_sections(int start, int end, struct seq_region_t **regions, int entry_count){
	List *li = ignore_reg_access_get_ign_reg_contained(start,end,regions,entry_count);
	check(li != NULL,"Error fetching contained ignore regions.");

	List *reg_for_analysis = List_create();
	seq_region_t *range = malloc(sizeof(struct seq_region_t));
	range->beg = start;
	LIST_FOREACH(li, first, next, cur){
		range->end = ((seq_region_t *) cur->value)->beg - 1;
		List_push(reg_for_analysis,range);
		range = malloc(sizeof(struct seq_region_t));
		range->beg = ((seq_region_t *) cur->value)->end + 1;
	}
	range->end = end;
	List_push(reg_for_analysis,range);
	List_clear_destroy(li);
	return reg_for_analysis;
error:
	List_clear_destroy(li);
	return NULL;

}

void ignore_reg_access_destroy_seq_region_t_arr(int entry_count, seq_region_t **regions){
	if(sizeof(regions) > 0){
		int i=0;
		for(i=0;i<entry_count;i++){
			if(regions[i] != NULL){
				free(regions[i]);
			}
		}
		free(regions);
	}
}
