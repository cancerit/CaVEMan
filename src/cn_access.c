/**   LICENSE
* Copyright (c) 2014-2015 Genome Research Ltd.
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

#include <cn_access.h>
#include <stdio.h>
#include <stdlib.h>
#include <dbg.h>
#include <List.h>
#include <assert.h>
#include <ignore_reg_access.h>

//Array of size 2 for normal and tumour.
List *cns[2] = {NULL,NULL};
static int max_cn = 10;

int cn_access_populate_cn(char *file_loc,uint8_t is_normal){
	assert(file_loc != NULL);
	if(cns[is_normal] != NULL) return 1;
	FILE *cn_file = NULL;
	int is_bed = 0;
	//Check for bed extension
	const char *ext = strrchr(file_loc, '.');
	if(ext && ext != file_loc && strcmp(ext+1,"bed")==0){
		is_bed = 1;
	}
	cn_file = fopen(file_loc,"r");
	check(cn_file>0,"Error trying to open copy number file for reading %s.",file_loc);
	List *li = List_create();
	int cop = 0;
	char rd[250];
	while(fgets(rd, 200, cn_file) != NULL){
		check(rd != NULL,"Invalid line read in ignored region file.");
		char *chr_nom = malloc(sizeof(char *));
		int beg,end;
		int chk = sscanf(rd,"%s\t%d\t%d\t%d",chr_nom,&beg,&end,&cop);
		check(chk == 4,"Incorrect line parsed from copy number file.\n");
		seq_region_t *reg = malloc(sizeof(struct seq_region_t));
		reg->chr_name = chr_nom;
		reg->beg = beg;
		if(is_bed==1){
			reg->beg = reg->beg + 1;
		}
		reg->end = end;
		if(cop > max_cn) {
			cop=max_cn;
		}
		reg->val = cop;
		List_push(li,reg);
	}
	cns[is_normal] = li;
	fclose(cn_file);
	return 1;
error:
	if(cn_file) fclose(cn_file);
	return 0;

}

int check_overlap(int reg_start, int reg_stop, uint32_t start, uint32_t stop){
	if(	(reg_start>=start && reg_stop<=stop) || //Region contained
					(start >= reg_start && start <= reg_stop) ||//overlaps start
					(stop >= reg_start && stop <= reg_stop)||//overlaps end
					(reg_start>=start && reg_start<=stop)||//Reg_start overlaps
					(reg_stop>=start && reg_stop<=stop)){//reg_end overlaps
			return 1;
	}
	return 0;
}

int8_t cn_access_get_mean_cn_for_range(char *file_loc,char *chr,uint32_t start,uint32_t stop,uint8_t is_normal){
	if(file_loc != NULL){
		int pop_cn_chk = cn_access_populate_cn(file_loc,is_normal);
		check(pop_cn_chk==1,"Error populating copy number from file %s.",file_loc);
	}

	int total_cn = 0;
	int count_cn_entries = 0;
	if(cns[is_normal] != NULL){
    LIST_FOREACH(cns[is_normal], first, next, cur){
    	int reg_start = ((seq_region_t *)cur->value)->beg;
    	int reg_stop = ((seq_region_t *)cur->value)->end;
      if(strcmp(((seq_region_t *)cur->value)->chr_name,chr) == 0
      						&& check_overlap(reg_start,reg_stop,start,stop) == 1){
        total_cn += ((seq_region_t *)cur->value)->val;
        count_cn_entries++;
      }
    }
	}

	int result = 2;//Default to cn 2
	if(count_cn_entries>0 && total_cn>0){
		result = total_cn/count_cn_entries;
	}
	return result;
error:
	return -1;

}

int8_t cn_access_get_copy_number_for_location(char *file_loc,char *chr,uint32_t pos, uint8_t is_normal){
	if(file_loc != NULL){
		int pop_cn_chk = cn_access_populate_cn(file_loc,is_normal);
		check(pop_cn_chk==1,"Error populating copy number from file %s.",file_loc);
	}else{
		return 0;
	}
	int cn = 0;
	if(cns[is_normal] != NULL){
    LIST_FOREACH(cns[is_normal], first, next, cur){
      if(strcmp(((seq_region_t *)cur->value)->chr_name,chr) == 0 && pos >= ((seq_region_t *)cur->value)->beg && pos <= ((seq_region_t *)cur->value)->end){
        cn = ((seq_region_t *)cur->value)->val;
        break;
      }
    }
	}else{
		sentinel("Somehow copy number for location %s:%d was not available.",chr,pos);
	}
	return cn;
error:
	return -1;
}

void cn_access_set_max_cn(int max_copy_number){
	max_cn = max_copy_number;
}

void clear_copy_number_store(){
	if(cns[0] != NULL){
		LIST_FOREACH(cns[0], first, next, cur){
			free(((seq_region_t *)cur->value)->chr_name);
		}
		List_clear_destroy(cns[0]);
		cns[0] = NULL;
	}

	if(cns[1] != NULL){
		LIST_FOREACH(cns[1], first, next, cur){
			free(((seq_region_t *)cur->value)->chr_name);
		}
		List_clear_destroy(cns[1]);
		cns[1] = NULL;
	}
	return;
}
