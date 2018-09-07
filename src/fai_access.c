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
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "dbg.h"
#include <fai_access.h>

int fai_access_get_name_from_index(int idx, char *index_file_name, char *chr_name, int *length){
	assert(index_file_name != NULL);
	assert(idx>0);
	//Open fai file
	FILE *fai = fopen(index_file_name,"r");
	check(fai != NULL,"Invalid line read\n");
	//read each chromosome until we have reached the index.
	int i=0;
	char rd[1000];
	while(fgets(rd, 1000, fai) != NULL){
		check(rd != NULL,"Invalid line read\n");
		i++;
		if(i == idx){
			int chk = sscanf(rd,"%s\t%d",chr_name,length);
			check(chk == 2,"Wrong number of entries (%d) found in fasta index file line %s",chk,rd);
			break;
		}
	}
	check(chr_name != NULL,"No line found in fai file for index %d\n", idx);
	//close file
	check(fclose(fai)==0,"Error closing fai file.");
	return 0;
error:
	if(fai)	fclose(fai);
	return -1;
}

int fai_access_get_count_length_all_contigs(char *fa_loc, int *count, int *total_len){
    char *chr_name = NULL;
    int length = 0;
    assert(fa_loc != NULL);
    //Open fai file
    chr_name = malloc(sizeof(char *));
	FILE *fai = fopen(fa_loc, "r");
	check(fai != NULL,"Invalid line read\n");
	//read each chromosome
    *count = 0;
    *total_len = 0;
    char rd[1000];
	while(fgets(rd, 1000, fai) != NULL){
		check(rd != NULL,"Invalid line read\n");
        *count = *count+1;
        fprintf(stderr,"********** %s\n", rd);
		int chk = sscanf(rd,"%s\t%d\t%*d\t%*d\t%*d",chr_name,&length);
		check(chk == 2,"Wrong number of entries (%d) found in fasta index file line %s",chk,rd);
        *total_len += strlen(chr_name);
	}
    free (chr_name);
	//close file
	check(fclose(fai)==0,"Error closing fai file.");
	return 0;
error:
	if(fai)	fclose(fai);
    if(chr_name) free(chr_name);
	return -1;
}

char *fai_access_get_ref_seqeuence_for_pos(char *fa_loc,char *char_nom,int start_one_based,int stop){
	assert(char_nom !=NULL);
	assert(fa_loc != NULL);
	char region[100];
	char *seq = NULL;
	faidx_t *fai = NULL;
	int chk = sprintf(region,"%s:%d-%d",char_nom,start_one_based,stop);
	check(chk>0,"Error formatting region.");
	int length = (stop-start_one_based)+1;
	fai = fai_load(fa_loc);
	check(fai != NULL,"Error opening FASTA index.");
	seq = fai_fetch(fai,region,&length);
	check(seq != NULL,"Error fetching reference sequence for region %s.",region);
	fai_destroy(fai);
	return seq;
error:
	if(fai) fai_destroy(fai);
	if(seq) free(seq);
	return NULL;
}
