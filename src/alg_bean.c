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

#include <assert.h>
#include <stdlib.h>
#include <dbg.h>
#include <alg_bean.h>
#include <bam_access.h>

#define ELEMENT_TYPE alg_bean_intrange
#define ELEMENTS_PER_NODE 16
#include <List.c>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

static float f1 = 2.63f;
static float f2 = 47.37f;
static float f3 = 23.68f;
static float f4 = 13.16f;
static float f5 = 13.16f;

char *chara = "A";
char *charc = "C";
char *charg = "G";
char *chart = "T";

const char *read_pos_id = "rd_pos";
const char *base_qual_id = "base_qual";
const char *map_qual_id = "map_qual";
const char *lane_id = "lane";



int alg_bean_create_default_file(FILE *file, char *norm, char *tum){
	assert(file != NULL);
	assert(norm!= NULL);
	assert(tum!= NULL);
	alg_bean_t *bean = alg_bean_generate_default_alg_bean(norm,tum);
	int rchk = alg_bean_write_file(file, bean);
	check(rchk==0,"Error writing alg bean.");
	alg_bean_destroy(bean);
	return 0;
error:
	if(bean) alg_bean_destroy(bean);
	return 1;
}

alg_bean_t *alg_bean_read_file(FILE *file){
	assert(file != NULL);
	struct alg_bean_t *bean = malloc(sizeof(struct alg_bean_t));
	char line [ 5000 ];
	while ( fgets(line,sizeof(line),file) != NULL ){
		char text_id[50];
		char list_txt [4950];
		int chk = sscanf(line,"%s\t%s",text_id,list_txt);
		check(chk>0,"No match for alg bean entries found in line %s.",line);
		if(strcmp(text_id,read_pos_id) == 0){
			bean->rd_pos = alg_bean_parse_float_list(list_txt);
			bean->rd_pos_size = List_count(bean->rd_pos);
		}else if(strcmp(text_id,base_qual_id)==0){
			alg_bean_intrange_List *base_qual_ranges = alg_bean_parse_int_range(list_txt);
			bean->base_qual = base_qual_ranges;
			bean->base_qual_size = List_count(base_qual_ranges);
		}else if(strcmp(text_id,map_qual_id)==0){
			alg_bean_intrange_List *map_qual_ranges = alg_bean_parse_int_range(list_txt);
			bean->map_qual = map_qual_ranges;
			bean->map_qual_size = List_count(map_qual_ranges);
		}else if(strcmp(text_id,lane_id)==0){
			String_List *lanes_list = alg_bean_parse_str_list(list_txt);
			bean->lane = lanes_list;
			bean->lane_size = List_count(lanes_list);
		}else{
			sentinel("Unrecognised id passed: %s.",text_id);
		}
	}

	//Check we have all 4 ids expected fulfilled
	check(bean->rd_pos != NULL && List_count(bean->rd_pos) > 0,"No results present for read position proportions.");
	check(bean->base_qual != NULL && List_count(bean->base_qual) > 0,"No results present for base quality ranges.");
	check(bean->map_qual != NULL && List_count(bean->map_qual) > 0,"No results present for mapping quality ranges.");
	check(bean->lane != NULL && List_count(bean->lane) > 0,"No results present for lane list.");

	//Now populate the last few that aren't manually allocated.
	//Called and ref base
	//Ref and called bases
	String_List *ref_base = String_List_create();
	String_List_push(ref_base,chara);
	String_List_push(ref_base,charc);
	String_List_push(ref_base,charg);
	String_List_push(ref_base,chart);
	bean->ref_base = ref_base;
	bean->call_base = ref_base;
	bean->ref_base_size = 4;
	bean->call_base_size = 4;

	//Read strand and read order (1st/2nd in pair)
	//Strands
	int a=1;
	int b=0;
	int_List *strand = int_List_create();
	int_List_push(strand,b);
	int_List_push(strand,a);
	bean->strand = strand;
	bean->read_order = strand;
	bean->strand_size = 2;
	bean->read_order_size = 2;

	return bean;
error:
	return NULL;
}

int alg_bean_get_index_for_str_arr(String_List *list,char *val){
	int i=0;
	LIST_FOR_EACH_ELEMENT(String, list, first, next, cur) {
	    if(strcmp(cur,val)==0){
	      return i;
	    }
	    i++;
	}
	return -1;
}

int alg_bean_get_index_for_intrange_arr(alg_bean_intrange_List *list,int val){
	int i=0;
	LIST_FOR_EACH_ELEMENT(alg_bean_intrange, list, first, next, cur) {
	    if(cur.from <= val && cur.to >= val){
	      return i;
	    }
	    i++;
	}
	return -1;
}

int alg_bean_get_index_for_char_arr(String_List *list,char *val){
	int i=0;
	LIST_FOR_EACH_ELEMENT(String, list, first, next, cur) {
	    if(strcmp(cur,val) == 0){
	      return i;
	    }
	    i++;
	}
	return -1;
}

int alg_bean_get_index_for_read_pos_prop_arr(float_List *list,int pos,int rd_len){
	alg_bean_intrange_List *lengths = alg_bean_intrange_List_create();
	int i=0;
	int last_stop = 1;
	LIST_FOR_EACH_ELEMENT(float, list, first, next, pct) {
	    int lng = (((float)rd_len/(float)100) * pct);
	    alg_bean_intrange range;
	    if(i==0){
	      range.from = 1;
	    }else{
	      range.from = last_stop + 1;
	    }
	    range.to = (range.from) + lng;
	    last_stop = range.to;
	    alg_bean_intrange_List_push(lengths,range);
	    i++;
	}
	int result = alg_bean_get_index_for_intrange_arr(lengths,pos);
	alg_bean_intrange_List_destroy(lengths);
	return result;
}

String_List *alg_bean_parse_str_list(char *txt){
	String_List *li = String_List_create();
	char *ftchar;
	ftchar = strtok(txt,";");
	while(ftchar != NULL){
		char *tmp = malloc(sizeof(char) *50);
		strcpy(tmp,ftchar);
		String_List_push(li,tmp);
		ftchar = strtok(NULL,";");
	}
	if(ftchar) free(ftchar);
	return li;
}

float_List *alg_bean_parse_float_list(char *txt){
	float_List *li = float_List_create();
	char *ftchar;
	ftchar = strtok(txt,";");
	while(ftchar != NULL){
		float tmp = atof(ftchar);
		float_List_push(li,tmp);
		ftchar = strtok(NULL,";");
	}
	free(ftchar);
	return li;
}

alg_bean_intrange_List *alg_bean_parse_int_range(char *txt){
	alg_bean_intrange_List *li = alg_bean_intrange_List_create();
	char *rge;
	rge = strtok(txt,";");
	while(rge != NULL){
	        alg_bean_intrange range;
		int chk = sscanf(rge,"%d-%d",&range.from,&range.to);
		check(chk==2,"Couldn't resolve range from text %s.",txt);
		alg_bean_intrange_List_push(li,range);
		rge = strtok(NULL,";");
	}
	return li;
error:
	return NULL;
}

int alg_bean_write_list_alg_bean_intrange(FILE *file, alg_bean_intrange_List *li){
	assert(file != NULL);
	assert(li != NULL);
	assert(List_count(li) > 0);
	int chk = 0;
	LIST_FOR_EACH_ELEMENT_MORE(alg_bean_intrange, li, first, next, cur, more) {
		chk = fprintf(file,"%d-%d",cur.from,cur.to);
		check(chk>0,"Error writing intrange.");
		if (more) {
			chk = fprintf(file,";");
			check(chk>0,"Error writing string separator.");
		}
	}
	chk = fprintf(file,"\n");
	check(chk!=0,"Error writing end.");
	return 0;
error:
	return -1;
}

int alg_bean_write_list_float(FILE *file, float_List *li){
	assert(file != NULL);
	assert(li != NULL);
	assert(List_count(li) > 0);
	int chk = 0;
	LIST_FOR_EACH_ELEMENT_MORE(float, li, first, next, cur, more) {
		chk = fprintf(file,"%.2f",cur);
		check(chk>0,"Error writing float.");
		if (more) {
			chk = fprintf(file,";");
			check(chk>0,"Error writing string separator.");
		}
	}
	chk = fprintf(file,"\n");
	check(chk!=0,"Error writing end.");
	return 0;
error:
	return -1;
}

int alg_bean_write_list_char(FILE *file, String_List *li){
	assert(file != NULL);
	assert(li != NULL);
	assert(List_count(li) > 0);
	int chk = 0;
	LIST_FOR_EACH_ELEMENT_MORE(String, li, first, next, cur, more) {
		chk = fprintf(file,"%s",cur);
		check(chk>=0,"Error writing string.");
		if (more) {
			chk = fprintf(file,";");
			check(chk>0,"Error writing string separator.");
		}
	}
	chk = fprintf(file,"\n");
	check(chk!=0,"Error writing end.");
	return 0;
error:
	return -1;
}

int alg_bean_write_file(FILE *file, alg_bean_t *bean){
	assert(file != NULL);
	assert(bean != NULL);
	//Check we've got a good bean
	check(bean->rd_pos != NULL,"alg_bean had NULL read position values.");
	check(bean->base_qual != NULL,"alg_bean had NULL base quality values.");
	check(bean->map_qual != NULL,"alg_bean had NULL mapping quality values.");
	check(bean->lane != NULL,"alg_bean had NULL lane values.");
	check(bean->read_order != NULL,"alg_bean had NULL read order values.");
	check(bean->ref_base != NULL,"alg_bean had NULL reference base values.");
	check(bean->call_base != NULL,"alg_bean had NULL called base values.");
	check(bean->strand != NULL,"alg_bean had NULL strand values.");
	//Write in order;

	//Read position - NEED WORK
	//check(1==0,"Read position format requires read and write methods.");
	int chk = fprintf(file,"%s\t","rd_pos");
	check(chk!=0,"Error writing read position ID.");
	chk = alg_bean_write_list_float(file, bean->rd_pos);
	check(chk==0,"Error when writing read position ranges to file.");

	//Base quality
	chk = fprintf(file,"%s\t","base_qual");
	check(chk!=0,"Error writing base quality ID.");
	chk = alg_bean_write_list_alg_bean_intrange(file, bean->base_qual);
	check(chk==0,"Error when writing base quality ranges to file.");

	//Mapping quality
	chk = fprintf(file,"%s\t","map_qual");
	check(chk!=0,"Error writing mapping quality ID.");
	chk = alg_bean_write_list_alg_bean_intrange(file, bean->map_qual);
	check(chk==0,"Error when writing mapping quality ranges to file.");

	//Lane
	chk = fprintf(file,"%s\t","lane");
	check(chk!=0,"Error writing lane ID.");
	chk = alg_bean_write_list_char(file, bean->lane);
	check(chk==0,"Error when writing lane ranges to file.");

	return 0;
error:
	return -1;
}

void alg_bean_destroy(alg_bean_t *bean){
	if(bean != NULL){
		if(bean->rd_pos != NULL){
			float_List_destroy(bean->rd_pos);
		}
		if(bean->base_qual != NULL){
			alg_bean_intrange_List_destroy(bean->base_qual);
		}
		if(bean->map_qual != NULL){
			alg_bean_intrange_List_destroy(bean->map_qual);
		}
		if(bean->lane != NULL){
		  LIST_FOR_EACH_ELEMENT(String, bean->lane, first, next, cur) {
		    free(cur);
		  }
		  String_List_destroy(bean->lane);
		}
		if(bean->ref_base){
			String_List_destroy(bean->ref_base);
		}
	 	if(bean->strand != NULL){
	 		int_List_destroy(bean->strand);
	 	}
	 	free(bean);
	}
	return;
}

alg_bean_t *alg_bean_generate_default_alg_bean(char *norm, char *tum){
	assert(norm!= NULL);
	assert(tum!= NULL);
	struct alg_bean_t *bn = malloc(sizeof(struct alg_bean_t));

	//Ref and called bases
	String_List *ref_base = String_List_create();
	String_List_push(ref_base,"A");
	String_List_push(ref_base,"C");
	String_List_push(ref_base,"G");
	String_List_push(ref_base,"T");
	bn->ref_base = ref_base;
	bn->call_base = ref_base;
	//Strands
	int a=1;
	int b=0;
	int_List *strand = int_List_create();
	int_List_push(strand,b);
	int_List_push(strand,a);
	bn->strand = strand;

	//Base quality
	//"0-10==11-15==16-20==21-30==31-200";
	alg_bean_intrange_List *base_q_list = alg_bean_intrange_List_create();
	alg_bean_intrange range_1;
	range_1.from = 0;
	range_1.to = 9;
	alg_bean_intrange range_2;
	range_2.from = 10;
	range_2.to = 19;
	alg_bean_intrange range_3;
	range_3.from = 20;
	range_3.to = 24;
	alg_bean_intrange range_4;
	range_4.from = 25;
	range_4.to = 29;
	alg_bean_intrange range_5;
	range_5.from = 30;
	range_5.to = 34;
	alg_bean_intrange range_8;
	range_8.from = 35;
	range_8.to = 39;
	alg_bean_intrange range_9;
	range_9.from = 40;
	range_9.to = 200;
	alg_bean_intrange_List_push(base_q_list,range_1);
	alg_bean_intrange_List_push(base_q_list,range_2);
	alg_bean_intrange_List_push(base_q_list,range_3);
	alg_bean_intrange_List_push(base_q_list,range_4);
	alg_bean_intrange_List_push(base_q_list,range_5);
	alg_bean_intrange_List_push(base_q_list,range_8);
	alg_bean_intrange_List_push(base_q_list,range_9);
	bn->base_qual = base_q_list;

	//Map quality
	//0-60==255
	alg_bean_intrange_List *map_q_list = alg_bean_intrange_List_create();
	alg_bean_intrange range_6;
	range_6.from = 0;
	range_6.to = 60;
	alg_bean_intrange range_7;
	range_7.from = 255;
	range_7.to = 255;
	alg_bean_intrange_List_push(map_q_list,range_6);
	alg_bean_intrange_List_push(map_q_list,range_7);
	bn->map_qual = map_q_list;

	//Chemistries
	//Read order now, so 0,1 for 1st and 2nd;
	bn->read_order = strand;

	//Lanes
	String_List *norm_lanes = bam_access_get_lane_list_from_header(norm,"1");
	String_List *tum_lanes = bam_access_get_lane_list_from_header(tum,"0");

	String_List *joined_lanes = String_List_create();
	alg_bean_hard_copy_char_list(joined_lanes,norm_lanes);
	alg_bean_hard_copy_char_list(joined_lanes,tum_lanes);

	//split the lanes into a string array.
	bn->lane = joined_lanes;
	{LIST_FOR_EACH_ELEMENT(String, norm_lanes, first, next, cur) {
	    free(cur);
	  }}
	String_List_destroy(norm_lanes);
	{LIST_FOR_EACH_ELEMENT(String, tum_lanes, first, next, cur) {
	  free(cur);
	  }}
	String_List_destroy(tum_lanes);

	//Read position
	float_List *float_list = float_List_create();
	float_List_push(float_list,f1);
	float_List_push(float_list,f2);
	float_List_push(float_list,f3);
	float_List_push(float_list,f4);
	float_List_push(float_list,f5);
	bn->rd_pos = float_list;
	/*alg_bean_intrange **rd_pos;
	 int rd_pos_size; */

	//Alg bean store file will be in the format NAMEOFCORVARIATE\trange1;range2;range3
	//Excluding read pos which needs to be relative to read length.
	//called base and ref base and read order and strand are not included as they're assumed to be in separate ranges anyway.

	return bn;
}

String_List *alg_bean_hard_copy_char_list(String_List *new_list, String_List *old){
  LIST_FOR_EACH_ELEMENT(String, old, first, next, cur) {
    char *tmp = malloc(sizeof(char) * (strlen(cur)+1));
    strcpy(tmp,cur);
    String_List_push(new_list,tmp);
  }
  return new_list;
}
