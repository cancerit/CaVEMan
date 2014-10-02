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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <List.h>
#include <List_algos.h>
#include <dbg.h>
#include <bam_access.h>
#include <time.h>

static file_holder *norm = NULL;
static file_holder *tum = NULL;
int counter = -1;
int include_sw = 0;
int include_dup = 0;
int include_se = 0;
int min_base_qual = 10;
uint8_t isnorm = 0;
char norm_char[5];

int bam_access_openbams(char *norm_file, char *tum_file){
	assert(norm_file != NULL);
	assert(tum_file != NULL);
	//Assign memory for the file name etc holding structs
	norm = malloc(sizeof(file_holder));
	tum = malloc(sizeof(file_holder));
	check_mem(norm);
	check_mem(tum);
	//Beginning and end of tmp struct for bam access
	norm->beg = 0; norm->end = 0x7fffffff;  // The max 32 bit integer.
	tum->beg = 0; tum->end = 0x7fffffff;  // The max 32 bit integer.
	//Open a file for read from compressed bam.
	norm->in = samopen(norm_file, "rb", 0);
	check(norm->in != 0,"Normal file %s failed to open.",norm_file);
	norm->idx = bam_index_load(norm_file);
	check(norm->idx != 0,"Normal index file %s failed to open.",norm_file);
	tum->in = samopen(tum_file, "rb", 0);
	check(tum->in != 0,"Tumour file %s failed to open.",tum_file);
	tum->idx = bam_index_load(tum_file);
	check(tum->idx != 0,"Normal index file %s failed to open.",tum_file);
	return 0;
error:
	if(norm->in) samclose(norm->in);
	if(tum->in) samclose(tum->in);
	return -1;
}

// callback for bam_fetch()
static int fetch_umnorm_counts_func(const bam1_t *b, void *data){
	//check Mapping Quality and not un mapped //4 // read unmapped
	if(b->core.qual == 0 || (b->core.flag & BAM_FUNMAP)){
		return 0;
	}
	//Bad read reasons:
	//8 // mate unmapped
	//256 // Non primary alignment
	//512 // fails platform/vendor checks
	//2048 is supplementary read
	if((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FQCFAIL) || (b->core.flag & 2048)){
		return 0;
	}
	//1024 is PCR/optical duplicate
	if((b->core.flag & BAM_FDUP)){
		return 0;
	}

	//We actually want this read
	bam_plbuf_t *pileup = (bam_plbuf_t*) data;
	bam_plbuf_push(b,pileup);
  return 0;
}

// callback for bam_plbuf_init()
static int pileup_umnorm_counts(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
	file_holder *tmp = (file_holder*)data;
  if ((pos+1) > tmp->beg && pos+1 <= tmp->end) {
  	int i=0;
   	for(i=0;i<n;i++){
   		const bam_pileup1_t *p = pl + i;
   		bam1_t *algn = p->b;
			uint8_t cbase = bam1_seqi(bam1_seq(algn),p->qpos);
			if(!(p->is_del) && bam1_qual(algn)[p->qpos] >= min_base_qual && (cbase == 1 || cbase == 2 || cbase == 4 || cbase == 8)){//check bases are ACGT
				int loc = (pos) - tmp->beg;
				if(tmp->base_counts[loc] == NULL || tmp->base_counts[loc] == 0){
					tmp->base_counts[loc] = calloc(4,sizeof(int));
					check_mem(tmp->base_counts[loc]);
				}
				char called_base = toupper(bam_nt16_rev_table[cbase]);
				int x=0;
				for(x=0;x<4;x++){
					if(called_base == tmp->bam_access_bases[x]){
						tmp->base_counts[loc][x]++;
						break;
					}
				}
			}
		}//End of iteration through each pileup read in this pos.
	}
	return 0;
error:
  return -1;
}

file_holder *bam_access_get_by_position_counts(char *norm_file, char *chr, uint32_t start, uint32_t end){
	//Open bam related stuff
	assert(norm_file != NULL);
	//Assign memory for the file name etc holding structs
	norm = malloc(sizeof(file_holder));
	check_mem(norm);
	//Beginning and end of tmp struct for bam access
	//Open a file for read from compressed bam.
	norm->in = samopen(norm_file, "rb", 0);
	check(norm->in != 0,"Normal file %s failed to open.",norm_file);
	norm->idx = bam_index_load(norm_file);
	check(norm->idx != 0,"Normal index file %s failed to open.",norm_file);
	norm->bam_access_bases = malloc(sizeof(char)*4);
	check_mem(norm->bam_access_bases);
	norm->bam_access_bases[0] = 'A';
	norm->bam_access_bases[1] = 'C';
	norm->bam_access_bases[2] = 'G';
	norm->bam_access_bases[3] = 'T';
	//Pileup and populate the list with position pileup_stats.
	char *region;
	char sta[20];
	region = malloc(sizeof(chr)+sizeof(":")+sizeof("-")+(sizeof(sta)*2));
	check_mem(region);
	sprintf(region,"%s:%d-%d",chr,start,end);
	norm->beg = start;
	norm->end = end;
	//One entry for each position covered
	norm->base_counts = calloc(((end-start)+1),sizeof(int *));
	check_mem(norm->base_counts);
	norm->base_counts_size = (end-start)+1;

	int ref;
	bam_plbuf_t *buf;
	// parse the region
	bam_parse_region(norm->in->header, region, &ref,
						 &norm->beg, &norm->end);
	check(ref >= 0,"Invalid tumour region: %s.",region);

	// initialize pileup
	buf = bam_plbuf_init(pileup_umnorm_counts, norm);

  bam_fetch(norm->in->x.bam, norm->idx, ref, norm->beg, norm->end, buf, fetch_umnorm_counts_func);
	bam_plbuf_push(0, buf); // finalize pileup
  bam_plbuf_destroy(buf);
	free(region);



	//Close bam related stuff
	if(norm->idx) bam_index_destroy(norm->idx);
	if(norm->in) samclose(norm->in);
	return norm;
error:
	if(norm->idx) bam_index_destroy(norm->idx);
	if(norm->in) samclose(norm->in);
	if(norm->base_counts){
		int k=0;
		for(k=0;k<norm->base_counts_size;k++){
			if(norm->base_counts[k] != 0) free(norm->base_counts[k]);
		}
		free(norm->base_counts);
	}
	if(norm) free(norm);

	return NULL;
}

char *bam_access_sample_name_platform_from_header(char *bam_file,char *sample, char *plat){
	assert(bam_file != NULL);
	samfile_t *bam = samopen(bam_file, "rb", 0);
	check(bam != 0,"Failed to open bam file to read sample name: %s.",bam_file);
	char *head_txt = bam->header->text;
	char *line;
	line = strtok(head_txt,"\n");
	while(line != NULL){
		if(strncmp("@RG",line,3)==0){
			char *tok;
			tok = strtok(line,"\t");
			while(tok != NULL){
				if(strcmp(plat,".")==0 && strncmp("PL:",tok,(sizeof(char)*3))==0){
					int chk = sscanf(tok,"%*[^:]:%s",plat);
					check(chk==1,"Error fetching platform\n");
				}else if(strncmp("SM:",tok,(sizeof(char)*3))==0){
					int chk=sscanf(tok,"%*[^:]:%s",sample);
					check(chk==1,"Error fetching sample\n");
				}
				tok = strtok(NULL,"\t");
			}
			break;
		}
		line = strtok(NULL,"\n");
	}
	check(plat!= NULL,"Platform was not found in RG line for VCF output.");
	check(sample!= NULL,"Sample name was not found in RG line for VCF output.");
	samclose(bam);
	return "";
error:
	if(bam) samclose(bam);
	return NULL;
}

List *bam_access_get_contigs_from_bam(char *bam_file, char *assembly, char *species){
	assert(bam_file != NULL);
	char *line = NULL;
	ref_seq_t *ref = NULL;
	List *conts = List_create();
	samfile_t *bam = samopen(bam_file, "rb", 0);
	check(bam != 0,"Failed to open bam file to read contigs: %s.",bam_file);
	char *head_txt = bam->header->text;
	line = strtok(head_txt,"\n");
	while(line != NULL){
		//First check it's a sequence line
		if(strncmp(line,"@SQ",3) == 0){
			ref = malloc(sizeof(struct ref_seq_t));
			//If we already have species and assembly don't look for them
			if(assembly != NULL && species != NULL){
				//Set assembly and species
				ref->ass = assembly;
				ref->spp = species;
				//Just match name and length
				int chk = sscanf(line,"@SQ\tSN:%s\tLN:%d\t",ref->name,&ref->length);
				if(chk!=2){
					free(ref);
					sentinel("Sequence name and length not found in sequence line %s.",line);
				}
				List_push(conts,ref);
			}else{
				//Look for species and assembly
				char spec[1000];
				char assem[1000];
				int chk = sscanf(line,"@SQ\tSN:%s\tLN:%d\tAS:%[A-Za-z0-9]\tSP:%s",ref->name,&ref->length,assem,spec);
				check(chk==4,"Sequence name, length, assembly and species not found in sequence line %s.",line);
				ref->ass = malloc(sizeof(char) * 100);
				check_mem(ref->ass);
				ref->spp = malloc(sizeof(char) * 100);
				check_mem(ref->spp);
				strcpy(ref->ass,assem);
				strcpy(ref->spp,spec);
				List_push(conts,ref);
			}
		}//End of checking for sequence line
		line = strtok(NULL,"\n");
	}//End of iterating through header lines.
	check(List_count(conts)==bam->header->n_targets,"Wrong number of ref sequences in list.");
	free(line);
	samclose(bam);
	return conts;
error:
	if(line) free(line);
	if(bam) samclose(bam);
	if(ref) free(ref);
	if(conts) List_clear_destroy(conts);
	return NULL;
}

void bam_access_closebams(){
	if(norm->idx) bam_index_destroy(norm->idx);
	if(norm->in) samclose(norm->in);
	if(tum->idx) bam_index_destroy(tum->idx);
	if(tum->in) samclose(tum->in);
	if(norm) free(norm);
	if(tum) free(tum);
	return;
}

int bam_access_get_count_for_region(char *chr_name, uint32_t start, uint32_t stop){
	assert(chr_name != NULL);
	int count = bam_access_get_count_with_bam(chr_name,start,stop,tum);
	check(count >= 0,"Invalid count returned from tum.");
	int cn = bam_access_get_count_with_bam(chr_name,start,stop,norm);
	check(cn >= 0,"Invalid count returned from normal.");
	count += cn;
	return count;
error:
	return -1;
}

// callback for bam_fetch()
static int fetch_count_func(const bam1_t *b, void *data){
	//check Mapping Quality and not un mapped //4 // read unmapped
	if(b->core.qual == 0 || (b->core.flag & BAM_FUNMAP)){
		return 0;
	}
	//Bad read reasons:
	//8 // mate unmapped
	//256 // Non primary alignment
	//512 // fails platform/vendor checks
	//2048 is supplementary read
	if((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FQCFAIL) || (b->core.flag & 2048)){
		return 0;
	}
	//1024 is PCR/optical duplicate
	if((include_dup == 0 && (b->core.flag & BAM_FDUP))){
		return 0;
	}
	//Proper pair and mate unmapped
	if(!(include_se == 0 && (b->core.flag & BAM_FPROPER_PAIR) && !(b->core.flag & BAM_FMUNMAP))){
		return 0;
	}
	//printf("XT DATA: %c\n",xt);
	//Now we check aux data for XT:M flags (the SW mapped marker from BWA)
	if(include_sw == 0){
		uint8_t *xt_data = bam_aux_get(b,"XT");
	 	if(xt_data != NULL && bam_aux2A(xt_data) == 'M'){
			return 0;
		}
	}
	//checkCigar(rec)

	counter++;
   return 0;
}

static int fetch_algo_func(const bam1_t *b, void *data){
  	bam_plbuf_t *pileup = (bam_plbuf_t*) data;

  //check Mapping Quality and not un mapped //4 // read unmapped
	if(b->core.qual == 0 || (b->core.flag & BAM_FUNMAP)){
		return 0;
	}
	//Bad read reasons:
	//8 // mate unmapped
	//256 // Non primary alignment
	//512 // fails platform/vendor checks
	//2048 is supplementary read
	if((b->core.flag & BAM_FSECONDARY) || (b->core.flag & BAM_FQCFAIL) || (b->core.flag & 2048)){
		return 0;
	}
	//1024 is PCR/optical duplicate
	if((include_dup == 0 && (b->core.flag & BAM_FDUP))){
		return 0;
	}
	//Proper pair and mate unmapped
	if(!(include_se == 0 && (b->core.flag & BAM_FPROPER_PAIR) && !(b->core.flag & BAM_FMUNMAP))){
		return 0;
	}
	//printf("XT DATA: %c\n",xt);
	//Now we check aux data for XT:M flags (the SW mapped marker from BWA)
	if(include_sw == 0){
		uint8_t *xt_data = bam_aux_get(b,"XT");
	 	if(xt_data != NULL && bam_aux2A(xt_data) == 'M'){
			return 0;
		}
	}
  	bam_plbuf_push(b,pileup);
  	return 0;
}

// callback for bam_plbuf_init()
static int pileup_count_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data) {
   return 0;
}

void List_insert_sorted(List *list, void *value, List_compare cmp){
	assert(list != NULL);
	ListNode *node = calloc(1,sizeof(ListNode));
	check_mem(node);
	node->value = value;
	if(list->last == NULL){//Empty list add a single entry
		list->first = node;
		list->last = node;
		list->count++;
		return;
	}else if(cmp(list->last->value,node->value)<=0){//Checking if this entry will tag on the end, hopefully quicker than iteration.
		list->last->next = node;
		node->prev = list->last;
		list->last = node;
		list->count++;
		return;
	}else if(cmp(list->first->value,value)>=0){//If the first value is more than the value passed we tack this on the beginning
		list->first->prev = node;
		node->next = list->first;
		list->first = node;
		list->count++;
		return;
	}else{
		 // Iterate backwards through the array.
		LIST_FOREACH(list, last, prev, cur) {
			if(cmp(cur->value, value) <= 0) {//If the current node's value is less than or equal to the value passed the value is entered here.
				if(cur->next == NULL){ //Last node being changed
					list->last = node;
				}else{//If we're in the right place to put the new value
								//(the new node must be less than or equal to the next node)
					node->next = cur->next;
					node->next->prev = node;
				}
				cur->next = node;
				node->prev = cur;
				list->count++;
				return;
			}
		}
	}
error:
	return;
}

int bam_access_compare_read_pos_t(const void *in_a, const void *in_b){
	const read_pos_t *a = in_a;
	const read_pos_t *b = in_b;

	if(a->ref_pos > b->ref_pos) return 1;
	if(a->ref_pos < b->ref_pos) return -1;
	return 0;
}

static int pileup_algo_unsorted_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pil, void *data) {
	//Finally check the base quality is more than or equal to the min base quality and it's not an 'N'.
   file_holder *tmp = (file_holder*)data;
   char *nom = malloc(sizeof(char) * 350);
   check_mem(nom);
   if ((pos+1) > tmp->beg && pos+1 <= tmp->end) {
   	int i=0;
   	for(i=0;i<n;i++){
   		const bam_pileup1_t *p = pil + i;
   		bam1_t *algn = p->b;
			uint8_t cbase = bam1_seqi(bam1_seq(algn),p->qpos);
			if(!(p->is_del) && bam1_qual(algn)[p->qpos] >= min_base_qual && (cbase == 1 || cbase == 2 || cbase == 4 || cbase == 8)){//check bases are ACGT
				//Now we add a new read pos struct to the list since the read is valid.
				read_pos_t *rp = malloc(sizeof(struct read_pos_t));
				check_mem(rp);
				rp->rd_len = bam_cigar2qlen(&algn->core,bam1_cigar(algn));
				rp->ref_pos = pos+1;
				rp->rd_pos = p->qpos+1;
				rp->called_base = cbase;
				rp->map_qual = algn->core.qual;
				rp->base_qual = bam1_qual(algn)[p->qpos];
				//Check strandedness
				if(algn->core.flag & BAM_FREVERSE){
					rp->strand = 1;
					rp->rd_pos = (rp->rd_len - rp->rd_pos) + 1;
				}else{
					rp->strand = 0;
				}
				//Check read order
				if(algn->core.flag & BAM_FREAD1){
					rp->read_order = 0;
				}else if(algn->core.flag & BAM_FREAD2){
					rp->read_order = 1;
				}
				nom = strcpy(nom,bam_aux2Z(bam_aux_get(algn,"RG")));
				nom = strcat(nom,"_");
				nom = strcat(nom,norm_char);
				int lane_i = alg_bean_get_index_for_str_arr(tmp->bean->lane,nom);
				check(lane_i>=0,"Error calculating lane index %s.",nom);
				rp->lane_i = lane_i;
				rp->normal = isnorm;
				List_push(tmp->reads,rp);
			}
		}//End of iteration through each pileup read in this pos.
	}
	free(nom);
	return 0;
error:
	if(nom) free(nom);
  return 0;
}

// callback for bam_plbuf_init()
static int pileup_algo_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pil, void *data) {

   //Finally check the base quality is more than or equal to the min base quality and it's not an 'N'.
   file_holder *tmp = (file_holder*)data;
   char *nom = malloc(sizeof(char) * 350);
   if ((pos+1) > tmp->beg && (pos+1) <= tmp->end) {
   	int i=0;
   	for(i=0;i<n;i++){
   		const bam_pileup1_t *p = pil + i;
			int qual = bam1_qual(p->b)[p->qpos];
			uint8_t c = bam1_seqi(bam1_seq(p->b), p->qpos);
			if(!(p->is_del) &&  qual >= min_base_qual && (c == 1 || c == 2 || c == 4 || c == 8)){
				//Now we add a new read pos struct to the list since the read is valid.
				read_pos_t *rp = malloc(sizeof(struct read_pos_t));
				check_mem(rp);
				rp->rd_len = bam_cigar2qlen(&p->b->core,bam1_cigar(p->b));
				rp->ref_pos = (int)pos+1;
				rp->rd_pos = p->qpos+1;
				rp->called_base = c;
				rp->map_qual = p->b->core.qual;
				rp->base_qual = qual;
				//Check strandedness
				if(p->b->core.flag & BAM_FREVERSE){
					rp->rd_pos = (rp->rd_len - rp->rd_pos) + 1;
					rp->strand = 1;
				}else{
					rp->strand = 0;
				}
				//Check read order
				if(p->b->core.flag & BAM_FREAD1){
					rp->read_order = 0;
				}else if(p->b->core.flag & BAM_FREAD2){
					rp->read_order = 1;
				}
				nom = strcpy(nom,bam_aux2Z(bam_aux_get(p->b,"RG")));
				nom = strcat(nom,"_");
				nom = strcat(nom,norm_char);
				int lane_i = alg_bean_get_index_for_str_arr(tmp->bean->lane,nom);
				check(lane_i>=0,"Error calculating lane index %s.",nom);
				rp->lane_i = lane_i;
				rp->normal = isnorm;
				List_insert_sorted(tmp->reads, rp, (List_compare)bam_access_compare_read_pos_t);
			}
		}//End of iteration through each pileup read in this pos.
		free(nom);
	}
	return 0;
error:
   return 0;
}

List *bam_access_get_sorted_reads_at_this_pos(char *chr_name, uint32_t start, uint32_t stop,
																uint8_t sorted, alg_bean_t *bean, file_holder* bams, uint8_t normal){
	//Pileup and populate the list with valid reads.
	char *region = NULL;
	char sta[20];
	region = malloc(sizeof(chr_name)+sizeof(":")+sizeof("-")+(sizeof(sta)*2));
	sprintf(region,"%s:%lu-%lu",chr_name,(long unsigned int)start,(long unsigned int)stop);
	bams->beg = start;
	bams->end = stop;
	bams->reads = List_create();
	bams->bean = bean;
	int ref;
	bam_plbuf_t *buf;
	isnorm = normal;
	sprintf(norm_char,"%i",(int)normal);
	// parse the tumour region
	bam_parse_region(bams->in->header, region, &ref,
						 &bams->beg, &bams->end);
	check(ref >= 0,"Invalid tumour region: %s.",region);

	// initialize pileup
	if(sorted==1){
		buf = bam_plbuf_init(pileup_algo_func, bams);
	}else{
		buf = bam_plbuf_init(pileup_algo_unsorted_func, bams);
	}
   bam_fetch(bams->in->x.bam, bams->idx, ref, bams->beg, bams->end, buf, fetch_algo_func);
	bam_plbuf_push(0, buf); // finalize pileup
   bam_plbuf_destroy(buf);
	free(region);
	return bams->reads;
error:
	if(region) free(region);
	if(bams->reads) List_clear_destroy(bams->reads);
	return NULL;
}

List *bam_access_get_reads_at_this_pos(char *chr_name, uint32_t start, uint32_t stop, uint8_t sorted, alg_bean_t *bean){
	assert(norm->in != NULL);
	assert(tum->in != NULL);
	List *norm_rds = NULL;
	List *tum_rds = NULL;

	norm_rds = bam_access_get_sorted_reads_at_this_pos(chr_name,start,stop,sorted,bean,norm,1);
	check(norm_rds != NULL,"Error retrieving normal reads.");
	tum_rds = bam_access_get_sorted_reads_at_this_pos(chr_name,start,stop,sorted,bean,tum,0);
	check(tum_rds != NULL,"Error retrieving tumour reads.");
	//Now merge the two sorted lists into one.
	List *merged_sorted = List_merge(norm_rds,tum_rds,bam_access_compare_read_pos_t);
	check(merged_sorted != NULL,"Error merging normal and tumour reads.");
	//Free up the old lists
	free(norm_rds);
	free(tum_rds);
	return merged_sorted;
error:
	if(norm_rds) List_clear_destroy(norm_rds);
	if(tum_rds) List_clear_destroy(tum_rds);
	if(merged_sorted) List_clear_destroy(merged_sorted);
	return NULL;
}

int bam_access_get_count_with_bam(char *chr_name, uint32_t start, uint32_t stop, file_holder *fh){
	counter = 0;
	tmpstruct_t tmp;
	tmp.beg = fh->beg; tmp.end = fh->end;
	tmp.in = fh->in;
	int ref;
	bam_plbuf_t *buf = NULL;
	char *region;
	char sta[20];
	region = malloc(sizeof(chr_name)+sizeof(":")+sizeof("-")+(sizeof(sta)*2));
	sprintf(region,"%s:%d-%d",chr_name,start,stop);
	bam_parse_region(tmp.in->header, region, &ref,
						 &tmp.beg, &tmp.end); // parse the region
	check(ref >= 0,"Invalid region: %s.",region);
	free(region);
	buf = bam_plbuf_init(pileup_count_func, &tmp); // initialize pileup
   bam_fetch(tmp.in->x.bam, fh->idx, ref, tmp.beg, tmp.end, buf, fetch_count_func);
   //bam_plbuf_push(0, buf); // finalize pileup
   //Do something with the buffer.
   bam_plbuf_destroy(buf);
   return counter;
error:
	if(region) free(region);
	if(buf) bam_plbuf_destroy(buf);
	return -1;
}

void bam_access_include_sw(int inc){
	include_sw = inc;
	return;
}

void bam_access_include_dup(int inc){
	include_dup = inc;
	return;
}

void bam_access_include_se(int inc){
	include_se = inc;
	return;
}

void bam_access_min_base_qual(int qual){
	min_base_qual = qual;
	return;
}

List *bam_access_get_lane_list_from_header(char *bam_loc, char *isnorm){
	assert(bam_loc != NULL);
	char *line = NULL;
	List *li = NULL;
	samfile_t *bam =  NULL;
	bam = samopen(bam_loc, "rb", 0);
	check(bam != 0,"Bam file %s failed to open to read header.",bam_loc);
	char *head_txt = bam->header->text;
	li = List_create();
	line = strtok(head_txt,"\n");
	while(line != NULL){
		//Check for a read group line
		if(strncmp(line,"@RG",3)==0){
			char *id = malloc(sizeof(char) * 100);
			char *lane = malloc(sizeof(char) * 150);
			int chk = sscanf(line,"@RG\tID:%s",id);
			if(chk==1){
				lane = strcpy(lane,id);
				lane = strcat(lane,"_");
				lane = strcat(lane,isnorm);
				int found = 0;
				LIST_FOREACH(li, first, next, cur){
					if(strcmp((char *)cur->value,lane)==0){
						found = 1;
					}
				}
				if(found==0){
					List_push(li,lane);
				}else{
					free(lane);
					free(id);
				}
			}else{
				free(id);
				free(lane);
				sentinel("ID and SM not found in RG line %s.",line);
			}
		}
		line = strtok(NULL,"\n");
	}
	if(line) free(line);
	samclose(bam);
	return li;
error:
	if(line) free(line);
	if(bam) samclose(bam);
	if(li) List_clear_destroy(li);
	return NULL;
}
