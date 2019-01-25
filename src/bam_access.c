/**   LICENSE
* Copyright (c) 2014-2016 Genome Research Ltd.
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
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <List.h>
#include <List_algos.h>
#include <dbg.h>
#include <bam_access.h>
#include <time.h>
#include "khash.h"

//New hash to store readname -> strand -> base
KHASH_SET_INIT_STR(rdnom)
KHASH_MAP_INIT_STR(strh,uint8_t)
KHASH_MAP_INIT_INT(rpos, read_pos_t *)
KHASH_MAP_INIT_INT(strd, uint8_t)
KHASH_MAP_INIT_STR(rdnom_strd, khash_t(strd) *)
KHASH_MAP_INIT_STR(rdnom_rp, khash_t(rpos) *)

static file_holder *norm = NULL;
static file_holder *tum = NULL;
int counter = -1;
int include_sw = 0;
int include_dup = 0;
int include_se = 0;
int min_base_qual = 10;
uint8_t isnorm = 0;
char norm_char[5];
int maxitercnt = 1000000000;

int bam_access_openbams(char *norm_file, char *tum_file, char *ref_file){
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
	norm->in = hts_open(norm_file, "r");
	check(norm->in != 0,"Normal file %s failed to open.",norm_file);
	norm->idx = sam_index_load(norm->in,norm_file);
	check(norm->idx != 0,"Normal index file %s failed to open.",norm_file);
	norm->head = sam_hdr_read(norm->in);
	hts_set_fai_filename(norm->in, ref_file);
	tum->in = hts_open(tum_file, "r");
	check(tum->in != 0,"Tumour file %s failed to open.",tum_file);
	tum->idx = sam_index_load(tum->in,tum_file);
	check(tum->idx != 0,"Normal index file %s failed to open.",tum_file);
	tum->head = sam_hdr_read(tum->in);
	hts_set_fai_filename(tum->in, ref_file);
	return 0;
error:
	if(norm->in) hts_close(norm->in);
	if(tum->in) hts_close(tum->in);
	if(tum->head) bam_hdr_destroy(tum->head);
	if(norm->head) bam_hdr_destroy(norm->head);
	if(norm->idx) hts_idx_destroy(norm->idx);
	if(tum->idx) hts_idx_destroy(tum->idx);
	return -1;
}

// callback for bam_plbuf_init()
static int pileup_blank(void *data, bam1_t *b) {
  return 0;
}

void bam_access_copy_read_pos(read_pos_t *old, read_pos_t *new){
    assert(old != NULL && new != NULL);
    new->ref_pos = old->ref_pos;
	new->rd_len = old->rd_len;
	new->normal = old->normal;
	new->read_order = old->read_order;
	new->strand = old->strand;
	new->called_base = old->called_base;
	new->rd_pos = old->rd_pos;
	new->base_qual = old->base_qual;
	new->map_qual = old->map_qual;
	new->lane_i = old->lane_i;
    return;
}

int bam_access_get_avg_readlength_from_bam(htsFile *sf){
  assert(sf != NULL);
  int read_count = 0;
  int32_t read_length_sum = 0;
  //Read first 100 lines of file and calculate average read length.
  bam_hdr_t *head = NULL;
	head = sam_hdr_read(sf);
  bam1_t *b = bam_init1();
  int ret;
  while ((ret = sam_read1(sf, head, b)) >= 0 && read_count < 100) {
		if((b->core.flag & BAM_FSECONDARY)
		|| (b->core.flag & BAM_FSUPPLEMENTARY)
		|| (b->core.flag & BAM_FQCFAIL)){
			continue;
		}
    read_count++;
    read_length_sum += b->core.l_qseq;
  }
  bam_destroy1(b);
  return (int)(read_length_sum/read_count);
}

int pos_counts_callback_old(uint32_t tid, uint32_t pos, int n_plp, const bam_pileup1_t *pil, void *data, int strand){
  file_holder *norm  = (file_holder *) data;
	khash_t(strh) *h;
	khiter_t k;
	h = kh_init(strh);
  int i=0;
  for(i=0;i<n_plp;i++){
    const bam_pileup1_t *p = pil + i;
    bam1_t *algn = p->b;

		int absent;
    uint8_t cbase = bam_seqi(bam_get_seq(algn),p->qpos);
		k = kh_put(strh, h, bam_get_qname(p->b), &absent);
		uint8_t pre_b;
		if(!absent){ //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
			k = kh_get(strh, h, bam_get_qname(p->b));
			pre_b = kh_val(h,k);
		}else{
			//Add the value to the hash
			kh_value(h, k) = cbase;
		}

    if(!(p->is_del) && bam_get_qual(algn)[p->qpos] >= min_base_qual && (cbase == 1 || cbase == 2 || cbase == 4 || cbase == 8)
															&& (absent || pre_b != cbase)){//check bases are ACGT and not same base in overlapping read]]
      int loc = (pos + 1) - norm->beg;
      if(norm->base_counts[loc] == NULL || norm->base_counts[loc] == 0){
        if(strand == 1 ){
          norm->base_counts[loc] = calloc(8,sizeof(int));
        }else{
          norm->base_counts[loc] = calloc(4,sizeof(int));
        }
        check_mem(norm->base_counts[loc]);
      }
      int is_rev = bam_is_rev(algn);
      char called_base = toupper(seq_nt16_str[cbase]);
      int x=0;
      for(x=0;x<4;x++){
        if(called_base == norm->bam_access_bases[x]){
          if(strand == 1 && is_rev == 1){
            norm->base_counts[loc][x+4]++; //Bump this to the reverse strand counts when we are using strand
          }else{
            norm->base_counts[loc][x]++;
          }
          break;
        }
      }
    }//End of checking base is ACGT & fits quality requirements
  }//Iteration through pileups at this position
	kh_destroy(strh, h);
  return 0;
error:
	kh_destroy(strh, h);
  return 1;
}

int pos_counts_callback(uint32_t tid, uint32_t pos, int n_plp, const bam_pileup1_t *pil, void *data, int strand){
    file_holder *norm  = (file_holder *) data;
    int *local_counts = NULL;
    khash_t(rdnom_strd) *rdnom_h;
	khiter_t k;
	rdnom_h = kh_init(rdnom_strd);

    int i=0;
    for(i=0;i<n_plp;i++){

        const bam_pileup1_t *p = pil + i;
        bam1_t *algn = p->b;

        int rd_name_missing;
        char *readname = bam_get_qname(p->b);
        //Get the index of the base at this position in this read
        uint8_t cbase = seq_nt16_int[bam_seqi(bam_get_seq(algn), p->qpos)];

        // Check whether this is a valid position
        if(!(p->is_del) && bam_get_qual(algn)[p->qpos] >= min_base_qual && cbase < 4 ) {
            k = kh_put(rdnom_strd, rdnom_h, readname, &rd_name_missing);
            khash_t(strd) *strd_h;
            //Check to see if this readname is already a key
            if(rd_name_missing){ // If the readname key doesn't yet exist
                strd_h = kh_init(strd);
                kh_value(rdnom_h, k) = strd_h;
            } else {
                //Retrieve the existing strand hash stored under the readname key
                k = kh_get(rdnom_strd, rdnom_h, readname);
                strd_h = kh_val(rdnom_h, k);
            }

            //Check that the strand doesn't already exist for this readname
            //In theory it is never possible for the strand to occur twice
            int strand_missing;
            int is_rev = bam_is_rev(algn);
            k = kh_put(strd, strd_h, is_rev, &strand_missing);
            check(strand_missing==1, "Repeated strand found for theis read %s.", readname);
            if(strand == 1 && is_rev == 1) cbase += 4; //Store a reverse strand count index
            kh_value(strd_h, k) = cbase; //Set the index to the value of this base
        }// End of checking if this is a valid base to use in the read

    }//End of iteration through each pileup opject

    //Initialise the base counts array. Local and stored
    int loc = (pos + 1) - norm->beg;
    if(norm->base_counts[loc] == NULL || norm->base_counts[loc] == 0){
        if(strand == 1 ){
            norm->base_counts[loc] = calloc(8,sizeof(int));
        }else{
            norm->base_counts[loc] = calloc(4,sizeof(int));
        }
        check_mem(norm->base_counts[loc]);
    }
    local_counts = calloc(4,sizeof(int));
    check_mem(local_counts);

    khash_t(strd) *strd_h;
    // Now we iterate through every read group hash entry 
    // Destroy the sub hashes as we go so cleanup is easier.
    kh_foreach_value(rdnom_h, strd_h, {
        int strand_fwd_missing = 0;
        int strand_rev_missing = 0;
        khiter_t k_fwd;
        khiter_t k_rev;

        //Check for forward strand
        k_fwd = kh_get(strd, strd_h, 0);
        strand_fwd_missing = (k_fwd == kh_end(strd_h));
        //Check for reverse strand
        k_rev = kh_get(strd, strd_h, 1);
        strand_rev_missing = (k_rev == kh_end(strd_h));

        //If we only have one or other strand
        if(strand_fwd_missing){
            uint8_t index = kh_val(strd_h, k_rev);
            //Add to base counts
            norm->base_counts[loc][index]++;
        } else if (strand_rev_missing){
            uint8_t index = kh_val(strd_h, k_fwd);
            //Add to base counts
            norm->base_counts[loc][index]++;
        } else { //We have both strands present at this position.
            uint8_t index_fwd = kh_val(strd_h, k_fwd);
            uint8_t index_rev = kh_val(strd_h, k_rev);
            uint8_t rev_corrected = index_rev;
            if(index_rev > 3){
                rev_corrected = index_rev - 4;
            }
            //If the two bases differ we add both to the counts
            if(rev_corrected != index_fwd){
                norm->base_counts[loc][index_rev]++;
                norm->base_counts[loc][index_fwd]++;
            }else{
                //Otherwise we add only one instance. Orientation decided later
                local_counts[index_fwd]++;
            }

        } //End of checking which strands we have for each readname
    }
    ); // End of kh_foreach. Iteration through each readname hash keys

    //Now we iterate through the local base counts and add according to what we have
    int x=0;
    for(x=0; x<4; x++){

        //While we have local_counts
        while(local_counts[x] > 0){
            if(strand == 1){ // If we have stranded counting
                //If the fwd strand is higher than the rev
                if(norm->base_counts[loc][x] > norm->base_counts[loc][x+4]){
                    norm->base_counts[loc][x+4]++;
                }else{
                    norm->base_counts[loc][4]++;
                }
            }else{ // not looking at stranded counting so we can ignore the fwd/rev
                norm->base_counts[loc][4]++;
            }
            local_counts[x]--;
        } // End of while we have local counts for this base

    } //End of iteration through each 

    //Iterate through the readname hash and clearup all hashes as values
    for (k = 0; k < kh_end(rdnom_h); ++k){
        if (kh_exist(rdnom_h, k)){
            kh_destroy(strd, kh_val(rdnom_h, k));
        }
    }
    //Cleanup the readname hash

    free(local_counts);
    kh_destroy(rdnom_strd, rdnom_h);
    return 0;

error:
    if(local_counts) free(local_counts);
    //Cleanup all hashes
    if(rdnom_h){
        for (k = 0; k < kh_end(rdnom_h); ++k){
            if (kh_exist(rdnom_h, k)){
                kh_destroy(strd, kh_val(rdnom_h, k));
            }
        }
        kh_destroy(rdnom_strd, rdnom_h);
    }

    return 1;
}

file_holder *bam_access_get_by_position_counts_stranded(char *norm_file, char *chr, uint32_t start, uint32_t end, int strand){
  //Open bam related stuff
	assert(norm_file != NULL);
	bam1_t *b = NULL;
	bam_plp_t buf = NULL;
	hts_itr_t *iter = NULL;
	//Assign memory for the file name etc holding structs
	norm = malloc(sizeof(file_holder));
	check_mem(norm);
	//Beginning and end of tmp struct for bam access
	//Open a file for read from compressed bam.
	norm->in = hts_open(norm_file, "r");
	check(norm->in != 0,"Normal file %s failed to open.",norm_file);
	norm->head = sam_hdr_read(norm->in);
	norm->idx = sam_index_load(norm->in,norm_file);
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
	region = malloc((sizeof((strlen(chr)+1)*sizeof(char)))+sizeof(":")+sizeof("-")+(sizeof(sta)*2));
	check_mem(region);
	sprintf(region,"%s:%d-%d",chr,start,end);
	norm->beg = start;
	norm->end = end;
	//One entry for each position covered
	norm->base_counts = calloc(((end-start)+1),sizeof(int *));
	check_mem(norm->base_counts);
	norm->base_counts_size = (end-start)+1;

	// initialize pileup

	buf = bam_plp_init(pileup_blank, (void *)norm);
	bam_plp_set_maxcnt(buf,maxitercnt);

  iter = sam_itr_querys(norm->idx, norm->head, region);

  b = bam_init1();

  int result;
  int tid, pos, n_plp = -1;
  const bam_pileup1_t *pil;
  while ((result = sam_itr_next(norm->in, iter, b)) >= 0) {
    if(b->core.qual == 0
      || (b->core.flag & BAM_FPROPER_PAIR) != BAM_FPROPER_PAIR
      || (b->core.flag & BAM_FUNMAP)
			|| (b->core.flag & BAM_FSECONDARY)
			|| (b->core.flag & BAM_FQCFAIL)
			|| (b->core.flag & BAM_FSUPPLEMENTARY)
			|| (b->core.flag & BAM_FDUP) ) continue;
    //Additional check for paired end orientation of reads (for this fix we asume paire end)
    if(!(b->core.flag & BAM_FREVERSE) == !(b->core.flag & BAM_FMREVERSE)) continue;
    bam_plp_push(buf, b);
    while ( (pil=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
      if(!((pos+1) >= norm->beg && (pos+1) <= norm->end)) continue;
      int res = pos_counts_callback(tid, pos, n_plp, pil, norm, strand);
      check(res==0,"Error running pileup callback");
    }//End of while we have pileup reads
  }
  if(result != -1){
    fprintf(stderr,"SAMTOOLS ERROR %d\n",result);
  }
	bam_plp_push(buf,0); // finalize pileup
  sam_itr_destroy(iter);
  //Now for the pileup method
  while ( (pil=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
    if(!((pos+1) >= norm->beg && (pos+1) <= norm->end)) continue;
    int res = pos_counts_callback(tid, pos, n_plp, pil, norm, strand);
    check(res==0,"Error running pileup callback");
  }//End of while we have pileup reads
  bam_destroy1(b);
  bam_plp_destroy(buf);

	free(region);

	//Close bam related stuff
	if(norm->idx) hts_idx_destroy(norm->idx);
	if(norm->in) hts_close(norm->in);
	return norm;
error:
	if(norm->idx) hts_idx_destroy(norm->idx);
	if(norm->in) hts_close(norm->in);
	if(iter) sam_itr_destroy(iter);
	if(b) bam_destroy1(b);
	if(norm->base_counts){
		int k=0;
		for(k=0;k<norm->base_counts_size;k++){
			if(norm->base_counts[k] != 0) free(norm->base_counts[k]);
		}
		free(norm->base_counts);
	}
	return NULL;
}

file_holder *bam_access_get_by_position_counts_with_strand(char *normFile, char *chr, uint32_t start, uint32_t end){
  return bam_access_get_by_position_counts_stranded(normFile, chr, start, end, 1);
}

file_holder *bam_access_get_by_position_counts(char *norm_file, char *chr, uint32_t start, uint32_t end){
	return bam_access_get_by_position_counts_stranded(norm_file, chr, start, end, 0);
}

char *bam_access_sample_name_platform_from_header(char *bam_file,char *sample, char *plat){
	assert(bam_file != NULL);
	htsFile *bam = hts_open(bam_file, "r");
	check(bam != 0,"Failed to open bam file to read sample name: %s.",bam_file);
	bam_hdr_t *header = sam_hdr_read(bam);
	char *head_txt = header->text;
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
	hts_close(bam);
	return "";
error:
	if(bam) hts_close(bam);
	return NULL;
}

int bam_access_parse_sq_line(char *line,char *species, char *assembly,char *name,uint32_t *length){
  char *tag = strtok(line,"\t");
  char *tmp = NULL;
  while(tag != NULL){
    int chk=0;
    tmp = malloc(sizeof(char) * 512);
    check_mem(tmp);
    chk = sscanf(tag,"SN:%[^\t\n]",tmp);
    if(chk>0){
      strcpy(name,tmp);
      tag = strtok(NULL,"\t");
      continue;
    }
    chk = sscanf(tag,"SP:%[^\t\n]",tmp);
    if(chk>0){
      strcpy(species,tmp);
      tag = strtok(NULL,"\t");
      continue;
    }
    chk = sscanf(tag,"AS:%[^\t\n]",tmp);
    if(chk>0){
      strcpy(assembly,tmp);
      tag = strtok(NULL,"\t");
      continue;
    }
    chk = sscanf(tag,"LN:%d",length);
    if(chk>0){
      tag = strtok(NULL,"\t");
      continue;
    }
    tag = strtok(NULL,"\t");
  }//End of tokenising SQ line
  free(tmp);
  return 0;
error:
  if(tmp) free(tmp);
  return 1;
}

List *bam_access_get_contigs_from_bam(char *bam_file, char *assembly, char *species){
	assert(bam_file != NULL);
	char *line = NULL;
	ref_seq_t *ref = NULL;
	char ** ptr = NULL;
	List *conts = List_create();
	char *tmp_line = NULL;
	htsFile *bam = hts_open(bam_file, "r");
	check(bam != 0,"Failed to open bam file to read contigs: %s.",bam_file);
	bam_hdr_t *header = sam_hdr_read(bam);
	char *head_txt = header->text;
	ptr = malloc(sizeof(char **));
  check_mem(ptr);
	line = strtok_r(head_txt,"\n",ptr);
	while(line != NULL){
		//First check it's a sequence line
		if(strncmp(line,"@SQ",3) == 0){
			ref = malloc(sizeof(struct ref_seq_t));
			check_mem(ref);
			ref->length = 0;
			tmp_line = malloc(sizeof(char) * (strlen(line)+1));
			check_mem(tmp_line);
			strcpy(tmp_line,line);
			//If we already have species and assembly don't look for them
			if(assembly != NULL && species != NULL){
				//Set assembly and species
				ref->ass = assembly;
				ref->spp = species;
				char *dummy = malloc(sizeof(char) * 100);
				check_mem(dummy);
				bam_access_parse_sq_line(tmp_line,dummy,dummy,ref->name,&(ref->length));
			}else{
				//Look for species and assembly as well as name and length
				ref->ass = malloc(sizeof(char) * 100);
				check_mem(ref->ass);
				ref->spp = malloc(sizeof(char) * 100);
				check_mem(ref->spp);
        bam_access_parse_sq_line(tmp_line,ref->spp,ref->ass,ref->name,&(ref->length));
			}
			check(ref->name!=NULL,"Sequence name not found/set in SQ line %s",line);
			check(ref->ass!=NULL,"Sequence assembly not found/set in SQ line %s",line);
			check(ref->spp!=NULL,"Sequence species not found/set in SQ line %s",line);
			check(ref->length>0,"Sequence length not found/set in SQ line %s",line);
			List_push(conts,ref);
		}//End of checking for sequence line
		line = strtok_r(NULL,"\n",ptr);
	}//End of iterating through header lines.
	check(List_count(conts)==header->n_targets,"Wrong number of ref sequences in list.");
	free(line);
	free(ptr);
	hts_close(bam);
	free(tmp_line);
	return conts;
error:
	if(line) free(line);
	if(bam) hts_close(bam);
	if(ptr) free(ptr);
	if(ref) free(ref);
	if(tmp_line) free(tmp_line);
	if(conts) List_clear_destroy(conts);
	return NULL;
}

void bam_access_closebams(){
	if(norm->idx) hts_idx_destroy(norm->idx);
	if(norm->in) hts_close(norm->in);
	if(tum->idx) hts_idx_destroy(tum->idx);
	if(tum->in) hts_close(tum->in);
	if(tum->head) bam_hdr_destroy(tum->head);
	if(norm->head) bam_hdr_destroy(norm->head);
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

int bam_access_check_bam_flags(const bam1_t *b){
	//check Mapping Quality and not un mapped //4 // read unmapped
	if(b->core.qual == 0 || (b->core.flag & BAM_FUNMAP)){
		return 0;
	}
	//Bad read reasons:
	//8 // mate unmapped
	//256 // Non primary alignment
	//512 // fails platform/vendor checks
	//2048 is supplementary read
	if((b->core.flag & BAM_FSECONDARY)
	    || (b->core.flag & BAM_FQCFAIL)
	    || (b->core.flag & BAM_FSUPPLEMENTARY)){
		return 0;
	}
	//1024 is PCR/optical duplicate
	if((include_dup == 0
	    && (b->core.flag & BAM_FDUP))){
		return 0;
	}
	//Proper pair and mate unmapped
	if(include_se == 0
	      && !((b->core.flag & BAM_FPROPER_PAIR)
	    && !(b->core.flag & BAM_FMUNMAP))){
		return 0;
	}
  if(include_se == 0){ // If we're looking for proper pairs only
    if(! (b->core.flag & BAM_FPROPER_PAIR)) return 0;
    if(!(b->core.flag & BAM_FREVERSE) == !(b->core.flag & BAM_FMREVERSE)) return 0;
  }
	//printf("XT DATA: %c\n",xt);
	//Now we check aux data for XT:M flags (the SW mapped marker from BWA)
	if(include_sw == 0){
		uint8_t *xt_data = bam_aux_get(b,"XT");
	 	if(xt_data != NULL && bam_aux2A(xt_data) == 'M'){
			return 0;
		}
	}
	return 1;
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

int reads_at_pos_callback_old(uint32_t tid, uint32_t pos, int n_plp, const bam_pileup1_t *pil, void *data, int sorted, int isnorm){
  file_holder* bams = (file_holder* )data;
	khash_t(strh) *h;
	khiter_t k;
  char *nom = malloc(sizeof(char) * 350);
	h = kh_init(strh);
  int i=0;
  for(i=0;i<n_plp;i++){
		const bam_pileup1_t *p = pil + i;
    int qual = bam_get_qual(p->b)[p->qpos];
    uint8_t c = bam_seqi(bam_get_seq(p->b), p->qpos);

		int absent;
    k = kh_put(strh, h, bam_get_qname(p->b), &absent);
		uint8_t pre_b;
		if(!absent){ //Read already processed to get base processed (we only increment if base is different between overlapping read pairs)
			k = kh_get(strh, h, bam_get_qname(p->b));
			pre_b = kh_val(h,k);
		}else{
			//Add the value to the hash
			kh_value(h, k) = c;
		}

    if(!(p->is_del) &&  qual >= min_base_qual && (c == 1 || c == 2 || c == 4 || c == 8)&& (absent || pre_b != c)){
      //Now we add a new read pos struct to the list since the read is valid.
      read_pos_t *rp = malloc(sizeof(struct read_pos_t));
      check_mem(rp);
      rp->rd_len = bam_cigar2qlen(p->b->core.n_cigar,bam_get_cigar(p->b));
      rp->ref_pos = pos+1;
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
      int lane_i = alg_bean_get_index_for_str_arr(bams->bean->lane,nom);
      check(lane_i>=0,"Error calculating lane index %s.",nom);
      rp->lane_i = lane_i;
      rp->normal = isnorm;
      if(sorted==1){
        List_insert_sorted(bams->reads, rp, (List_compare)bam_access_compare_read_pos_t);
      }else{
        List_push(bams->reads,rp);
      }
    }//End of if this is a useful read, ACGT, and within qual boundaries.
  }//End iterating through pileups at this position
  free(nom);
	kh_destroy(strh, h);
  return 0;
error:
	kh_destroy(strh, h);
  return 1;
}

int reads_at_pos_callback(uint32_t tid, uint32_t pos, int n_plp, const bam_pileup1_t *pil, void *data, int sorted, int isnorm){
    file_holder* bams = (file_holder* ) data;
    int *local_counts = NULL;
    char *nom = NULL;
    khash_t(rdnom_rp) *h_rd_nom;
    khiter_t k;
    local_counts = calloc(8,sizeof(int));
    check_mem(local_counts);
    nom = malloc(sizeof(char) * 350);
    check_mem(nom);
    h_rd_nom = kh_init(rdnom_rp);
    char *readname = NULL;

    int i=0;
    for(i=0;i<n_plp;i++){
        const bam_pileup1_t *p = pil + i;
        char *readname = bam_get_qname(p->b);
        int qual = bam_get_qual(p->b)[p->qpos];
        uint8_t c = bam_seqi(bam_get_seq(p->b), p->qpos);

        //Check if this is a useable position
        if(!(p->is_del) && qual >= min_base_qual && seq_nt16_int[c] < 4){
            //Populate all the read position information
            read_pos_t *rp = malloc(sizeof(struct read_pos_t));
            check_mem(rp);
            rp->rd_len = bam_cigar2qlen(p->b->core.n_cigar,bam_get_cigar(p->b));
            rp->ref_pos = pos+1;
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
            int lane_i = alg_bean_get_index_for_str_arr(bams->bean->lane,nom);
            check(lane_i>=0,"Error calculating lane index %s.",nom);
            rp->lane_i = lane_i;
            rp->normal = isnorm;

            // Now check to see if the readname already exists in the hash
            int rd_name_missing;
            khash_t(rpos) *h_strd;
            k = kh_put(rdnom_rp, h_rd_nom, readname, &rd_name_missing);
            if(rd_name_missing){
                h_strd = kh_init(rpos);
                kh_value(h_rd_nom, k) = h_strd;
            }else{
                k = kh_get(rdnom_rp, h_rd_nom, readname);
                h_strd = kh_val(h_rd_nom, k);
            }

            // We have the strand subhash. 
            // Check to see if we already have a strand entry (should never happen)
            int strand_missing;
            k = kh_put(rpos, h_strd, rp->strand, &strand_missing);
            check(strand_missing==1, "Error, strand was already found for this readname %s.", readname);
            kh_value(h_strd, k) = rp; // Add the read position to the subhash

        } // End of if this is a usable position

    }//End of iteration through each pileup object

    khash_t(rpos) *h_strd;
    kh_foreach_value(h_rd_nom, h_strd, {
        int strand_fwd_missing = 0;
        int strand_rev_missing = 0;
        khiter_t k_fwd;
        khiter_t k_rev;

        k_fwd = kh_get(rpos, h_strd, 0);
        strand_fwd_missing = (k_fwd == kh_end(h_strd));
        k_rev = kh_get(rpos, h_strd, 1);
        strand_rev_missing = (k_rev == kh_end(h_strd));

        // If either of the strands are missing we can automatically include them
        if (strand_fwd_missing) {
            read_pos_t *rp_rev = malloc(sizeof(struct read_pos_t));
            check_mem(rp_rev);
            bam_access_copy_read_pos(kh_val(h_strd,k_rev), rp_rev);
            kh_val(h_strd,k_rev) = NULL;
            uint8_t c_idx = seq_nt16_int[rp_rev->called_base];
            c_idx = c_idx+4;
            local_counts[c_idx]++;
            if(sorted==1){
                List_insert_sorted(bams->reads, rp_rev, (List_compare)bam_access_compare_read_pos_t);
            }else{
                List_push(bams->reads,rp_rev);
            }       
        } else if (strand_rev_missing) {
            read_pos_t *rp_fwd = malloc(sizeof(struct read_pos_t));
            check_mem(rp_fwd);
            bam_access_copy_read_pos(kh_val(h_strd,k_fwd), rp_fwd);
            kh_val(h_strd,k_fwd) = NULL;
            uint8_t c_idx = seq_nt16_int[rp_fwd->called_base];
            local_counts[c_idx]++;
            if(sorted==1){
                List_insert_sorted(bams->reads, rp_fwd, (List_compare)bam_access_compare_read_pos_t);
            }else{
                List_push(bams->reads, rp_fwd);
            }
        } else { //If both strands exist at this position and the bases differ we use both
            if(kh_val(h_strd, k_fwd)->called_base != kh_val(h_strd, k_rev)->called_base){ //If the bases differ
                read_pos_t *rp_fwd = malloc(sizeof(struct read_pos_t));
                check_mem(rp_fwd);
                bam_access_copy_read_pos(kh_val(h_strd, k_fwd), rp_fwd);
                kh_val(h_strd,k_fwd) = NULL;
                uint8_t c_idx_zero = seq_nt16_int[rp_fwd->called_base];
                read_pos_t *rp_rev = malloc(sizeof(struct read_pos_t));
                check_mem(rp_rev);
                bam_access_copy_read_pos(kh_val(h_strd, k_rev), rp_rev);
                kh_val(h_strd,k_rev) = NULL;
                uint8_t c_idx_one = seq_nt16_int[rp_rev->called_base];

                if(sorted==1){
                    List_insert_sorted(bams->reads, rp_fwd, (List_compare)bam_access_compare_read_pos_t);
                    List_insert_sorted(bams->reads, rp_rev, (List_compare)bam_access_compare_read_pos_t);
                }else{
                    List_push(bams->reads,rp_fwd);
                    List_push(bams->reads,rp_rev);
                }
                local_counts[c_idx_zero]++;
                local_counts[c_idx_one+4]++;
            }

        } // End of if there are missing strands
    }
    ); // End of iteration through each readname

    // Second iteration, we use only those readnames with both strands and the same base
    kh_foreach_value(h_rd_nom, h_strd, {
        int strand_fwd_missing = 0;
        int strand_rev_missing = 0;
        khiter_t k_fwd;
        khiter_t k_rev;

        k_fwd = kh_get(rpos, h_strd, 0);
        strand_fwd_missing = (k_fwd == kh_end(h_strd));
        k_rev = kh_get(rpos, h_strd, 1);
        strand_rev_missing = (k_rev == kh_end(h_strd));

        if(strand_rev_missing == 0 && strand_fwd_missing == 0){

            if(kh_val(h_strd, k_rev) && kh_val(h_strd, k_fwd)){
                
                read_pos_t *rp_fwd = malloc(sizeof(struct read_pos_t));
                check_mem(rp_fwd);
                bam_access_copy_read_pos(kh_val(h_strd,k_fwd), rp_fwd);
                kh_val(h_strd,k_fwd) = NULL;
                uint8_t c_idx_fwd = seq_nt16_int[rp_fwd->called_base];

                read_pos_t *rp_rev = malloc(sizeof(struct read_pos_t));
                check_mem(rp_rev);
                bam_access_copy_read_pos(kh_val(h_strd,k_rev), rp_rev);
                kh_val(h_strd,k_rev) = NULL;
                uint8_t c_idx_rev = seq_nt16_int[rp_rev->called_base];
                if(c_idx_rev == c_idx_fwd){
                    c_idx_rev = c_idx_rev + 4;
                    if(local_counts[c_idx_rev] < local_counts[c_idx_fwd]){
                        if(sorted==1){
                            List_insert_sorted(bams->reads, rp_rev, (List_compare)bam_access_compare_read_pos_t);
                        }else{
                            List_push(bams->reads,rp_rev);
                        }
                        local_counts[c_idx_rev]++;
                    }else{
                        if(sorted==1){
                            List_insert_sorted(bams->reads, rp_fwd, (List_compare)bam_access_compare_read_pos_t);
                        }else{
                            List_push(bams->reads,rp_fwd);
                        }
                        local_counts[c_idx_fwd]++;
                    }
                }//End of ensuring the two bases match
            }
        }


    }); // End of second iteration through readnames
    free(readname);
    free(nom);
    free(local_counts);
    //Readname hash and subhashes
    if(h_rd_nom){
        for (k = 0; k < kh_end(h_rd_nom); ++k){
            if (kh_exist(h_rd_nom, k)){
                kh_destroy(rpos, kh_val(h_rd_nom, k));
            }
        }
        kh_destroy(rdnom_rp, h_rd_nom);
    }
    
    return 0;
error:    
    if(readname) free(readname);
    if(local_counts) free(local_counts);
    if(nom) free(nom);
    //Readname hash and subhashes
    if(h_rd_nom){
        for (k = 0; k < kh_end(h_rd_nom); ++k){
            if (kh_exist(h_rd_nom, k)){
                kh_destroy(rpos, kh_val(h_rd_nom, k));
            }
        }
        kh_destroy(rdnom_rp, h_rd_nom);
    }
    return 1;
}

int reads_at_pos_callback_olap_old(uint32_t tid, uint32_t pos, int n_plp, const bam_pileup1_t *pil, void *data, int sorted, int isnorm){
    file_holder* bams = (file_holder* ) data;
    int *local_counts = NULL;
    char *nom = NULL;
    khash_t(rdnom_rp) *h_rd_nom;
	khiter_t k_rd_nom_iter;
    khash_t(rpos) *h_strd;
    khiter_t k_strd_iter;
    h_rd_nom = kh_init(rdnom_rp);
    local_counts = calloc(8,sizeof(int));
    check_mem(local_counts);
    nom = malloc(sizeof(char) * 350);
    check_mem(nom);
    int i=0;
    for(i=0;i<n_plp;i++){
        const bam_pileup1_t *p = pil + i;
        int qual = bam_get_qual(p->b)[p->qpos];
        uint8_t c = bam_seqi(bam_get_seq(p->b), p->qpos);
        int rd_name_missing;
        if(!(p->is_del) && qual >= min_base_qual && seq_nt16_int[c] < 4){
            fprintf(stderr,"%s\n", bam_get_qname(p->b));
            //Build a new read pos struct and add to the map.
            read_pos_t *rp = malloc(sizeof(struct read_pos_t));
            check_mem(rp);
            rp->rd_len = bam_cigar2qlen(p->b->core.n_cigar,bam_get_cigar(p->b));
            rp->ref_pos = pos+1;
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
            int lane_i = alg_bean_get_index_for_str_arr(bams->bean->lane,nom);
            check(lane_i>=0,"Error calculating lane index %s.",nom);
            rp->lane_i = lane_i;
            rp->normal = isnorm;
            //Logic to put this into the readname map
            k_rd_nom_iter = kh_put(rdnom_rp, h_rd_nom, bam_get_qname(p->b), &rd_name_missing);
            //Check to see if we already have an entry for this readname
            if(rd_name_missing){ // If the readname hash doesn't exist
                h_strd = kh_init(rpos);
                kh_value(h_rd_nom, k_rd_nom_iter) = h_strd;
            }else{ //Otherwise fetch the strand hash
                k_rd_nom_iter = kh_get(rdnom_rp, h_rd_nom, bam_get_qname(p->b));
                h_strd = kh_val(h_rd_nom,k_rd_nom_iter);
            }

            //Check to see if we already have a strand entry (should never happen)
            int strand_missing;
            k_strd_iter = kh_put(rpos, h_strd, rp->strand, &strand_missing);
            if(!strand_missing){
                sentinel("Found strand already present for this read");
            }
            kh_value(h_strd, k_strd_iter) = rp;
        } //If this is a usable base
    } // End of iteration through each pileup object
    //Iterate through each readname hashkey and only use the ones hitting a single strand
    kh_foreach_value(h_rd_nom, h_strd, {
        int strand_zero_is_missing = 0;
        int strand_one_is_missing = 0;

        k_strd_iter = kh_get(rpos, h_strd, 0);
        strand_zero_is_missing = (k_strd_iter == kh_end(h_strd));
        k_strd_iter = kh_get(rpos, h_strd, 1);
        strand_one_is_missing = (k_strd_iter == kh_end(h_strd));

        if(strand_zero_is_missing){ //If we've only covered one or other of the strands...
            //Get strand one rp object
            k_rd_nom_iter = kh_get(strd, h_strd, 1);
            read_pos_t *rp_one = kh_val(h_strd,k_rd_nom_iter);       
            uint8_t c_idx = seq_nt16_int[rp_one->called_base];
            c_idx = c_idx+4;
            local_counts[c_idx]++;
            if(sorted==1){
                List_insert_sorted(bams->reads, rp_one, (List_compare)bam_access_compare_read_pos_t);
            }else{
                List_push(bams->reads,rp_one);
            }
        } else if(strand_one_is_missing) {
            //Get strand zero rp object
            k_rd_nom_iter = kh_get(strd, h_strd, 0);
            read_pos_t *rp_zero = kh_val(h_strd,k_rd_nom_iter);
            uint8_t c_idx = seq_nt16_int[rp_zero->called_base];
            local_counts[c_idx]++;

            if(sorted==1){
                List_insert_sorted(bams->reads, rp_zero, (List_compare)bam_access_compare_read_pos_t);
            }else{
                List_push(bams->reads,rp_zero);
            }
        } else {// Only use both strands if the bases are different in this section
            k_rd_nom_iter = kh_get(strd, h_strd, 1);
            read_pos_t *rp_one = kh_val(h_strd,k_rd_nom_iter);       
            uint8_t c_idx_one = seq_nt16_int[rp_one->called_base];
            c_idx_one = c_idx_one+4;

            k_rd_nom_iter = kh_get(strd, h_strd, 0);
            read_pos_t *rp_zero = kh_val(h_strd,k_rd_nom_iter);
            uint8_t c_idx_zero = seq_nt16_int[rp_zero->called_base];

            if(rp_zero->called_base != rp_one->called_base){
                if(sorted==1){
                    List_insert_sorted(bams->reads, rp_zero, (List_compare)bam_access_compare_read_pos_t);
                    List_insert_sorted(bams->reads, rp_one, (List_compare)bam_access_compare_read_pos_t);
                }else{
                    List_push(bams->reads,rp_zero);
                    List_push(bams->reads,rp_one);
                }
                local_counts[c_idx_zero]++;
                local_counts[c_idx_one]++;
            }
        }
    });


    //Second iteration, only using those entries where we have both strands, in combination
    // with the counts we have built up.
    kh_foreach_value(h_rd_nom, h_strd, {
        int strand_zero_is_missing = 0;
        int strand_one_is_missing = 0;

        khiter_t k_strd_iter_zero = kh_get(strd, h_strd, 0);
        strand_zero_is_missing = (k_strd_iter_zero == kh_end(h_strd));
        khiter_t k_strd_iter_one = kh_get(strd, h_strd, 1);
        strand_one_is_missing = (k_strd_iter_one == kh_end(h_strd));
        if(strand_zero_is_missing == 0 && strand_one_is_missing == 0){ //If we've covered both strands...
            k_strd_iter_zero = kh_get(strd, h_strd, 1);
            read_pos_t *rp_one = kh_val(h_strd,k_strd_iter_zero);   
            uint8_t c_idx_one = seq_nt16_int[rp_one->called_base];
            c_idx_one = c_idx_one+4;
            k_strd_iter_zero = kh_get(strd, h_strd, 0);
            read_pos_t *rp_zero = kh_val(h_strd,k_strd_iter_zero);
            uint8_t c_idx_zero = seq_nt16_int[rp_zero->called_base];
            //If the bases are the same
            if(rp_zero->called_base == rp_one->called_base){
                if(local_counts[c_idx_one] < local_counts[c_idx_zero]){
                    if(sorted==1){
                        List_insert_sorted(bams->reads, rp_one, (List_compare)bam_access_compare_read_pos_t);
                    }else{
                        List_push(bams->reads,rp_one);
                    }
                    local_counts[c_idx_one]++;
                } else {
                    if(sorted==1){
                        List_insert_sorted(bams->reads, rp_zero, (List_compare)bam_access_compare_read_pos_t);
                    }else{
                        List_push(bams->reads,rp_zero);
                    }
                    local_counts[c_idx_zero]++;
                }
            }
        }
    });
    fprintf(stderr, "TEST8\n");
    free(nom);
    free(local_counts);
    fprintf(stderr, "TEST9\n");
    kh_destroy(rpos, h_strd);
    fprintf(stderr, "TEST10\n");
    kh_destroy(rdnom_rp, h_rd_nom);
    fprintf(stderr, "TEST11\n");
    return 0;
error:
    if(nom) free(nom);
    if(local_counts) free(local_counts);
	kh_destroy(rpos, h_strd);
    kh_destroy(rdnom_rp, h_rd_nom);
    return 1;
}

List *bam_access_get_sorted_reads_at_this_pos(char *chr_name, uint32_t start, uint32_t stop, uint8_t sorted, alg_bean_t *bean, file_holder* bams, uint8_t normal){
	//Pileup and populate the list with valid reads.
	char *region = NULL;
	bam_plp_t buf = NULL;
	bam1_t *b = NULL;
	hts_itr_t *iter = NULL;

	char sta[20];
	region = malloc(sizeof(chr_name)+sizeof(":")+sizeof("-")+(sizeof(sta)*2));
	sprintf(region,"%s:%lu-%lu",chr_name,(long unsigned int)start,(long unsigned int)stop);
	bams->beg = start;
	bams->end = stop;
	bams->reads = List_create();
	bams->bean = bean;
	isnorm = normal;
	sprintf(norm_char,"%i",(int)normal);
  // initialize pileup
  buf = bam_plp_init(pileup_blank,(void *)bams);
  bam_plp_set_maxcnt(buf,maxitercnt);
  int count = 0;
  b = bam_init1();
  iter = sam_itr_querys(bams->idx, bams->head, region);
  int result;
  int tid, pos, n_plp = -1;
  const bam_pileup1_t *pil;
  while ((result = sam_itr_next(bams->in, iter, b)) >= 0) {
    if(b->core.qual == 0
        || (b->core.flag & BAM_FUNMAP)
        || (b->core.flag & BAM_FSECONDARY)
        || (b->core.flag & BAM_FQCFAIL)
        || (b->core.flag & BAM_FSUPPLEMENTARY)){
      continue;
		}
		if((include_dup == 0 && (b->core.flag & BAM_FDUP))){
		  continue;
	  }
    if(!(include_se == 0 
          && (b->core.flag & BAM_FPROPER_PAIR) 
          && !(b->core.flag & BAM_FMUNMAP))
      ){
      continue;
    }
    if(include_se == 0){ // If we're looking for proper pairs only
      if(! (b->core.flag & BAM_FPROPER_PAIR)) continue;
      if(!(b->core.flag & BAM_FREVERSE) == !(b->core.flag & BAM_FMREVERSE)) continue;
    }
    if(include_sw == 0){
      uint8_t *xt_data = bam_aux_get(b,"XT");
      if(xt_data != NULL && bam_aux2A(xt_data) == 'M'){
        continue;
      }
    }
	  count++;
    bam_plp_push(buf, b);
    while ( (pil=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
      if(!((pos+1) >= norm->beg &&
                    (pos+1) <= norm->end)){
          continue;
      }
      int res = reads_at_pos_callback(tid, pos, n_plp, pil, bams, sorted, isnorm);
      check(res==0,"Error running callback");
    }

  }//End of while iterator to populate pileup

  sam_itr_destroy(iter);
	bam_plp_push(buf,0); // finalize pileup


  while ( (pil=bam_plp_next(buf, &tid, &pos, &n_plp)) > 0) {
    if(!((pos+1) >= norm->beg 
        && (pos+1) <= norm->end)) continue;

    int res = reads_at_pos_callback(tid, pos, n_plp, pil, bams, sorted, isnorm);
    check(res==0,"Error running callback");
  }//End of iterating through pileups
	bam_destroy1(b);
  bam_plp_destroy(buf);
	free(region);
	return bams->reads;
error:
  if(iter) sam_itr_destroy(iter);
	if(b) bam_destroy1(b);
	if(buf) bam_plp_destroy(buf);
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
	bam1_t* b = NULL;
	hts_itr_t *iter = NULL;
	char *region;
	char sta[20];
	region = malloc(sizeof(chr_name)+sizeof(":")+sizeof("-")+(sizeof(sta)*2));
	sprintf(region,"%s:%d-%d",chr_name,start,stop);
	fprintf(stderr,"REGION: %s\n",region);
  b = bam_init1();
  iter = sam_itr_querys(fh->idx, fh->head, region);
  int result;
  while ((result = sam_itr_next(fh->in, iter, b)) >= 0){

      //check Mapping Quality and not un mapped //4 // read unmapped
    if(b->core.qual == 0 || (b->core.flag & BAM_FUNMAP)
        || (b->core.flag & BAM_FSECONDARY)
        || (b->core.flag & BAM_FQCFAIL)
        || (b->core.flag & BAM_FSUPPLEMENTARY)){
      continue;
    }
    //1024 is PCR/optical duplicate
    if((include_dup == 0 && (b->core.flag & BAM_FDUP))){
      continue;
    }
    //Proper pair and mate unmapped
    if(!(include_se == 0 && (b->core.flag & BAM_FPROPER_PAIR) && !(b->core.flag & BAM_FMUNMAP))){
      continue;
    }
    if(include_se == 0){ // If we're looking for proper pairs only
      if(! (b->core.flag & BAM_FPROPER_PAIR)) continue;
      if(!(b->core.flag & BAM_FREVERSE) == !(b->core.flag & BAM_FMREVERSE)) continue;
    }
    //printf("XT DATA: %c\n",xt);
    //Now we check aux data for XT:M flags (the SW mapped marker from BWA)
    if(include_sw == 0){
      uint8_t *xt_data = bam_aux_get(b,"XT");
      if(xt_data != NULL && bam_aux2A(xt_data) == 'M'){
        continue;
      }
    }
    counter++;

  }//End of iteration through reads in this region
  sam_itr_destroy(iter);
  bam_destroy1(b);
  free(region);
  return counter;
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
	char ** ptr = NULL;
	char *tmp_line = NULL;
	htsFile *bam =  NULL;
	int rg_found = 0;
	bam = hts_open(bam_loc, "r");
	check(bam != 0,"Bam file %s failed to open to read header.",bam_loc);
	bam_hdr_t *header = sam_hdr_read(bam);
	char *head_txt = header->text;
	li = List_create();
	line = strtok(head_txt,"\n");
	while(line != NULL){
		//Check for a read group line
		if(strncmp(line,"@RG",3)==0){
			ptr = malloc(sizeof(char **));
			check_mem(ptr);
			char *id = malloc(sizeof(char) * 100);
			char *lane = malloc(sizeof(char) * 150);
			tmp_line = strtok_r(line,"\t",ptr);
			while(tmp_line != NULL){
				int chk = sscanf(tmp_line,"ID:%s",id);
				if(chk==1){
					rg_found = 1;
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

				}// If this is a match for RG
				tmp_line = strtok_r(NULL,"\t",ptr);
			}
		}//End of if this is an RG line
		line = strtok(NULL,"\n");
	}

	check(rg_found==1,"No RG lines with IDs found in header of bam file %s.",bam_loc);

	if(line) free(line);
	hts_close(bam);
	return li;
error:
	if(line) free(line);
	if(bam) hts_close(bam);
	if(li) List_clear_destroy(li);
	return NULL;
}

htsFile *bam_access_populate_file(const char *bam_loc, const char *ref_file){
	htsFile *sf = NULL;
	sf = hts_open(bam_loc,"r");
	check(sf != 0,"File %s failed to open.",bam_loc);
	hts_set_fai_filename(sf, ref_file);
	return sf;
error:
	if(sf) hts_close(sf);
	return NULL;
}

hts_idx_t *bam_access_populate_file_index(samFile *sf, const char *bam_loc){
	hts_idx_t *idx = NULL;
	idx = sam_index_load(sf,bam_loc);
	check(idx != NULL,"Index %s failed to open.",bam_loc);
	return idx;
error:
	if(idx) hts_idx_destroy(idx);
	return NULL;
}

hts_itr_t *bam_access_get_hts_itr(htsFile *sf, hts_idx_t *idx, const char *chr, uint32_t from, uint32_t to){
	hts_itr_t *iter = NULL;
	bam_hdr_t *head = NULL;
	head = sam_hdr_read(sf);
	int tid = bam_name2id(head, chr);
	free(head);
	iter = sam_itr_queryi(idx, tid, from, to);
	check(iter!=0,"Error fetching iterator %s:%d-%d.",chr,from,to);
	return iter;
error:
	if(iter) hts_itr_destroy(iter);
	return NULL;
}
