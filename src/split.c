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
#include <getopt.h>
#include <string.h>
#include <assert.h>
#include <dbg.h>
#include <split.h>
#include <file_tests.h>
#include <alg_bean.h>
#include <fai_access.h>
#include <split_access.h>
#include <bam_access.h>
#include <config_file_access.h>
#include <cn_access.h>
#include "khash.h"

//New hash to store unique readlengths
KHASH_MAP_INIT_INT(rdlenkhash, uint8_t)

static int includeSW = 0;
static int includeSingleEnd = 0;
static int includeDups = 0;
static unsigned int max_read_count = 350000;
static unsigned int read_length_base = 100;
static char tum_bam_file[512];
static char norm_bam_file[512];
static char *config_file = "caveman.cfg.ini";
static char results[512];// = "results";
static char ref_idx[512];// = "";
static char list_loc[512];// = "splitList";
static char alg_bean_loc[512];// = "alg_bean";
static char version[50];// = "alg_bean";
static char ignore_regions_file[512];// = NULL;
static char norm_cn_loc[512];
static char tum_cn_loc[512];
static int idx = 0;

void split_print_usage (int exit_code){
    printf ("Usage: caveman split -i jobindex [-f path] [-c int] [-m int] [-e int] \n\n");
    printf("-i  --index [int]                 Job index (e.g. from $LSB_JOBINDEX)\n\n");
    printf("Optional\n");
    printf("-f  --config-file [file]          Path to the config file produced by setup [default:'%s'].\n",config_file);
    printf("-e  --read-count [int]            Guide for maximum read count in a section [default:%d]\n",max_read_count);
    printf("-h    help                          Display this usage information.\n");
  exit(exit_code);
}

int split_setup_options(int argc, char *argv[]){
    const struct option long_opts[] =
    {
        {"config-file", required_argument, 0, 'f'},
        {"increment", required_argument, 0, 'c'},
        {"max-read-count",required_argument , 0, 'm'},
        {"read-count", required_argument, 0, 'e'},
        {"index", required_argument, 0, 'i'},
        {"help", no_argument, 0, 'h'},
        { NULL, 0, NULL, 0}
   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

    //Iterate through options
    while((iarg = getopt_long(argc, argv, "f:i:e:h",long_opts, &index)) != -1){
        switch(iarg){
            case 'h':
                split_print_usage(0);
                break;

            case 'f':
                config_file = optarg;
                break;

            case 'i':
                if(sscanf(optarg, "%i", &idx) != 1){
                    sentinel("Error parsing -i argument '%s'. Should be an integer > 0",optarg);
                }
                break;

            case 'e':
                if(sscanf(optarg, "%i", &max_read_count) != 1){
                    sentinel("Error parsing -e argument '%s'. Should be an integer > 0",optarg);
                }
                break;

            case '?':
                split_print_usage (1);
                break;

            default:
                split_print_usage (1);

       }; // End of args switch statement

    }//End of iteration through options

    //Do some checking to ensure required arguments were passed
    if(idx == 0){
        split_print_usage(1);
    }

    if(check_exist(config_file) != 1){
        printf("Config file %s does not appear to exist. Have you run caveman setup?\n",config_file);
        split_print_usage(1);
    }

    return 0;

error:
    split_print_usage(1);
    return -1;
}

int round_divide_integer(int dividend, int divisor){
    if(dividend == 0 || divisor == 0){
        return 1;
    }
    return (dividend + (divisor / 2)) / divisor;
}

uint32_t min(uint32_t one, uint32_t two){
    if(one<=two){
        return one;
    }else if(two<one){
        return two;
    }
    return one;
}

int split_main(int argc, char *argv[]){

    htsFile *sf_norm = NULL;
    hts_idx_t *idx_norm = NULL;
    htsFile *sf_tum = NULL;
    hts_idx_t *idx_tum = NULL;
    hts_itr_t *iter_norm = NULL;
    hts_itr_t *iter_tum = NULL;
    bam1_t *norm_read = NULL;
    bam1_t *tum_read = NULL;
    seq_region_t **ignore_regs = NULL;
    FILE *alg_bean_file = NULL;
    alg_bean_t *alg = NULL;
    char *read_len_pos_arr_file = NULL;
    khash_t(rdlenkhash) *rdlen_h;
    int ignore_reg_count = 0;
    FILE *output_rp = NULL;

    int is_err = split_setup_options(argc,argv);
    check(is_err==0,"Error parsing options");

    char *chr_name = malloc(sizeof(char *));
    //Open the config file and do relevant things
    FILE *config = fopen(config_file,"r");
    check(config != NULL,"Failed to open config file for reading. Have you run caveman-setup?");

    int cfg = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
                                results,list_loc,&includeSW,&includeSingleEnd,&includeDups,version,norm_cn_loc,tum_cn_loc);


    check(strcmp(version,CAVEMAN_VERSION)==0,"Stored version in %s %s and current code version %s did not match.",config_file,version,CAVEMAN_VERSION);

    check(cfg==0,"Error parsing config file.");
    bam_access_include_sw(includeSW);
    bam_access_include_se(includeSingleEnd);
    bam_access_include_dup(includeDups);

    //Open reference file and read in chromosomes - getting chr name and length for this index
    int chr_length = 0;
    int chk = 0;
    chk = fai_access_get_name_from_index(idx, ref_idx, chr_name, &chr_length);
    check(chk==0, "Error encountered trying to get chromosome name and length from FASTA index file.");

    printf("Found chr: %s of length: %d at index %d\n",chr_name,chr_length,idx);

     //Open a file to write sections, named according to CHR.
    char *fname = malloc(strlen(chr_name) + strlen(list_loc) + 3);
    check_mem(fname);
    //Create filename here through name concatenation.
    strcpy(fname,list_loc);
    strcat(fname,".");
    strcat(fname,chr_name);
    FILE *output = fopen(fname,"w");
    free(fname);
    check(output != NULL, "Error opening file %s for write.",fname);

    //Load in a set of ignore regions from tsv format, only require this chromosome.
    ignore_reg_count = ignore_reg_access_get_ign_reg_count_for_chr(ignore_regions_file,chr_name);
    check(ignore_reg_count >= 0,"Error trying to check the number of ignored regions for this chromosome.");

    printf("Found %d ignored regions for chromosome %s.\n",ignore_reg_count,chr_name);

    //Now create a store for said regions.
    ignore_regs = malloc(sizeof(struct seq_region_t *) *  ignore_reg_count);
    check_mem(ignore_regs);
    check(ignore_reg_access_get_ign_reg_for_chr(ignore_regions_file,chr_name,ignore_reg_count,ignore_regs)==0,"Error fetching ignored regions from file.");

    //Check there's not a whole chromosome block.
    if(!(ignore_reg_count == 1 && ignore_regs[0]->beg == 1 && ignore_regs[0]->end >= chr_length)){
        //Initialise readlength hash
        rdlen_h = kh_init(rdlenkhash);
        check(rdlen_h != NULL,"Memory allocation error when initialising hash kh_init(rdlenkhash).");

        //No chromosome block, so carry on.
        uint32_t start = 1;
        uint32_t stop = chr_length;
        uint64_t rd_count = 0;

        sf_norm = bam_access_populate_file(norm_bam_file,ref_idx);
        check(sf_norm!=NULL,"Error populating file norm seq file %s.",norm_bam_file);
        idx_norm = bam_access_populate_file_index(sf_norm, norm_bam_file);
        check(idx_norm!=NULL,"Error populating index for norm seq file %s.",norm_bam_file);
        sf_tum = bam_access_populate_file(tum_bam_file,ref_idx);
        check(sf_tum!=NULL,"Error populating file for tum seq file %s.",tum_bam_file);
        idx_tum = bam_access_populate_file_index(sf_tum, tum_bam_file);
        check(idx_tum!=NULL,"Error populating index for tum seq file %s.",tum_bam_file);

        //read the first 100 reads and get an idea of average read length.
        int avg_read_len_norm = bam_access_get_avg_readlength_from_bam(sf_norm);
        int avg_read_len_tum = bam_access_get_avg_readlength_from_bam(sf_tum);
        //Use a comparison of average read length to read_length_base in order to calculate a useful split size.
        float avg_read_len = ((float)avg_read_len_norm + (float)avg_read_len_tum) / (float)2;
        //Adjust max read count according to difference between avg_read_len and read_length_base
        float proportion_rd_length =  (float)read_length_base / avg_read_len;
        max_read_count = (int)((float)max_read_count * proportion_rd_length);

        hts_close(sf_norm);
        hts_idx_destroy(idx_norm);
        hts_close(sf_tum);
        hts_idx_destroy(idx_tum);

        sf_norm = bam_access_populate_file(norm_bam_file,ref_idx);
        check(sf_norm!=NULL,"Error populating file norm seq file %s.",norm_bam_file);
        idx_norm = bam_access_populate_file_index(sf_norm, norm_bam_file);
        check(idx_norm!=NULL,"Error populating index for norm seq file %s.",norm_bam_file);
        sf_tum = bam_access_populate_file(tum_bam_file,ref_idx);
        check(sf_tum!=NULL,"Error populating file for tum seq file %s.",tum_bam_file);
        idx_tum = bam_access_populate_file_index(sf_tum, tum_bam_file);
        check(idx_tum!=NULL,"Error populating index for tum seq file %s.",tum_bam_file);

        iter_norm = bam_access_get_hts_itr(sf_norm, idx_norm, chr_name, start, stop);
        check(iter_norm!=NULL,"Error fetching normal iterator or section %s:%d-%d.",chr_name,start,stop);
        iter_tum = bam_access_get_hts_itr(sf_tum, idx_tum, chr_name, start, stop);
        check(iter_tum!=NULL,"Error fetching tumour iterator or section %s:%d-%d.",chr_name,start,stop);

        //Setup a read for iteration
        norm_read = bam_init1();
        tum_read = bam_init1();
        int iter_n_status = 0;
        int iter_t_status = 0;
        uint32_t sect_start = 1;
        uint32_t sect_stop = 0;
        uint32_t curr_n_pos = 0;
        uint32_t curr_t_pos = 0;

        khiter_t k;
        //Have both iterators, now need to iterate through each in sync so we don't get ahead of the stops.
        while(iter_n_status>=0 || iter_t_status>=0){ //Keep iterating until both iterators are out of reads.
            while(curr_n_pos<=curr_t_pos && iter_n_status>=0 && iter_t_status>=0 && rd_count<=max_read_count){ //While the positions aren't equal and tumour has reads left. Normal jumps ahead
                iter_n_status = sam_itr_next(sf_norm,iter_norm,norm_read);
                check(iter_n_status>=-1,"Error detected (%d) when trying to iterate through region.",iter_n_status);
                curr_n_pos = norm_read->core.pos;
                if(iter_n_status>=0 && bam_access_check_bam_flags(norm_read) == 1 && ignore_reg_access_get_ign_reg_overlap(curr_n_pos,ignore_regs,ignore_reg_count) == NULL){
                    int read_len = norm_read->core.l_qseq;
                    int rd_len_missing;
                    k = kh_put(rdlenkhash, rdlen_h, read_len, &rd_len_missing);
                    if(rd_len_missing){ // If the readname key doesn't yet exist
                        kh_value(rdlen_h, k) = read_len;
                    }
                    rd_count++;
                }
            }//End of this iteration through normal reads

            while(curr_t_pos<=curr_n_pos && iter_t_status>=0 && iter_n_status>=0){ //While the positions aren't equal and normal has reads left
                iter_t_status = sam_itr_next(sf_tum,iter_tum,tum_read);
                check(iter_t_status>=-1,"Error detected (%d) when trying to iterate through region.",iter_t_status);
                curr_t_pos = tum_read->core.pos;
                if(iter_t_status>=0 && bam_access_check_bam_flags(tum_read) == 1 && ignore_reg_access_get_ign_reg_overlap(curr_t_pos,ignore_regs,ignore_reg_count) == NULL){
                    int read_len = tum_read->core.l_qseq;
                    int rd_len_missing;
                    k = kh_put(rdlenkhash, rdlen_h, read_len, &rd_len_missing);
                    if(rd_len_missing){ // If the readname key doesn't yet exist
                        kh_value(rdlen_h, k) = read_len;
                    }
                    rd_count++;
                }
            }//End of this iteration through tumour reads

            //An extra section for where one or the other iterator is out of reads (we still need to count for the mstep).
            if(iter_n_status<0 && iter_t_status>=0){ //No more normal reads
                while(iter_t_status>=0 && rd_count<=max_read_count){
                    iter_t_status = sam_itr_next(sf_tum,iter_tum,tum_read);
                    check(iter_t_status>=-1,"Error detected (%d) when trying to iterate through region.",iter_t_status);
                    curr_t_pos = tum_read->core.pos;
                    if(iter_t_status>=0 && bam_access_check_bam_flags(tum_read) == 1 && ignore_reg_access_get_ign_reg_overlap(curr_t_pos,ignore_regs,ignore_reg_count) == NULL){
                        int read_len = tum_read->core.l_qseq;
                        int rd_len_missing;
                        k = kh_put(rdlenkhash, rdlen_h, read_len, &rd_len_missing);
                        if(rd_len_missing){ // If the readname key doesn't yet exist
                            kh_value(rdlen_h, k) = read_len;
                        }
                        rd_count++;
                    }
                }
            }//End of iteration through tumour reads where only tumour reads remain

            if(iter_t_status<0 && iter_n_status>=0){ //No more tumour reads
                while(iter_n_status>=0 && rd_count<=max_read_count){
                    iter_n_status = sam_itr_next(sf_norm,iter_norm,norm_read);
                    check(iter_n_status>=-1,"Error detected (%d) when trying to iterate through region.",iter_n_status);
                    curr_n_pos = norm_read->core.pos;
                    if(iter_n_status>=0 && bam_access_check_bam_flags(norm_read) == 1 && ignore_reg_access_get_ign_reg_overlap(curr_n_pos,ignore_regs,ignore_reg_count) == NULL){
                        int read_len = norm_read->core.l_qseq;
                        int rd_len_missing;
                        k = kh_put(rdlenkhash, rdlen_h, read_len, &rd_len_missing);
                        if(rd_len_missing){ // If the readname key doesn't yet exist
                            kh_value(rdlen_h, k) = read_len;
                        }
                        rd_count++;
                    }
                }
            }//End of iteration through normal reads where only normal reads remain

            //Reads have equal start positions, check the count.
            if(rd_count>=max_read_count){
                //Set old stop position (Min of curr_t_pos & curr_n_pos)
                sect_stop = min(curr_t_pos,curr_n_pos);
                seq_region_t *reg = ignore_reg_access_get_ign_reg_overlap(sect_stop+1,ignore_regs,ignore_reg_count);
                if(reg != NULL){
                    sect_stop = reg->end+1;
                }
                //This is the position on which to separate the split sections so print it.
                if(sect_stop>0 && sect_stop >= sect_start) {
                    split_access_print_section(output,chr_name,sect_start,sect_stop);
                }
                //printf("Found %d reads for %s:%d-%d\n",rd_count,chr_name,sect_start,sect_stop);
                //Reset read count
                rd_count=1;//Set as 1 due to the way the loop records read counts .
                //Set new start position
                sect_start = sect_stop+1;
            }//End of checking if we've hit our read cutoff.

        }//End of moving through both normal and tumour iterators

        //Read in alg_bean
        //Read in the alg bean
        alg_bean_file = fopen(alg_bean_loc,"r");
        check(alg_bean_file != 0 ,"Error trying to open alg_bean file: %s.",alg_bean_loc);
        alg = alg_bean_read_file(alg_bean_file);
        check(alg != NULL,"Error reading alg_bean from file.");
        check(fclose(alg_bean_file)==0,"Error closing alg bean file.");

        //Iterate through readlen hashkeys and output a format giving all readlength options for readpos boundaries
        //This can be read in by mstep and estep

        //Open a file to write sections, named according to CHR.
        char *read_len_pos_arr_file = malloc(strlen(chr_name) + strlen(list_loc) + 6);
        check_mem(read_len_pos_arr_file);
        //Create filename here through name concatenation.
        char *dupdir;
        dupdir = strdup(list_loc);
        strcpy(read_len_pos_arr_file,dirname(dupdir));
        strcat(read_len_pos_arr_file,"/readpos.");
        strcat(read_len_pos_arr_file,chr_name);
        output_rp = fopen(read_len_pos_arr_file,"w");
        check(output_rp != NULL, "Error opening file %s for write.",read_len_pos_arr_file);
        free(read_len_pos_arr_file);

        for (k = kh_begin(rdlen_h); k != kh_end(rdlen_h); ++k){  // traverse
            // test if a bucket contains data
            if (kh_exist(rdlen_h, k)){
                //Convert alg bean read length proportions to a 'by readlength' set of ranges
                int len = kh_value(rdlen_h, k);
                List *lengths = alg_bean_get_position_list_from_read_pos_proportion_arr(alg->rd_pos,len);
                //Output readlength splits and read length itself to file
                fprintf(output_rp,"%d\t",len);
                LIST_FOREACH(lengths, first, next, cur){
                    alg_bean_intrange *range =  (alg_bean_intrange *) cur->value;
                    fprintf(output_rp,"%d-%d;",range->from,range->to);
                }
                fprintf(output_rp,"\n");
                List_clear_destroy(lengths);
            }
        }
        fclose(output_rp);
        

        //No more reads left so we must print the last section.
        kh_destroy(rdlenkhash, rdlen_h);
        split_access_print_section(output,chr_name,sect_start,chr_length);
        bam_destroy1(norm_read);
        bam_destroy1(tum_read);
        hts_close(sf_norm);
        hts_idx_destroy(idx_norm);
        hts_close(sf_tum);
        hts_idx_destroy(idx_tum);
        hts_itr_destroy(iter_norm);
        hts_itr_destroy(iter_tum);
    }//End of checking if this is a valid contig to split.


    return 0;
error:
    if(rdlen_h) kh_destroy(rdlenkhash, rdlen_h);
    if(output_rp) fclose(output_rp);
    if(norm_read) bam_destroy1(norm_read);
    if(tum_read) bam_destroy1(tum_read);
    if(sf_norm) hts_close(sf_norm);
    if(idx) hts_idx_destroy(idx_norm);
    if(sf_tum) hts_close(sf_tum);
    if(idx) hts_idx_destroy(idx_tum);
    if(iter_norm) hts_itr_destroy(iter_norm);
    if(iter_tum) hts_itr_destroy(iter_tum);
    return -1;
}
