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

static int includeSW = 0;
static int includeSingleEnd = 0;
static int includeDups = 0;
static unsigned int increment = 250000;
static unsigned int read_count = 1000000;
static double maxPropRdCount = 1.5;
static char tum_bam_file[512];
static char norm_bam_file[512];
static char *config_file = "caveman.cfg.ini";
static char results[512];// = "results";
static char ref_idx[512];// = "";
static char list_loc[512];// = "splitList";
static char alg_bean_loc[512];// = "alg_bean";
static char ignore_regions_file[512];// = NULL;
static int idx;


void split_print_usage (int exit_code){
	printf ("Usage: caveman split -i jobindex [-f path] [-c int] [-m int] [-e int] \n\n");
  	printf("-i  --index [int]                 Job index (e.g. from $LSB_JOBINDEX)\n\n");
	printf("Optional\n");
	printf("-f  --config-file [file]          Path to the config file produced by setup [default:'%s'].\n",config_file);
	printf("-c  --increment [int]             Increment to use when deciding split sizes [default:%d]\n",increment);
	printf("-m  --max-read-count [double]     Proportion of read-count to allow as a max in a split section [default:%f]\n",maxPropRdCount);
	printf("-e  --read-count [int]            Guide for maximum read count in a section [default:%d]\n",read_count);
	printf("-h	help                          Display this usage information.\n");
  exit(exit_code);
}

void split_setup_options(int argc, char *argv[]){
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
   while((iarg = getopt_long(argc, argv, "f:i:m:c:e:h",
                            								long_opts, &index)) != -1){
   	switch(iarg){
   		case 'h':
         	split_print_usage(0);
         	break;

      	case 'f':
      		config_file = optarg;
      		break;

      	case 'i':
      		idx = atoi(optarg);
      		break;

      	case 'm':
      		maxPropRdCount = atof(optarg);
      		break;

      	case 'c':
      		increment = atoi(optarg);
      		break;

      	case 'e':
      		read_count = atoi(optarg);
      		break;

			case '?':
            split_print_usage (1);
            break;

      	default:
      		split_print_usage (1);

   	}; // End of args switch statement

   }//End of iteration through options

   //Do some checking to ensure required arguments were passed
   if(idx == NULL || idx == 0){
   	split_print_usage(1);
   }

   if(check_exist(config_file) != 1){
   	printf("Config file %s does not appear to exist. Have you run caveman setup?\n",config_file);
   	split_print_usage(1);
   }

   return;
}

int split_main(int argc, char *argv[]){
	split_setup_options(argc,argv);

	//Open the config file and do relevant things
	FILE *config = fopen(config_file,"r");
	check(config != NULL,"Failed to open config file for reading. Have you run caveman-setup?");

	int cfg = config_file_access_read_config_file(config,tum_bam_file,norm_bam_file,ref_idx,ignore_regions_file,alg_bean_loc,
								results,list_loc,&includeSW,&includeSingleEnd,&includeDups);

	check(cfg==0,"Error parsing config file.");
   bam_access_include_sw(includeSW);
   bam_access_include_se(includeSingleEnd);
   bam_access_include_dup(includeDups);

   //Open reference file and read in chromosomes - getting chr name and length for this index
   int chr_length = 0;
   char *chr_name = malloc(sizeof(char *));
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
   int ignore_reg_count = ignore_reg_access_get_ign_reg_count_for_chr(ignore_regions_file,chr_name);
   check(ignore_reg_count >= 0,"Error trying to check the number of ignored regions for this chromosome.");

   printf("Found %d ignored regions for chromosome %s.\n",ignore_reg_count,chr_name);

   //Now create a store for said regions.
   struct seq_region_t **ignore_regs;
   ignore_regs = malloc(sizeof(struct seq_region_t *) *  ignore_reg_count);
   check_mem(ignore_regs);
   check(ignore_reg_access_get_ign_reg_for_chr(ignore_regions_file,chr_name,ignore_reg_count,ignore_regs)==0,"Error fetching ignored regions from file.");

   //Check there's not a whole chromosome block.
   if(!(ignore_reg_count == 1 && ignore_regs[0]->beg == 1 && ignore_regs[0]->end >= chr_length)){
		//No chromosome block, so carry on.
		//Open bam file and iterate through chunks until we reach the cutoff.
		chk = bam_access_openbams(norm_bam_file,tum_bam_file);
		check(chk == 0,"Error trying to open bam files.");

		int sect_start = 1;
		int sect_stop = 1;
		int rdCount = 0;
		int last_stop = sect_start-1;
		//Iterate through sections, checking for overlap with ignored regions until we reach the cutoff.
		while(sect_start<chr_length){
			List *ign_this_sect = List_create();
			if(sect_stop == 1){
				sect_stop = (sect_start + increment);
			}
			if(sect_stop > chr_length){
				sect_stop = chr_length;
			}
			//Check if stop is in an ignored region
			List *contained = ignore_reg_access_get_ign_reg_contained(sect_start,sect_stop,ignore_regs,ignore_reg_count);
			if(List_count(contained) > 0){
				LIST_FOREACH(contained, first, next, cur){
					List_push(ign_this_sect,(seq_region_t *)cur->value);
				}
			}
			seq_region_t *overlap = ignore_reg_access_get_ign_reg_overlap(sect_stop,ignore_regs,ignore_reg_count);
			if(overlap != NULL){
				sect_stop = overlap->end;
				List_push(ign_this_sect,overlap);
			}
			rdCount += get_read_counts_with_ignore(ign_this_sect,last_stop+1,sect_stop,chr_name);
			List_destroy(contained);
			check(rdCount >= 0,"Problem retrieving reads for section.");
			if(rdCount < read_count){
				if(sect_stop >= chr_length){
					split_access_print_section(output,chr_name,sect_start,sect_stop);
					sect_start = sect_stop+1;
					last_stop = sect_stop;
					sect_stop = sect_start + increment;
					List_clear_destroy(ign_this_sect);
					rdCount = 0;
					continue;
				}else{
					last_stop = sect_stop;
					sect_stop += increment;
					List_clear_destroy(ign_this_sect);
				}
			}else if(rdCount > (read_count * maxPropRdCount)){
				sect_stop = shrink_section_to_size(chr_name,sect_start,sect_stop,ignore_regs,ignore_reg_count);
				check(sect_stop > 0,"Error resizing over sized section.");
				split_access_print_section(output,chr_name,sect_start,sect_stop);
				sect_start = sect_stop+1;
				last_stop = sect_stop;
				sect_stop = sect_start + increment;
				List_clear_destroy(ign_this_sect);
				rdCount = 0;
				continue;
			}else if(rdCount >= read_count && rdCount <= (read_count * maxPropRdCount)){
				split_access_print_section(output,chr_name,sect_start,sect_stop);
				sect_start = sect_stop+1;
				last_stop = sect_stop;
				sect_stop = sect_start + increment;
				List_clear_destroy(ign_this_sect);
				rdCount = 0;
				continue;
			}else{
				sentinel("We shouldn't have reached this.");
			}
		}
		//Close bams
		bam_access_closebams();
   }

	ignore_reg_access_destroy_seq_region_t_arr(ignore_reg_count, ignore_regs);
	free(chr_name);
	//Close output file
	fclose(output);
  return 0;
error:
	if(chr_name) free(chr_name);
	ignore_reg_access_destroy_seq_region_t_arr(ignore_reg_count, ignore_regs);
	bam_access_closebams();
	return -1;
}

int shrink_section_to_size(char *chr_name,int sect_start, int sect_stop, struct seq_region_t **ignore_regs, int ignore_reg_count){
	int read_num = (read_count * maxPropRdCount) + 1;
	int new_inc = increment;
	new_inc /= 16;
	while(read_num > (int)(read_count * maxPropRdCount)){
		List *ign_this_sect = List_create();
		sect_stop -= new_inc;
		//Check we haven't gone back past the start of the last RG, if we have, try again!
		if(sect_stop <= sect_start){
			sect_stop = sect_start + new_inc;
			new_inc = new_inc / 2;
			sect_stop -= new_inc;
		}
		//Check if stop is in an ignored region
		List *contained = ignore_reg_access_get_ign_reg_contained(sect_start,sect_stop,ignore_regs,ignore_reg_count);
		check(contained != NULL, "Problem fetching contained regions for %d-%d.",sect_start,sect_stop);
		if(List_count(contained) > 0){
			LIST_FOREACH(contained, first, next, cur){
				List_push(ign_this_sect,(seq_region_t *)cur->value);
			}
		}
		seq_region_t *overlap = ignore_reg_access_get_ign_reg_overlap(sect_stop,ignore_regs,ignore_reg_count);
		if(overlap != NULL){
			sect_stop = overlap->end;
			List_push(ign_this_sect,overlap);
		}
		read_num = get_read_counts_with_ignore(ign_this_sect,sect_start,sect_stop,chr_name);
		free(contained);
		check(read_num >= 0,"Problem retrieving reads for section.");
		List_clear_destroy(ign_this_sect);
	}
	return sect_stop;
error:
	return -1;
}

int get_read_counts_with_ignore(List *ignore,int start, int stop, char *chr){
	int total = 0;
	if(List_count(ignore) > 0){
		int tmp_start = start;
		int tmp_stop = stop;
		LIST_FOREACH(ignore,first,next,cur){
			if(!(((seq_region_t *)cur->value)->beg >= start && ((seq_region_t *)cur->value)->end <= stop)) continue;
			int sta = ((seq_region_t *)cur->value)->beg;
			int sto = ((seq_region_t *)cur->value)->end;
			tmp_stop = sta;
			total+= bam_access_get_count_for_region(chr,tmp_start,tmp_stop - 1);
			tmp_start = sto + 1;
			tmp_stop = stop;
		}
		total+= bam_access_get_count_for_region(chr,tmp_start,tmp_stop);
	}else{
		total = bam_access_get_count_for_region(chr,start,stop);
	}
	return total;
}
