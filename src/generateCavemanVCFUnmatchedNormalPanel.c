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
#include <dbg.h>
#include <file_tests.h>
#include <bam_access.h>
#include <split_access.h>
#include <output.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <fai_access.h>

static int min_base_qual = 0;
static char *bam_file_locs;// = NULL;
static char *split_file_loc;// = NULL;
static char *out_file_loc;// = NULL;
static char *ref_file_loc;// = NULL;
static char *species;// = NULL;
static char *species_vers;// = NULL;
static char *sample_names;// = NULL;
static char *platform = "Illumina";
static char *protocol = "WGS";
static int idx = 0;
static const char WHOLE_GENOME_PROTOCOL[] = "WGS";
static const char EXOME_PROTOCOL[] = "WGS";
static const char ILLUMINA_PLATFORM[] = "Illumina";
static const char *PERMITTED_PROTOCOLS[2] = {WHOLE_GENOME_PROTOCOL,EXOME_PROTOCOL};
static const char *PERMITTED_PLATFORMS[1] = {ILLUMINA_PLATFORM};
static const int MAX_SECTION_SIZE = 2000000;

static const char *VCF_FILE_FORMAT = "##fileformat=VCFv4.1";
static const char *VCF_FILE_DATE = "##fileDate=%s\n";
static const char *VCF_SOURCE = "##source_";
static const char *VCF_SOURCE_2 = ".1=generateCavemanVCFUnmatchedNormalPanel";
static const char *VCF_INFO = "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n";
static const char *VCF_FORMAT_GENO = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
static const char *VCF_A_ALL = "##FORMAT=<ID=AZ,Number=1,Type=Integer,Description=\"Count of A alleles at this position\">\n";
static const char *VCF_C_ALL = "##FORMAT=<ID=CZ,Number=1,Type=Integer,Description=\"Count of C alleles at this position\">\n";
static const char *VCF_G_ALL = "##FORMAT=<ID=GZ,Number=1,Type=Integer,Description=\"Count of G alleles at this position\">\n";
static const char *VCF_T_ALL = "##FORMAT=<ID=TZ,Number=1,Type=Integer,Description=\"Count of T alleles at this position\">\n";
static const char *VCF_SAMPLE = "##SAMPLE=<ID=%s,Description=\"UnmatchedNormal\",Platform=%s,Protocol=%s>\n";
static const char *VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n";
static const char *VCF_VAR_LINE_START = "%s\t%d\t.\t%c\t.\t.\t.\t";
static const char *VCF_VAR_LINE_MID = "DP=%d\tGT:AZ:CZ:GZ:TZ\t";
static const char *VCF_VAR_LINE = "0|0:%d:%d:%d:%d";

typedef struct sample_bam{
	char *sample_name;
	char *bam_file;
	struct file_holder *holder;
}sample_bam;

typedef sample_bam* sample_bam_ptr;

#define ELEMENT_TYPE sample_bam_ptr
#define ELEMENTS_PER_NODE 8
#include <List.h>
#include <List.c>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

void gen_panel_print_usage (int exit_code){
	printf("Usage: generateCavemanVCFUnmatchedNormalPanel -i jobindex [-f path] [-c int] [-m int] [-e int] \n\n");
  printf("-i  --index [int]                 Job index (e.g. from $LSB_JOBINDEX)\n\n");
  printf("-b  --bam-files [file]            Path to the bam files to be analysed (comma separated list, in the same order as sample-names).\n");
  printf("-o  --out-file [file]             Filename to write results to.\n");
  printf("-r  --reference [file]            Path to reference file (genome.fa).\n");
  printf("-l  --split-sections [file]       Path to the bed file of split sections to analyse.\n");
  printf("-s  --species [string]            Species name (used in VCF output).\n");
  printf("-v  --spp-vers [string]           Species version (used in VCF output).\n");
  printf("-n  --sample-names [string]       Name of samples being analysed (comma separated list, in the same order as bam-files).\n");
  printf("-p  --platform [string]           Platform sample was generated on (Illumina). [default:%s]\n",platform);
  printf("-t  --protocol [string]           Protocol used to generate sample data (WGS|WXS). [default:%s]\n",protocol);

	printf("Optional\n");

	printf("-q  --min-base-qual [int]         Minimum base quality [default:%d]\n",min_base_qual);
	//printf("-m  --man [int]                   Full documentation.\n");
	printf("-h	help                          Display this usage information.\n");
  exit(exit_code);
}

void gen_panel_setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
             	{"index", required_argument, 0, 'i'},
             	{"bam-file", required_argument, 0, 'b'},
             	{"out-file",required_argument , 0, 'o'},
             	{"reference", required_argument, 0, 'r'},
             	{"split-sections", required_argument, 0, 'l'},
             	{"species", required_argument, 0, 's'},
             	{"spp-vers", required_argument, 0, 'v'},
             	{"sample-names", required_argument, 0, 'n'},
             	{"platform", required_argument, 0, 'p'},
             	{"protocol", required_argument, 0, 't'},
             	{"min-base-qual", required_argument, 0, 'q'},
             	{"help", no_argument, 0, 'h'},
             	{ NULL, 0, NULL, 0}
   }; //End of declaring opts

   int index = 0;
   int iarg = 0;

   //Iterate through options
   while((iarg = getopt_long(argc, argv, "i:b:o:r:l:s:v:n:p:t:q:h",
                            								long_opts, &index)) != -1){
   	switch(iarg){
   		case 'h':
      	gen_panel_print_usage(0);
        break;

			case 'b':
				bam_file_locs = optarg;
				break;

			case 'i':
				idx = atoi(optarg);
				break;

			case 'o':
				out_file_loc = optarg;
				break;

			case 'r':
				ref_file_loc = optarg;
				break;

			case 'l':
				split_file_loc = optarg;
				break;

			case 's':
				species = optarg;
				break;

			case 'v':
				species_vers = optarg;
				break;

			case 'n':
				sample_names = optarg;
				break;

			case 'p':
				platform = optarg;
				break;

			case 't':
				protocol = optarg;
				break;

			case '?':
        gen_panel_print_usage (1);
        break;

      default:
      	gen_panel_print_usage (1);

   	}; // End of args switch statement

   }//End of iteration through options

   //Do some checking to ensure required arguments were passed
   if(idx == 0){
   	gen_panel_print_usage(1);
   }

   if(check_exist(split_file_loc) != 1){
   	printf("Split sections file %s does not appear to exist.\n",split_file_loc);
   	gen_panel_print_usage(1);
   }

   if(check_exist(ref_file_loc) != 1){
   	printf("Reference file %s does not appear to exist.\n",ref_file_loc);
   	gen_panel_print_usage(1);
   }

   int i=0;
   int found=0;
   for(i=0;i<2;i++){
   	if(strcmp(PERMITTED_PROTOCOLS[i],protocol)==0){
   		found = 1;
   	}
   }
   if(found == 0){
   	printf("Protocol %s is not contained in the valid protocols list.\n",protocol);
   }

   i=0;
   found=0;
   for(i=0;i<1;i++){
   	if(strcmp(PERMITTED_PLATFORMS[i],platform)==0){
   		found = 1;
   	}
   }
   if(found == 0){
   	printf("Platform %s is not contained in the valid platforms list.\n",platform);
   }
   return;
}

int gen_panel_write_sample_lines(sample_bam_ptr_List *samples,FILE *vcf_out) {
	LIST_FOR_EACH_ELEMENT(sample_bam_ptr, samples, first, next, cur) {
		check(fprintf(vcf_out,VCF_SAMPLE,cur->sample_name,platform,protocol)>0,"Error writing file format line to VCF");
	}
	return 0;
 error:
	return 1;
}

int gen_panel_write_VCF_header(sample_bam_ptr_List *samples,FILE *vcf_out) {
	//Print header to VCF output file
	char *contigs = NULL;
	check(fprintf(vcf_out,"%s\n",VCF_FILE_FORMAT)>0,"Error writing file format line to VCF");

	char date[50];
	time_t current_time;
	time(&current_time);
	strftime(date,sizeof(date),"%Y%m%d",localtime(&current_time));

	check(fprintf(vcf_out,VCF_FILE_DATE,date)>0,"Error writing file date line to VCF");
	check(fprintf(vcf_out,"%s%s%s\n",VCF_SOURCE,date,VCF_SOURCE_2)>0,"Error writing source line to VCF");
	contigs = output_generate_reference_contig_lines(List_first(samples)->bam_file, species_vers, species);
	check(contigs!=NULL,"Error calculating contig lines for VCF");
	check(fprintf(vcf_out,"%s",contigs)>0,"Error writing contig lines to VCF");
	check(fprintf(vcf_out,"%s",VCF_INFO)>0,"Error writing info line to VCF");
	check(fprintf(vcf_out,"%s",VCF_FORMAT_GENO)>0,"Error writing format line to VCF");
	check(fprintf(vcf_out,"%s",VCF_A_ALL)>0,"Error writing A allele format line to VCF");
	check(fprintf(vcf_out,"%s",VCF_C_ALL)>0,"Error writing C allele format line to VCF");
	check(fprintf(vcf_out,"%s",VCF_G_ALL)>0,"Error writing G allele format line to VCF");
	check(fprintf(vcf_out,"%s",VCF_T_ALL)>0,"Error writing T allele format line to VCF");
	check(gen_panel_write_sample_lines(samples,vcf_out)==0,"Error writing sample lines to VCF");
	check(fprintf(vcf_out,VCF_HEADER,sample_names)>0,"Error writing file format line to VCF");
	free(contigs);
	return 0;
 error:
	if(contigs) free(contigs);
	return 1;
}

sample_bam_ptr_List *gen_panel_get_list_of_samples_and_locs() {
	//Create list of sample names & bam files, at the same time converting sample_names to a tab separated string.
	sample_bam_ptr_List *samples = NULL;
	char *sampchar = NULL;
	char *bamchar = NULL;
	samples = sample_bam_ptr_List_create();
	//Iterate through list of sample_names using comma as separator
	char *new_sample_names = malloc((strlen(sample_names)+1) * sizeof(char));
	check_mem(new_sample_names);
	sampchar = strtok(sample_names,",");
	while(sampchar != NULL){
		struct sample_bam *this_samp = malloc(sizeof(struct sample_bam));
		check_mem(this_samp);
		this_samp->sample_name = malloc(sizeof(char) * (strlen(sampchar) + 1));
		check_mem(this_samp->sample_name);
		strcpy(this_samp->sample_name,sampchar);
		check(this_samp->sample_name != NULL,"Error copying sample name to object.");
		sampchar = strtok(NULL,",");
		sample_bam_ptr_List_push(samples,this_samp);
	}
	free(sampchar);

	bamchar = strtok(bam_file_locs,",");
	LIST_FOR_EACH_ELEMENT_MORE(sample_bam_ptr, samples, first, next, cur, more) {
		if(bamchar == NULL){
			sentinel("List of sample names and bam locations are not equal in length.");
		}
		cur->bam_file = malloc(sizeof(char) * (strlen(bamchar) + 1));
		check_mem(cur->bam_file);
		char *tmp = strcpy(cur->bam_file,bamchar);
		check(tmp != NULL,"Error copying bam location to object.");
		if(check_exist(cur->bam_file) != 1){
			printf("Bam file %s does not appear to exist.\n",cur->bam_file);
			gen_panel_print_usage(1);
		}
		char *chk=NULL;
		chk=strcat(new_sample_names,cur->sample_name);
		check(chk!=NULL,"Error appending new sample names.");
		if (more) {
			strcat(new_sample_names,"\t");
		}
		bamchar = strtok(NULL,",");
	}
	strcpy(sample_names,new_sample_names);
	free(new_sample_names);
	free(bamchar);
	return samples;
 error:
	if(new_sample_names) free(new_sample_names);
	if(bamchar) free(bamchar);
	if(samples) sample_bam_ptr_List_destroy(samples);
	return NULL;
}

void gen_panel_clear_pileups(sample_bam_ptr_List *samples) {
	LIST_FOR_EACH_ELEMENT(sample_bam_ptr, samples, first, next, this) {
		if(this->holder){
			if(this->holder->base_counts){
				int i=0;
				for(i=0;i<this->holder->base_counts_size;i++){
					if(this->holder->base_counts[i]) free(this->holder->base_counts[i]);
				}
				free(this->holder->base_counts);
			}
			free(this->holder->bam_access_bases);
			free(this->holder);
		}
	}
	return;
}

int gen_panel_generate_pileups_for_segment(char *ref_file_loc, char *chr_name, int start, int end, sample_bam_ptr_List *samples, FILE *vcf_out) {
	char *ref_seq = fai_access_get_ref_seqeuence_for_pos(ref_file_loc,chr_name,start,end);
	check(ref_seq != NULL,"Error retrieving reference sequence for section %s:%d-%d.",chr_name,start,end);
	fprintf(stdout,"PILEUP REGION: %s:%d-%d\n",chr_name,start,end);
	//Pileup and counts.
	//Iterate through each sample for these locations and write that sample to output line.
	LIST_FOR_EACH_ELEMENT(sample_bam_ptr, samples, first, next, cur) {
		cur->holder = bam_access_get_by_position_counts(cur->bam_file, chr_name, start, end);
		check(cur->holder != NULL,"Error accessing by position counts for sample %s.",cur->sample_name);
	}

	//Now output each position
	int i=0;
	for(i=0;i<((end-start)+1);i++){
		char ref_base = toupper(ref_seq[i]);
		if(!(ref_base == 'A' || ref_base == 'C' || ref_base == 'G' || ref_base == 'T')) continue;
		int sum = 0;
		int write = fprintf(vcf_out,VCF_VAR_LINE_START,
				    chr_name,
				    start+i,
				    ref_base);
		check(write>0,"Error writing VCF beginning of variant line for position %d.",start+i);
		char *last_of_line = NULL;
		int sizeoflist = List_count(samples);
		last_of_line = malloc(sizeof(char) * ((50*sizeoflist)+1));
		check_mem(last_of_line);
		strcpy(last_of_line,"");
		LIST_FOR_EACH_ELEMENT_MORE(sample_bam_ptr, samples, first, next, this, more) {
			if(this->holder->base_counts[i] != 0){
				sum += (this->holder->base_counts[i][0]+
					this->holder->base_counts[i][1]+
					this->holder->base_counts[i][2]+
					this->holder->base_counts[i][3]);
				char tmp[128];
				int chk = 0;
				chk = sprintf(	tmp,VCF_VAR_LINE,
						this->holder->base_counts[i][0],
						this->holder->base_counts[i][1],
						this->holder->base_counts[i][2],
						this->holder->base_counts[i][3]);
				check(chk>0,"Error copying allele counts string.");
				strcat(last_of_line,tmp);
			}else{
				strcat(last_of_line,"-");
			}
			if (more) {
				strcat(last_of_line,"\t");
			}
		}
		strcat(last_of_line,"\n");
		write = 0;
		write = fprintf(vcf_out,VCF_VAR_LINE_MID,sum);
		check(write>0,"Error writing VCF middle of variant line for position %d.",start+i);
		write = 0;
		write = fprintf(vcf_out,"%s",last_of_line);
		check(write>0,"Error writing VCF end of variant line for position %d.",start+i);
		free(last_of_line);
	}


	gen_panel_clear_pileups(samples);
	free(ref_seq);
	return 0;

 error:
	gen_panel_clear_pileups(samples);
	if(ref_seq) free(ref_seq);
	return 1;

}

void gen_panel_destroy_sample_list(sample_bam_ptr_List *samples) {
	LIST_FOR_EACH_ELEMENT(sample_bam_ptr, samples, first, next, cur) {
		free(cur->bam_file);
		free(cur->sample_name);
		free(cur);
	}
	sample_bam_ptr_List_destroy(samples);
	return;
}

int main(int argc, char *argv[]){
	gen_panel_setup_options(argc,argv);

	//Get split section from file given the index.
	char chr_name[50];
	int start_zero_based = 0;
	int stop = 0;
	sample_bam_ptr_List *samples = NULL;

	split_access_get_section_from_index(split_file_loc,chr_name,&start_zero_based,&stop,idx);
	check(stop > 0,"Error fetching region from split file.");
	check(chr_name != NULL, "Error fetching region from split file.");

	samples = gen_panel_get_list_of_samples_and_locs();
	check(samples != NULL,"Error retrieving sample and bam location list from command line.");

	//Open output file
	//Open the config file and do relevant things
	FILE *vcf_out = fopen(out_file_loc,"w");
	check(vcf_out != NULL,"Failed to open VCF file '%s' for writing.",out_file_loc);

	check(gen_panel_write_VCF_header(samples,vcf_out)==0,"Error writing VCF header to file.");

	//We have our location. Let's pileup and output to file.
	bam_access_min_base_qual(min_base_qual);
	printf("Section for analysis %s:%d-%d\n",chr_name,start_zero_based+1,stop);

	int start = start_zero_based+1;
	int end = (start + MAX_SECTION_SIZE) - 1;
	if(end > stop) end = stop;
	while(end<=stop && start<=stop){
		check(gen_panel_generate_pileups_for_segment(ref_file_loc, chr_name, start, end, samples, vcf_out)==0,"Error generating pileups for region %s:%d-%d.",chr_name,start_zero_based+1,stop);
		start = end+1;
		end = (start + MAX_SECTION_SIZE) - 1;
		if(end > stop) end = stop;
	}
	//Close VCF output file.
	fclose(vcf_out);
	gen_panel_destroy_sample_list(samples);
	return 0;
error:
	if(samples) gen_panel_destroy_sample_list(samples);
	return 1;
}
