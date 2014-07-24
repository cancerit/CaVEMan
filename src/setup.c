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
#include <alg_bean.h>
#include <config_file_access.h>

static int includeSW = 0;
static int includeSingleEnd = 0;
static int includeDups = 0;
static char *tum_bam_file = NULL;
static char *norm_bam_file = NULL;
static char *results = "results";
static char *ref_idx = NULL;
static char *list_loc = "splitList";
static char *alg_bean_loc = "alg_bean";
static char *ignore_regions_file = NULL;
static char *CaVEManfg_ini = "caveman.cfg.ini";


void setup_print_usage (int exit_code){
	printf ("Usage: caveman setup -t tum.bam -n norm.bam -r reference.fa.fai -g ignore_regions.tab [-f path] [-l path] [-a path] [-wzu]\n\n");
	printf("-t  --tumour-bam [file]             Location of tumour bam\n");
	printf("-n  --normal-bam [file]             Location of normal bam\n");
	printf("-r  --reference-index [file]        Location of reference fasta index\n");
	printf("-g  --ignore-regions-file [file]    Location of tsv ignore regions file\n\n");
	printf("Optional\n");
	printf("-c  --config-file [file]            File to write caveman run config file [default:'%s']\n",CaVEManfg_ini);
	printf("-f  --results-folder [file]         Folder to write results [default:'%s']\n",results);
	printf("-l  --split-file [file]             File to write list of split sections [default:'%s']\n",list_loc);
	printf("-a  --alg-bean-file [file]          Location to write alg-bean [default:'%s']\n",alg_bean_loc);
	printf("-w  --include-smith-waterman        Include SW mapped reads in the analysis\n");
	printf("-z  --include-single-end            Use single end reads for this analysis\n");
	printf("-u  --include-duplicates            Include reads marked as duplicates in the analysis\n");
	printf("-h	--help                          Display this usage information.\n");
  exit(exit_code);
}

void setup_setup_options(int argc, char *argv[]){
	const struct option long_opts[] =
	{
		{"tumour-bam", required_argument, 0, 't'},
		{"normal-bam", required_argument, 0, 'n'},
		{"results-folder", required_argument, 0, 'f'},
		{"split-file", required_argument, 0, 'l'},
		{"alg-bean-file",required_argument, 0, 'a'},
		{"reference-index", required_argument, 0, 'r'},
		{"include-smith-waterman", no_argument, 0, 'w'},
		{"include-single-end", no_argument, 0, 'z'},
		{"include-duplicates", no_argument, 0, 'u'},
		{"config-file",required_argument, 0,'c'},
		{"ignore-regions-file", required_argument, 0, 'g'},
		{"help", no_argument, 0, 'h'},
		{ NULL, 0, NULL, 0}
   }; //End of declaring opts
   
   int index = 0;
   int iarg = 0;
   
   //Iterate through options
   while((iarg = getopt_long(argc, argv, "t:n:r:f:a:l:g:c:wzuh",
                            								long_opts, &index)) != -1){
   	switch(iarg){
   		case 'h':
         	setup_print_usage(0);
         	break;

         case 'c':
         	CaVEManfg_ini = optarg;
         	break;
     			
         case 't':
         	tum_bam_file = optarg;
            break;
         	
         case 'n':
         	norm_bam_file = optarg;
            break;
            
         case 'g':
         	ignore_regions_file = optarg;
            break;
         
         case 'r':
         	ref_idx = optarg;
         	break;
      	
      	case 'f':
      		results = optarg;
      		break;
      		
      	case 'a':
      		alg_bean_loc = optarg;
      		break;
      		
      	case 'l':
      		list_loc = optarg;
      		break;    	

      	case 'w':
      		includeSW = 1;
      		break;
      		
      	case 'z':
      		includeSingleEnd = 1;
      		break;
      		
      	case 'u':
      		includeDups = 1;
      		break;
	
			case '?':
            setup_print_usage (1);
            break;
	
      	default:
      		setup_print_usage (1);
                           
   	}; // End of args switch statement
   	
   }//End of iteration through options
   
   //Do some checking to ensure required arguments were passed
   if(tum_bam_file == NULL || norm_bam_file == NULL || ref_idx  == NULL || ignore_regions_file == NULL ){
   	setup_print_usage(1);
   }
   return;
}

int setup_main(int argc, char *argv[]){
	setup_setup_options(argc,argv);
	
	//Create and write the config file in the current directory. 
	FILE *config_read;
	//Try reading to see if we already have one
	if((config_read = fopen(CaVEManfg_ini,"r")) == 0){
		FILE *config_out = fopen(CaVEManfg_ini,"w");
		check(config_out != NULL,"Error trying to open config file location for write: %s.",CaVEManfg_ini);
		int res = config_file_access_write_config_file(config_out, tum_bam_file, norm_bam_file, ref_idx, 
			ignore_regions_file, alg_bean_loc, results, list_loc, includeSW, includeSingleEnd, includeDups);
		check(res==0,"Problem encountered when writing new config file to to %s.",CaVEManfg_ini);
		res = fclose(config_out);	
		check(res==0,"Error closing config file.");
	}else{
		int res = fclose(config_read);
		check(res==0,"Error closing config file reader.");
		printf("Config file file: '%s' already exists.\nDelete it if you want a new one created, or CaVEMan will reuse this one.\n",CaVEManfg_ini);
	}
	
	//Create the alg bean
	FILE *bean_read;
	if((bean_read = fopen(alg_bean_loc,"r")) == 0){
		FILE *bean_out = fopen(alg_bean_loc,"w");
		check(bean_out != NULL,"Error trying to open alg_bean file location for write: %s.",alg_bean_loc);
		int res = alg_bean_create_default_file(bean_out,norm_bam_file,tum_bam_file);
		check(res==0,"Problem encountered when writing new alg bean to %s.",alg_bean_loc);
		res = fclose(bean_out);
		check(res==0,"Error closing alg bean file.");
	}else{
		int res = fclose(bean_read);
		check(res==0,"Error closing alg bean file reader.");
		printf("Alg bean file: '%s' already exists.\nDelete it if you want a new one created, or CaVEMan will reuse this one.\n",alg_bean_loc);
	}
	
	return 0;
	
error:
	return -1;
		
}