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
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <config_file_access.h>
#include <dbg.h>

const char *CWD="CWD";
const char *MUT_TUM = "MUT_BAM";
const char *NORM_TUM = "NORM_BAM";
const char *REF_INDEX = "REF_IDX";
const char *IGN_FILE = "IGNORE";
const char *ALG_FILE = "ALG_FILE";
const char *RES_DIR = "RESULT_DIR";
const char *SPLIT_FILE = "SPLIT_FILE";
const char *SW_KEY = "SW";
const char *SE_KEY = "SE";
const char *DUP_KEY = "DUP";
const char *VERSION_KEY = "VER";
const char *NORM_CN_KEY = "NORMCN";
const char *TUM_CN_KEY = "TUMCN";


int config_file_access_read_config_file(FILE *file, char *tum_bam_file, char *norm_bam_file, char *ref_idx,
			char *ignore_regions_file, char *alg_bean_loc, char *results, char *list_loc, int *includeSW,
			int *includeSingleEnd, int *includeDups, char *version, char *norm_cn, char *tum_cn)
{

	int found_norm_cn = 0;
	int found_tum_cn = 0;
	char line [ 3074 ];
	while ( fgets(line,sizeof(line),file) != NULL ){
		char key[ 16];
		char value [2048];
		int chk = sscanf(line,"%[^=]=%s",key,value);
		check(chk==2,"Couldn't resolve key - value pair from text %s in config file.",line);
		if(strcmp(CWD,key)==0){
			size_t size = 500;
			char *cur_wd = (char *) malloc(size);
			check_mem(cur_wd);
			char *chk = getcwd(cur_wd,size);
			check(chk!=NULL,"Error retrieving CWD when reading config file.");
			if(strcmp(cur_wd,value)!=0){
				sentinel("Your current working directory '%s' is not the same as the directory you setup CaVEMan: '%s'. Please change to that directory and try again.",
										cur_wd,value);
			}
		}else if(strcmp(MUT_TUM,key)==0){
			strcpy(tum_bam_file,value);
		}else if(strcmp(NORM_TUM,key)==0){
			strcpy(norm_bam_file,value);
		}else if(strcmp(REF_INDEX,key)==0){
			strcpy(ref_idx,value);
		}else if(strcmp(IGN_FILE,key)==0){
			strcpy(ignore_regions_file,value);
		}else if(strcmp(ALG_FILE,key)==0){
			strcpy(alg_bean_loc,value);
		}else if(strcmp(RES_DIR,key)==0){
			strcpy(results,value);
		}else if(strcmp(SPLIT_FILE,key)==0){
			strcpy(list_loc,value);
		}else if(strcmp(SW_KEY,key)==0){
			if(atoi(value) != 0){
				*includeSW = 1;
			}else{
				*includeSW = 0;
			}
		}else if(strcmp(SE_KEY,key)==0){
			if(atoi(value) != 0){
				*includeSingleEnd = 1;
			}else{
				*includeSingleEnd = 0;
			}
		}else if(strcmp(DUP_KEY,key)==0){
			if(atoi(value) != 0){
				*includeDups = 1;
			}else{
				*includeDups = 0;
			}
		}else if(strcmp(VERSION_KEY,key)==0){
			strcpy(version,value);
		}else if(strcmp(NORM_CN_KEY,key)==0){
			strcpy(norm_cn,value);
		}else if(strcmp(TUM_CN_KEY,key)==0){
			strcpy(tum_cn,value);
		}else{
			sentinel("Unrecognised key in config file '%s'.",key);
		}
	}

	if(found_norm_cn==0){
		norm_cn = NULL;
	}
	if(found_tum_cn==0){
		tum_cn = NULL;
	}
	return 0;
error:
	return -1;
}

int resolve_real_path(char *to_be_resolved,char *to_be_allocated){
	char *ptr = realpath(to_be_resolved,to_be_allocated);
	check(ptr!=NULL,"Checking real path was assigned for %s.",to_be_resolved);
	return 0;
error:
	return -1;
}

int config_file_access_write_config_file(FILE *file, char *tum_bam_file, char *norm_bam_file, char *ref_idx,
							char *ignore_regions_file, char *alg_bean_loc, char *results, char *list_loc, int includeSW,
							int includeSingleEnd, int includeDups, char *norm_cn, char *tum_cn)
{
	assert(tum_bam_file != NULL && norm_bam_file != NULL && ref_idx != NULL && ignore_regions_file != NULL
				&& alg_bean_loc != NULL && results != NULL && list_loc != NULL);
	char *curr_wd = NULL;

	//Potentially non existant, so create it first
	int dir = mkdir(results,S_IRWXU);
	if(dir==0){
		printf("Created results directory '%s'.\n",results);
	}

	size_t size = 500;
	curr_wd = (char *) malloc(size);
	check_mem(curr_wd);
	char *chk = getcwd(curr_wd,size);
	check(chk!=NULL,"Error retrieving cwd for config file.");
	int res = fprintf(file,"%s=%s\n",CWD,curr_wd);
	check(res>=0,"Error writing cwd to config file.");
	res = fprintf(file,"%s=%s\n",MUT_TUM,tum_bam_file);
	check(res>=0,"Error writing mutant bam file loc to config file.");
	res = fprintf(file,"%s=%s\n",NORM_TUM,norm_bam_file);
	check(res>=0,"Error writing normal bam file loc to config file.");
	res = fprintf(file,"%s=%s\n",REF_INDEX,ref_idx);
	check(res>=0,"Error writing reference index file loc to config file.");
	res = fprintf(file,"%s=%s\n",IGN_FILE,ignore_regions_file);
	check(res>=0,"Error writing ignore regions file loc to config file.");
	res = fprintf(file,"%s=%s\n",ALG_FILE,alg_bean_loc);
	check(res>=0,"Error writing alg bean file loc to config file.");
	res = fprintf(file,"%s=%s\n",RES_DIR,results);
	check(res>=0,"Error writing results directory to config file.");
	res = fprintf(file,"%s=%s\n",SPLIT_FILE,list_loc);
	check(res>=0,"Error writing split sections file loc to config file.");
	res = fprintf(file,"%s=%d\n",SW_KEY,includeSW);
	check(res>=0,"Error writing include SW to config file.");
	res = fprintf(file,"%s=%d\n",SE_KEY,includeSingleEnd);
	check(res>=0,"Error writing inclue single end to config file.");
	res = fprintf(file,"%s=%d\n",DUP_KEY,includeDups);
	check(res>=0,"Error writing include duplicates file to config file.");
	if(norm_cn != NULL){
		res = fprintf(file,"%s=%s\n",NORM_CN_KEY,norm_cn);
		check(res>=0,"Error writing normal cn file to config file.");
	}
	if(tum_cn != NULL){
		res = fprintf(file,"%s=%s\n",TUM_CN_KEY,tum_cn);
		check(res>=0,"Error writing tumour cn file to config file.");
	}
	//Finally print the version for version checking.
	res = fprintf(file,"%s=%s\n",VERSION_KEY,CAVEMAN_VERSION);
	check(res>=0,"Error writing version information to config file.");
	return 0;
error:
	if(curr_wd) free(curr_wd);
	return -1;
}
