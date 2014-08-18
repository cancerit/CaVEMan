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

#include "minunit.h"
#include <bam_access.h>
#include <alg_bean.h>
#include <List.h>

char *bam_file = "testData/mt.bam";
char *test_wt_bam = "testData/test_wt.bam";
char *test_mt_bam = "testData/test_mt.bam";
char *mut_wt_bam = "testData/testing_wt.bam";
char *mut_mt_bam = "testData/testing_mt.bam";

char *test_bam_access_openbams_close_bams(){
	bam_access_openbams(bam_file,bam_file);
	bam_access_closebams();
	return NULL;
}

int bam_access_compare_read_pos_t_test(const void *in_a, const void *in_b){
	const read_pos_t *a = in_a;
	const read_pos_t *b = in_b;

	if(a->ref_pos > b->ref_pos) return 1;
	if(a->ref_pos < b->ref_pos) return -1;
	return 0;
}

char *test_list_algos_List_insert_sorted(){
	List *li = List_create();
	read_pos_t *p1 = malloc(sizeof(read_pos_t));
	read_pos_t *p2 = malloc(sizeof(read_pos_t));
	read_pos_t *p3 = malloc(sizeof(read_pos_t));
	read_pos_t *p4 = malloc(sizeof(read_pos_t));
	read_pos_t *p5 = malloc(sizeof(read_pos_t));
	read_pos_t *p6 = malloc(sizeof(read_pos_t));
	read_pos_t *p7 = malloc(sizeof(read_pos_t));
	p1->ref_pos = 6;
	p2->ref_pos = 7;
	p3->ref_pos = 2;
	p4->ref_pos = 4;
	p5->ref_pos = 3;
	p6->ref_pos = 5;
	p7->ref_pos = 1;
	List_insert_sorted(li, p1, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p2, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p3, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p4, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p5, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p6, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p7, (List_compare) bam_access_compare_read_pos_t_test);
	mu_assert(li != NULL, "NULL list found.");
	mu_assert(li->first != NULL, "NULL list->first found.");
	mu_assert(li->last != NULL, "NULL list->last found.");
	mu_assert(List_count(li) == 7,"Incorrect number of elements in list.");
	int i=1;
	LIST_FOREACH(li, first, next, cur) {
		mu_assert(((read_pos_t *)cur->value)->ref_pos == i,"Incorrect number in sorted order.");
		i++;
	}
	List_destroy(li);
	return NULL;
}

char *test_bam_access_get_count_for_region(){
	bam_access_openbams(bam_file,bam_file);
	char *chr_name = "22";
	int start = 17619559;
	int stop = 17619559;
	int count_got = bam_access_get_count_for_region(chr_name,start,stop);
	mu_assert(count_got==2,"Incorrect number of reads returned.");
	bam_access_closebams();

	//Now try with more reads expected...
	chr_name = "1";
	start = 124945;
	stop = 124945;
	count_got = 0;
	bam_access_openbams(test_wt_bam,test_mt_bam);
	count_got = bam_access_get_count_for_region(chr_name,start,stop);
	mu_assert(count_got==86,"Incorrect number of reads returned in second lookup.");
	bam_access_closebams();
	return NULL;
}

char *test_bam_access_get_reads_at_this_pos(){
	bam_access_openbams(bam_file,bam_file);
	char *chr_name = "22";
	int start = 17619559;
	int stop = 17619559;
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(bam_file,bam_file);
	List *got = bam_access_get_reads_at_this_pos(chr_name, start, stop,1,test_bean);
	mu_assert(List_count(got)==2,"Wrong number of reads fetched from bam file.");
	bam_access_closebams();

	alg_bean_destroy(test_bean);
	chr_name = "1";
	start = 124945;
	stop = 124945;
	bam_access_openbams(test_wt_bam,test_mt_bam);
	test_bean = alg_bean_generate_default_alg_bean(test_wt_bam,test_mt_bam);
	got = bam_access_get_reads_at_this_pos(chr_name, start, stop,1,test_bean);
	mu_assert(List_count(got)==72,"Wrong number of reads fetched from bam file.");
	bam_access_closebams();

	alg_bean_destroy(test_bean);

	//1	192462357	958902214	C	A	.	PASS	DP=77;MP=1.0e+00;GP=2.2e-06;TG=CC/AAC;TP=5.9e-01;SG=CC/AAA;SP=4.1e-01;VT=Sub	GT:AA:CA:GA:TA:PM	0/0:0:25:0:0:0.0e+00	0/1:37:15:0:0:7.1e-01
	bam_access_openbams(mut_wt_bam,mut_mt_bam);
	test_bean = alg_bean_generate_default_alg_bean(mut_wt_bam,mut_mt_bam);
	got = bam_access_get_reads_at_this_pos("1", 192462357, 192462357,1,test_bean);
	mu_assert(List_count(got)==77,"Wrong number of reads fetched from know mutant bam files.");
	int norm_count = 0;
	int tum_count = 0;
	LIST_FOREACH(got, first, next, cur){
		read_pos_t *rp = (read_pos_t *)cur->value;
		if(rp->normal==1){
			mu_assert(rp->called_base=='C',"Wrong called base in normal.");
			norm_count++;
		}else{
			mu_assert(rp->called_base=='C'||rp->called_base=='A',"Wrong tumour called base.");
			tum_count++;
		}
	}
	mu_assert(norm_count==25,"Incorrect numer of normal reads fetched.");
	mu_assert(tum_count==52,"Incorrect numer of tumour reads fetched.");
	alg_bean_destroy(test_bean);
	return NULL;
}

char *test_bam_access_get_lane_list_from_header(){
	List *lanes = bam_access_get_lane_list_from_header(bam_file,"0");
	mu_assert(List_count(lanes)==2,"Wrong number of lanes fetched from bam file.");
	return NULL;
}

char *test_bam_access_sample_name_platform_from_header(){
	//char *bam_file,char *sample, char *plat);
	char *sample = malloc(sizeof(char) * 250);
	char *plat = malloc(sizeof(char) * 250);
	strcpy(plat,".");
	char *exp_sampl = "TUMOURa";
	char *exp_plat = "HiSeq";
	//@RG	ID:1288335	PL:HiSeq	PU:9413_2	LB:TUMOURa 6766555_28085	PI:453	MI:603	DS:short	PG:1288335	SM:TUMOURa	CN:SANGER
	char *res = bam_access_sample_name_platform_from_header(test_mt_bam,sample,plat);
	mu_assert(res!=NULL,"Error trying to fetch sample/platform from bam header.");
	mu_assert(strcmp(exp_sampl,sample)==0,"Samples don't match from fetch.");
	mu_assert(strcmp(exp_plat,plat)==0,"Platforms don't match from fetch.");
	return NULL;
}

char *test_bam_access_get_contigs_from_bam(){
	List *contigs = bam_access_get_contigs_from_bam(test_mt_bam, NULL, NULL);
	mu_assert(List_count(contigs)==25,"Wrong number of contigs fetched from header.");
	return NULL;
}

char *test_bam_access_get_contigs_from_bam_no_spp(){
	List *contigs = bam_access_get_contigs_from_bam(test_mt_bam, "ASSEMBLY", "SPP");
	mu_assert(List_count(contigs)==25,"Wrong number of contigs fetched from header.");
	return NULL;
}


char *all_tests() {
   mu_suite_start();
   mu_run_test(test_bam_access_openbams_close_bams);
   mu_run_test(test_bam_access_get_reads_at_this_pos);
   mu_run_test(test_bam_access_get_count_for_region);
   mu_run_test(test_bam_access_get_lane_list_from_header);
   mu_run_test(test_list_algos_List_insert_sorted);
   mu_run_test(test_bam_access_sample_name_platform_from_header);
   mu_run_test(test_bam_access_get_contigs_from_bam);
   mu_run_test(test_bam_access_get_contigs_from_bam_no_spp);
   return NULL;
}

RUN_TESTS(all_tests);
