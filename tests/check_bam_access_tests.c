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

#include <stdlib.h>
#include <check.h>
#include <bam_access.h>
#include <alg_bean.h>
#include <ctype.h>
#include <List.h>
#include <ctype.h>

#include "check_alg_bean_tests.h"

char *bam_file = "../testData/mt.bam";
char *cram_file = "../testData/mt.cram";
char *test_wt_bam = "../testData/test_wt.bam";
char *test_mt_bam = "../testData/test_mt.bam";
char *test_wt_cram = "../testData/test_wt.cram";
char *test_mt_cram = "../testData/test_mt.cram";
char *mut_wt_bam = "../testData/testing_wt.bam";
char *mut_mt_bam = "../testData/testing_mt.bam";
char *mut_wt_cram="../testData/testing_wt.cram";
char *mut_mt_cram = "../testData/testing_mt.cram";
char *bam_access_test_fai_out = NULL;
int exp_contig_length = 252;

int bam_access_compare_read_pos_t_test(const void *in_a, const void *in_b){
	const read_pos_t *a = in_a;
	const read_pos_t *b = in_b;

	if(a->ref_pos > b->ref_pos) return 1;
	if(a->ref_pos < b->ref_pos) return -1;
	return 0;
}

START_TEST (test_bam_access_openbams_close_bams){
    ck_assert_int_eq(bam_access_openbams(bam_file,bam_file,bam_access_test_fai_out),0);
	bam_access_closebams();
}
END_TEST

START_TEST (test_bam_access_openbams_close_bams_cram){
    ck_assert_int_eq(bam_access_openbams(cram_file,cram_file,bam_access_test_fai_out),0);
	bam_access_closebams();
}
END_TEST

START_TEST (test_list_algos_List_insert_sorted){
    List *li = List_create();
	read_pos_t *p1 = malloc(sizeof(read_pos_t));
	read_pos_t *p2 = malloc(sizeof(read_pos_t));
	read_pos_t *p3 = malloc(sizeof(read_pos_t));
	read_pos_t *p4 = malloc(sizeof(read_pos_t));
	read_pos_t *p5 = malloc(sizeof(read_pos_t));
	read_pos_t *p6 = malloc(sizeof(read_pos_t));
	read_pos_t *p7 = malloc(sizeof(read_pos_t));
	unsigned long int on = 1;
	unsigned long int tw = 2;
	unsigned long int th = 3;
	unsigned long int fo = 4;
	unsigned long int fi = 5;
	unsigned long int si = 6;
	unsigned long int se = 7;
	p1->ref_pos = si;
	p2->ref_pos = se;
	p3->ref_pos = tw;
	p4->ref_pos = fo;
	p5->ref_pos = th;
	p6->ref_pos = fi;
	p7->ref_pos = on;
	List_insert_sorted(li, p1, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p2, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p3, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p4, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p5, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p6, (List_compare) bam_access_compare_read_pos_t_test);
	List_insert_sorted(li, p7, (List_compare) bam_access_compare_read_pos_t_test);
    ck_assert(li!=NULL);
    ck_assert(li->first!=NULL);
    ck_assert(li->last!=NULL);
    ck_assert_int_eq(List_count(li), 7);
	unsigned long int i=1;
	LIST_FOREACH(li, first, next, cur) {
        ck_assert_uint_eq(((read_pos_t *)cur->value)->ref_pos,i);
		i++;
	}
	List_destroy(li);
    free(p1);
	free(p2);
	free(p3);
	free(p4);
	free(p5);
	free(p6);
	free(p7);
}
END_TEST

START_TEST (test_bam_access_get_count_for_region){
    bam_access_openbams(bam_file,bam_file,NULL);
	char *chr_name = "22";
	int start = 17619559;
	int stop = 17619559;
	int count_got = bam_access_get_count_for_region(chr_name,start,stop);
	ck_assert_int_eq(count_got,2);
	bam_access_closebams();

	//Now try with more reads expected...
	chr_name = "1";
	start = 124945;
	stop = 124945;
	count_got = 0;
	bam_access_openbams(test_wt_bam,test_mt_bam,bam_access_test_fai_out);
	count_got = bam_access_get_count_for_region(chr_name,start,stop);
	ck_assert_int_eq(count_got,86);
	bam_access_closebams();
}
END_TEST

START_TEST (test_bam_access_get_count_for_region_cram){
    bam_access_openbams(cram_file,cram_file,NULL);
	char *chr_name = "22";
	int start = 17619559;
	int stop = 17619559;
	int count_got = bam_access_get_count_for_region(chr_name,start,stop);
	ck_assert_int_eq(count_got,2);
	bam_access_closebams();

	//Now try with more reads expected...
	chr_name = "1";
	start = 124945;
	stop = 124945;
	count_got = 0;
	bam_access_openbams(test_wt_cram,test_mt_cram,bam_access_test_fai_out);
	count_got = bam_access_get_count_for_region(chr_name,start,stop);
	ck_assert_int_eq(count_got,86);
	bam_access_closebams();
}
END_TEST

START_TEST (test_bam_access_get_reads_at_this_pos){
    bam_access_openbams(bam_file,bam_file,bam_access_test_fai_out);
	char *chr_name = "22";
	int start = 17619559;
	int stop = 17619559;
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(bam_file,bam_file);
	List *got = bam_access_get_reads_at_this_pos(chr_name, start, stop,1,test_bean);
	ck_assert_int_eq(List_count(got),2);
	bam_access_closebams();

	alg_bean_destroy(test_bean);
	chr_name = "1";
	start = 124945;
	stop = 124945;
	bam_access_openbams(test_wt_bam,test_mt_bam,bam_access_test_fai_out);
	test_bean = alg_bean_generate_default_alg_bean(test_wt_bam,test_mt_bam);
	got = bam_access_get_reads_at_this_pos(chr_name, start, stop,1,test_bean);
	ck_assert_int_eq(List_count(got),72);
	bam_access_closebams();

	alg_bean_destroy(test_bean);

	//1	192462357	958902214	C	A	.	PASS	DP=77;MP=1.0e+00;GP=2.2e-06;TG=CC/AAC;TP=5.9e-01;SG=CC/AAA;SP=4.1e-01;VT=Sub	GT:AA:CA:GA:TA:PM	0/0:0:25:0:0:0.0e+00	0/1:37:15:0:0:7.1e-01
	bam_access_openbams(mut_wt_bam,mut_mt_bam,bam_access_test_fai_out);
	test_bean = alg_bean_generate_default_alg_bean(mut_wt_bam,mut_mt_bam);
	got = bam_access_get_reads_at_this_pos("1", 192462357, 192462357,1,test_bean);
	ck_assert_int_eq(List_count(got),77);
	int norm_count = 0;
	int tum_count = 0;
	LIST_FOREACH(got, first, next, cur){
		read_pos_t *rp = (read_pos_t *)cur->value;
		if(rp->normal==1){
			ck_assert_int_eq(toupper(seq_nt16_str[rp->called_base]),'C');
			norm_count++;
		}else{
			ck_assert_msg(toupper(seq_nt16_str[rp->called_base])=='C'||
                        toupper(seq_nt16_str[rp->called_base])=='A',"Wrong tumour called base.");
			tum_count++;
		}
	}
	ck_assert_int_eq(norm_count,25);
	ck_assert_int_eq(tum_count,52);
	alg_bean_destroy(test_bean);
}
END_TEST

START_TEST (test_bam_access_get_reads_at_this_pos_cram){
    bam_access_openbams(cram_file,cram_file,bam_access_test_fai_out);
	char *chr_name = "22";
	int start = 17619559;
	int stop = 17619559;
	fprintf(stderr,"generate alg_bean cram_file\n");
	alg_bean_t *test_bean = alg_bean_generate_default_alg_bean(cram_file,cram_file);
	fprintf(stderr,"generate alg_bean get reads\n");
	List *got = bam_access_get_reads_at_this_pos(chr_name, start, stop,1,test_bean);
	fprintf(stderr,"generate alg_bean got reads\n");
	ck_assert_int_eq(List_count(got),2);
	bam_access_closebams();

	alg_bean_destroy(test_bean);
	chr_name = "1";
	start = 124945;
	stop = 124945;
	bam_access_openbams(test_wt_cram,test_mt_cram,bam_access_test_fai_out);
	fprintf(stderr,"generate alg_bean test_mt_cram\n");
	test_bean = alg_bean_generate_default_alg_bean(test_wt_cram,test_mt_cram);
	fprintf(stderr,"generate alg_bean get reads\n");
	got = bam_access_get_reads_at_this_pos(chr_name, start, stop,1,test_bean);
	fprintf(stderr,"generate alg_bean got reads\n");
	ck_assert_int_eq(List_count(got),72);
	bam_access_closebams();

	alg_bean_destroy(test_bean);

	//1	192462357	958902214	C	A	.	PASS	DP=77;MP=1.0e+00;GP=2.2e-06;TG=CC/AAC;TP=5.9e-01;SG=CC/AAA;SP=4.1e-01;VT=Sub	GT:AA:CA:GA:TA:PM	0/0:0:25:0:0:0.0e+00	0/1:37:15:0:0:7.1e-01
	bam_access_openbams(mut_wt_cram,mut_mt_cram,bam_access_test_fai_out);
	test_bean = alg_bean_generate_default_alg_bean(mut_wt_cram,mut_mt_cram);
	got = bam_access_get_reads_at_this_pos("1", 192462357, 192462357,1,test_bean);
	ck_assert_int_eq(List_count(got),77);
	int norm_count = 0;
	int tum_count = 0;
	LIST_FOREACH(got, first, next, cur){
		read_pos_t *rp = (read_pos_t *)cur->value;
		if(rp->normal==1){
			ck_assert_int_eq(toupper(seq_nt16_str[rp->called_base]),'C');
			norm_count++;
		}else{
			ck_assert_msg(toupper(seq_nt16_str[rp->called_base])=='C'||
                        toupper(seq_nt16_str[rp->called_base])=='A',"Wrong tumour called base.");
			tum_count++;
		}
	}
	ck_assert_int_eq(norm_count,25);
	ck_assert_int_eq(tum_count,52);
	alg_bean_destroy(test_bean);
}
END_TEST

START_TEST (test_bam_access_get_lane_list_from_header){
    List *lanes = bam_access_get_lane_list_from_header(bam_file,"0");
	ck_assert_int_eq(List_count(lanes),2);
}
END_TEST

START_TEST (test_bam_access_get_lane_list_from_header_cram){
	List *lanes = bam_access_get_lane_list_from_header(cram_file,"0");
	ck_assert_int_eq(List_count(lanes),2);
}
END_TEST

START_TEST (test_bam_access_sample_name_platform_from_header){
    //char *bam_file,char *sample, char *plat);
	char *sample = malloc(sizeof(char) * 250);
	char *plat = malloc(sizeof(char) * 250);
	strcpy(plat,".");
	char *exp_sampl = "TUMOURa";
	char *exp_plat = "HiSeq";
	//@RG	ID:1288335	PL:HiSeq	PU:9413_2	LB:TUMOURa 6766555_28085	PI:453	MI:603	DS:short	PG:1288335	SM:TUMOURa	CN:SANGER
	char *res = bam_access_sample_name_platform_from_header(test_mt_bam,sample,plat);
	ck_assert(res!=NULL);
	ck_assert_str_eq(exp_sampl,sample);
	ck_assert_str_eq(exp_plat,plat);
    free(sample);
    free(plat);
}
END_TEST

START_TEST (test_bam_access_sample_name_platform_from_header_cram){
    //char *bam_file,char *sample, char *plat);
	char *sample = malloc(sizeof(char) * 250);
	char *plat = malloc(sizeof(char) * 250);
	strcpy(plat,".");
	char *exp_sampl = "TUMOURa";
	char *exp_plat = "HiSeq";
	//@RG	ID:1288335	PL:HiSeq	PU:9413_2	LB:TUMOURa 6766555_28085	PI:453	MI:603	DS:short	PG:1288335	SM:TUMOURa	CN:SANGER
	char *res = bam_access_sample_name_platform_from_header(test_mt_cram,sample,plat);
	ck_assert(res!=NULL);
	ck_assert_str_eq(exp_sampl,sample);
	ck_assert_str_eq(exp_plat,plat);
    free(sample);
    free(plat);
}
END_TEST

START_TEST (test_bam_access_get_contigs_from_bam){
    int test = 0;
    char *assembly = NULL;
    char *spp = NULL;
	List *contigs = bam_access_get_contigs_from_bam(test_mt_bam, assembly, spp, &test);
	ck_assert_int_eq(List_count(contigs),25);
    ck_assert_msg(strcmp(((ref_seq_t *)(contigs->first->value))->ass,"37")==0,"assembly not correct. Exp: %s\tgot: %s","37",assembly);
    ck_assert_msg(strcmp(((ref_seq_t *)(contigs->first->value))->spp,"HUMAN")==0,"species not correct. Exp: %s\tgot: %s","37",assembly);
    ck_assert_int_eq(test,exp_contig_length);
}
END_TEST

START_TEST (test_bam_access_get_contigs_from_bam_cram){
    int test = 0;
	List *contigs = bam_access_get_contigs_from_bam(test_mt_cram, NULL, NULL, &test);
	ck_assert_int_eq(List_count(contigs),25);
    ck_assert_int_eq(test,exp_contig_length);
}
END_TEST

START_TEST (test_bam_access_get_contigs_from_bam_no_spp){
    int test = 0;
	List *contigs = bam_access_get_contigs_from_bam(test_mt_bam, "ASSEMBLY", "SPP", &test);
	ck_assert_int_eq(List_count(contigs),25);
    ck_assert_int_eq(test,exp_contig_length);
}
END_TEST

START_TEST (test_bam_access_get_contigs_from_bam_no_spp_cram){
    int test = 0;
	List *contigs = bam_access_get_contigs_from_bam(test_mt_cram, "ASSEMBLY", "SPP", &test);
	ck_assert_int_eq(List_count(contigs),25);
    ck_assert_int_eq(test,exp_contig_length);
}
END_TEST

START_TEST (test_bam_access_check_bam_flags){
    bam1_t *b = bam_init1();
	b->core.qual = 60;
	b->core.flag = 99; //mapped, proper pair, mate rev strand, 1st in pair
	//Check pass
	ck_assert_int_eq(bam_access_check_bam_flags(b),1);
	//check fail on qual
	b->core.qual = 0;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	b->core.qual = 60;
	//Check fail on unmapped
	b->core.flag = 101;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	//Check fail on secondary mapping
	b->core.flag = 355;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	//Check fail on QCfail
	b->core.flag = 611;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	//Check fail on 2048 flag
	b->core.flag = 2147;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	//check fail on duplicate
	b->core.flag = 1123;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	//check pass duplicate
	bam_access_include_dup(1);
	b->core.flag = 1123;
	ck_assert_int_eq(bam_access_check_bam_flags(b),1);
	//check pass single end
	b->core.flag = 16;
	bam_access_include_dup(0);
	bam_access_include_se(1);
	ck_assert_int_eq(bam_access_check_bam_flags(b),1);
	//Check fail single end on proper pair
	b->core.flag = 81;
	bam_access_include_se(0);
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	//check fail single end on mate unmapped
	b->core.flag = 91;
	ck_assert_int_eq(bam_access_check_bam_flags(b),0);
	bam_destroy1(b);
}
END_TEST

START_TEST (test_bam_access_populate_file_index){
    samFile *sf = NULL;
	hts_idx_t *idx = NULL;
	sf = bam_access_populate_file(test_mt_bam, bam_access_test_fai_out );
	ck_assert(sf!=NULL);
	idx = bam_access_populate_file_index(sf,test_mt_bam );
	ck_assert(idx!=NULL);
	sam_close(sf);
	hts_idx_destroy(idx);
}
END_TEST

START_TEST (test_bam_access_populate_file_index_cram){
    samFile *sf = NULL;
	hts_idx_t *idx = NULL;
	sf = bam_access_populate_file(test_mt_cram , bam_access_test_fai_out );
	ck_assert(sf!=NULL);
	idx = bam_access_populate_file_index(sf,test_mt_cram );
	ck_assert(idx!=NULL);
	sam_close(sf);
	hts_idx_destroy(idx);
}
END_TEST

START_TEST (test_bam_access_get_hts_itr){
    htsFile *sf = NULL;
	hts_idx_t *idx = NULL;
	hts_itr_t *itr = NULL;
	char *chr = "22";
	uint32_t from = 17619559;
	uint32_t to = 17619559;
	sf = bam_access_populate_file(test_mt_bam, bam_access_test_fai_out  );
	ck_assert(sf!=NULL);
	idx = bam_access_populate_file_index(sf,test_mt_bam );
	ck_assert(idx!=NULL);
	itr = bam_access_get_hts_itr(sf, idx, chr, from, to);
	ck_assert(itr!=NULL);
	hts_close(sf);
	hts_idx_destroy(idx);
	hts_itr_destroy(itr);
}
END_TEST

START_TEST (test_bam_access_get_hts_itr_cram){
    htsFile *sf = NULL;
	hts_idx_t *idx = NULL;
	hts_itr_t *itr = NULL;
	char *chr = "22";
	uint32_t from = 17619559;
	uint32_t to = 17619559;
	sf = bam_access_populate_file(cram_file, bam_access_test_fai_out  );
	ck_assert(sf!=NULL);
	idx = bam_access_populate_file_index(sf,cram_file );
	ck_assert(idx!=NULL);
	itr = bam_access_get_hts_itr(sf, idx, chr, from, to);
	ck_assert(itr!=NULL);
	hts_close(sf);
	hts_idx_destroy(idx);
	hts_itr_destroy(itr);
}
END_TEST

Suite *check_bam_access_tests_suite(void){
    Suite *s;
    TCase *tc_bam_access;

    s = suite_create("bam access");

    /* Core test case */
    tc_bam_access = tcase_create("bam access testing");

    tcase_add_test(tc_bam_access, test_bam_access_openbams_close_bams);
    tcase_add_test(tc_bam_access, test_bam_access_openbams_close_bams_cram);
    tcase_add_test(tc_bam_access, test_list_algos_List_insert_sorted); 
    tcase_add_test(tc_bam_access, test_bam_access_get_count_for_region); 
    tcase_add_test(tc_bam_access, test_bam_access_get_count_for_region_cram);
    tcase_add_test(tc_bam_access, test_bam_access_get_reads_at_this_pos);
    tcase_add_test(tc_bam_access, test_bam_access_get_reads_at_this_pos_cram);
    tcase_add_test(tc_bam_access, test_bam_access_get_lane_list_from_header);
    tcase_add_test(tc_bam_access, test_bam_access_get_lane_list_from_header_cram);
    tcase_add_test(tc_bam_access, test_bam_access_sample_name_platform_from_header);
    tcase_add_test(tc_bam_access, test_bam_access_sample_name_platform_from_header_cram);
    tcase_add_test(tc_bam_access, test_bam_access_get_contigs_from_bam);
    tcase_add_test(tc_bam_access, test_bam_access_get_contigs_from_bam_cram);
    tcase_add_test(tc_bam_access, test_bam_access_get_contigs_from_bam_no_spp);
    tcase_add_test(tc_bam_access, test_bam_access_get_contigs_from_bam_no_spp_cram);
    tcase_add_test(tc_bam_access, test_bam_access_check_bam_flags);
    tcase_add_test(tc_bam_access, test_bam_access_populate_file_index);
    tcase_add_test(tc_bam_access, test_bam_access_populate_file_index_cram);
    tcase_add_test(tc_bam_access, test_bam_access_get_hts_itr);
    tcase_add_test(tc_bam_access, test_bam_access_get_hts_itr_cram);

    suite_add_tcase(s, tc_bam_access);

    return s;
}
