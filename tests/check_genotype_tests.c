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
#include <genotype.h>
#include <List.h>
#include <string.h>

#include "check_genotype_tests.h"

START_TEST(test_genotype_calculate_genotypes){
	char *ref_base = "A";
	int cn = 2;
	List *genos = genotype_calculate_genotypes(cn,ref_base);
	ck_assert_msg(genos != NULL, "Genotypes list was NULL");
	ck_assert_msg(List_count(genos) == 7, "Wrong number of genotypes returned in cn 2 test.");
	cn = 4;
	List *sec_genos = genotype_calculate_genotypes(4,ref_base);
	ck_assert_msg(sec_genos != NULL, "Genotypes list was NULL");
	ck_assert_msg(List_count(sec_genos) == 13, "Wrong number of genotypes returned in cn 4 test.");
	genotype_clear_genotype_cache();
}
END_TEST

START_TEST(test_genotype_hard_copy_genotype_t_list){
	List *new_list = List_create();
	List *old = List_create();
	genotype_t *geno = genotype_init_genotype();
	int count = 10;
	char base = 'C';
	genotype_set_base_count(geno, base, count);
	ck_assert_msg(geno->c_count==10,"Wrong number of C bases recorded");
	List_push(old,geno);
	genotype_hard_copy_genotype_t_list(new_list,old);
	ck_assert_msg(List_count(new_list) == 1, "Wrong number of elements in list.");
	LIST_FOREACH(new_list, first, next, cur){
		//Only one element...
		ck_assert_msg(((genotype_t *)cur->value)->c_count == 10, "Wrong c count in copied value.");
	}
	List_clear_destroy(new_list);
	List_clear_destroy(old);
}
END_TEST

START_TEST(test_genotype_add_base_to_count){
	int chk = 0;
	genotype_t *geno = genotype_init_genotype();
	chk = genotype_add_base_to_count(geno,'C');
	ck_assert_msg(chk==0, "Error adding base to genotype count.");
	ck_assert_msg(geno->c_count == 1, "Incorrect count of C after first add.");
	chk = genotype_add_base_to_count(geno,'C');
	ck_assert_msg(chk==0, "Error adding base to genotype count.");
	ck_assert_msg(geno->c_count == 2, "Incorrect count of C after second add.");
	chk = genotype_add_base_to_count(geno,'T');
	ck_assert_msg(chk==0, "Error adding base to genotype count.");
	ck_assert_msg(geno->t_count == 1, "Incorrect count of T after third add.");
	ck_assert_msg(geno->c_count == 2, "Incorrect count of C after third add.");
	ck_assert_msg(geno->a_count == 0, "Incorrect count of A after third add.");
	ck_assert_msg(geno->g_count == 0, "Incorrect count of G after third add.");
	free(geno);
}
END_TEST

START_TEST(test_genotype_set_base_count){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 2;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->t_count == 2, "Incorrectly set T count.");
	free(geno);
}
END_TEST

START_TEST(test_genotype_get_base_count){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->a_count == 3, "Incorrectly set A count.");
	int got = genotype_get_base_count(geno,base);
	ck_assert_msg(got==base_count,"Incorrect base count returned by function.");
	free(geno);
}
END_TEST

START_TEST(test_genotype_init_genotype){
	genotype_t *geno = genotype_init_genotype();
	ck_assert_msg(geno->t_count == 0, "Incorrect count of T after init.");
	ck_assert_msg(geno->c_count == 0, "Incorrect count of C after init.");
	ck_assert_msg(geno->a_count == 0, "Incorrect count of A after init.");
	ck_assert_msg(geno->g_count == 0, "Incorrect count of G after init.");
	free(geno);
}
END_TEST

START_TEST(test_genotype_copy_genotype){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 2;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->t_count == 2, "Incorrectly set T count.");
	genotype_t *new = genotype_copy_genotype(geno);
	ck_assert_msg(new->a_count == 3, "Incorrectly set A count.");
	ck_assert_msg(new->t_count == 2, "Incorrectly set T count.");
	free(new);
	free(geno);
}
END_TEST

START_TEST(test_genotype_get_var_base_proportion){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->a_count == 1, "Incorrectly set A count.");
	base = 'T';
	base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->t_count == 1, "Incorrectly set T count.");
	double expected = 0.5;
	double got = genotype_get_var_base_proportion(geno, base, 2);
	ck_assert_msg(got==expected,"Incorrect variant base proportion calculated.");
	free(geno);
	geno = genotype_init_genotype();
	base = 'A';
	base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->t_count == 1, "Incorrectly set T count.");
	expected = 0.75;
	got = genotype_get_var_base_proportion(geno, base, 4);
	ck_assert_msg(got==expected,"Incorrect variant base proportion calculated.");
	free(geno);
}
END_TEST

START_TEST(test_genotype_get_genotype_t_as_string){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	ck_assert_msg(geno->t_count == 1, "Incorrectly set T count.");
	char *expected = "AAAT";
	char *got = genotype_get_genotype_t_as_string(geno);
	ck_assert_msg(strcmp(expected,got)==0,"Wrong string returned");
	free(geno);
	free(got);
}
END_TEST

START_TEST(test_genotype_generate_genotype_list_for_cn_and_ref_base){
	//int norm_cn, int tum_cn, char *ref_base
	int norm_cn = 2;
	int tum_cn = 2;
	//Make sure we aren't using the cache...
	char *ref_base = "A";
	genotype_store_t *store = genotype_generate_genotype_list_for_cn_and_ref_base(norm_cn, tum_cn, ref_base);
	ck_assert_msg(store != NULL, "Genotype store not generated properly");
	ck_assert_msg(List_count(store->normal_genos)==7,"Wrong number of normal genotypes.");
	ck_assert_msg(List_count(store->tumour_genos)==7,"Wrong number of tumour genotypes.");

	genotype_destroy_genotype_store(store);
	genotype_clear_genotype_cache();

	tum_cn = 3;
	store = genotype_generate_genotype_list_for_cn_and_ref_base(norm_cn, tum_cn, ref_base);
	ck_assert_msg(store != NULL, "Genotype store not generated properly");
	ck_assert_msg(store->somatic_count==9,"Wrong number of somatic genotypes.");
	ck_assert_msg(store->hom_count==3,"Wrong number of hom snp genotypes.");
	ck_assert_msg(store->het_count==12,"Wrong number of het snp genotypes.");
	ck_assert_msg(List_count(store->tumour_genos)==10,"Wrong number of tumour genotypes.");

	genotype_destroy_genotype_store(store);
	genotype_clear_genotype_cache();

	tum_cn = 4;
	store = genotype_generate_genotype_list_for_cn_and_ref_base(norm_cn, tum_cn, ref_base);
	ck_assert_msg(store != NULL, "Genotype store not generated properly");
	ck_assert_msg(List_count(store->normal_genos)==7,"Wrong number of normal genotypes.");
	ck_assert_msg(List_count(store->tumour_genos)==13,"Wrong number of tumour genotypes.");
	ck_assert_msg(store->somatic_count==12,"Wrong number of somatic genotypes.");

	genotype_destroy_genotype_store(store);
	genotype_clear_genotype_cache();

}
END_TEST

Suite * check_genotype_tests_suite(void){
    Suite *s;
    TCase *tc_genotype;

    s = suite_create("genotype_tests");

    /* Core test case */
    tc_genotype = tcase_create("genptype testing");

    tcase_add_test(tc_genotype, test_genotype_init_genotype);
    tcase_add_test(tc_genotype, test_genotype_calculate_genotypes);
    tcase_add_test(tc_genotype, test_genotype_add_base_to_count);
    tcase_add_test(tc_genotype, test_genotype_hard_copy_genotype_t_list);
    tcase_add_test(tc_genotype, test_genotype_set_base_count);
    tcase_add_test(tc_genotype, test_genotype_get_base_count);
    tcase_add_test(tc_genotype, test_genotype_copy_genotype);
    tcase_add_test(tc_genotype, test_genotype_get_var_base_proportion);
    tcase_add_test(tc_genotype, test_genotype_get_genotype_t_as_string);
    tcase_add_test(tc_genotype, test_genotype_generate_genotype_list_for_cn_and_ref_base);

    suite_add_tcase (s, tc_genotype);

    return s;
}
