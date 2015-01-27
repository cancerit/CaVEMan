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
#include <genotype.h>
#include <string.h>

char *test_genotype_calculate_genotypes(){
	char *ref_base = "A";
	int cn = 2;
	genotype_t_ptr_List *genos = genotype_calculate_genotypes(cn,ref_base);
	mu_assert(genos != NULL, "Genotypes list was NULL");
	mu_assert(List_count(genos) == 7, "Wrong number of genotypes returned in cn 2 test.");
	cn = 4;
	genotype_t_ptr_List *sec_genos = genotype_calculate_genotypes(4,ref_base);
	mu_assert(sec_genos != NULL, "Genotypes list was NULL");
	mu_assert(List_count(sec_genos) == 13, "Wrong number of genotypes returned in cn 4 test.");
	genotype_clear_genotype_cache();
	return NULL;
}

char *test_genotype_hard_copy_genotype_t_list(){
	genotype_t_ptr_List *new_list = genotype_t_ptr_List_create(); 
	genotype_t_ptr_List *old = genotype_t_ptr_List_create();
	genotype_t *geno = genotype_init_genotype();
	int count = 10;
	char base = 'C';
	genotype_set_base_count(geno, base, count);
	mu_assert(geno->c_count==10,"Wrong number of C bases recorded");
	genotype_t_ptr_List_push(old,geno);
	genotype_hard_copy_genotype_t_list(new_list,old);
	mu_assert(List_count(new_list) == 1, "Wrong number of elements in list.");
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, new_list, first, next, cur){
		//Only one element...
		mu_assert(cur->c_count == 10, "Wrong c count in copied value.");
	}
	List_clear_destroy(genotype_t_ptr, new_list);
	List_clear_destroy(genotype_t_ptr, old);
	return NULL;
}

char *test_genotype_add_base_to_count(){
	genotype_t *geno = genotype_init_genotype();
	genotype_add_base_to_count(geno,'C');
	mu_assert(geno->c_count == 1, "Incorrect count of C after first add.");
	genotype_add_base_to_count(geno,'C');
	mu_assert(geno->c_count == 2, "Incorrect count of C after second add.");
	genotype_add_base_to_count(geno,'T');
	mu_assert(geno->t_count == 1, "Incorrect count of T after third add.");
	mu_assert(geno->c_count == 2, "Incorrect count of C after third add.");
	mu_assert(geno->a_count == 0, "Incorrect count of A after third add.");
	mu_assert(geno->g_count == 0, "Incorrect count of G after third add.");
	free(geno);
	return NULL;
}

char *test_genotype_set_base_count(){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 2;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->t_count == 2, "Incorrectly set T count.");
	free(geno);
	return NULL;
}

char *test_genotype_get_base_count(){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->a_count == 3, "Incorrectly set A count.");
	int got = genotype_get_base_count(geno,base);
	mu_assert(got==base_count,"Incorrect base count returned by function.");
	free(geno);
	return NULL;
}

char *test_genotype_init_genotype(){
	genotype_t *geno = genotype_init_genotype();
	mu_assert(geno->t_count == 0, "Incorrect count of T after init.");
	mu_assert(geno->c_count == 0, "Incorrect count of C after init.");
	mu_assert(geno->a_count == 0, "Incorrect count of A after init.");
	mu_assert(geno->g_count == 0, "Incorrect count of G after init.");
	free(geno);
	return NULL;
}

char *test_genotype_copy_genotype(){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 2;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->t_count == 2, "Incorrectly set T count.");
	genotype_t *new = genotype_copy_genotype(geno);
	mu_assert(new->a_count == 3, "Incorrectly set A count.");
	mu_assert(new->t_count == 2, "Incorrectly set T count.");
	free(new);
	free(geno);
	return NULL;
}

char *test_genotype_get_var_base_proportion(){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->a_count == 1, "Incorrectly set A count.");
	base = 'T';
	base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->t_count == 1, "Incorrectly set T count.");
	double expected = 0.5;
	double got = genotype_get_var_base_proportion(geno, base, 2);
	mu_assert(got==expected,"Incorrect variant base proportion calculated.");
	free(geno);
	geno = genotype_init_genotype();
	base = 'A';
	base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->t_count == 1, "Incorrectly set T count.");
	expected = 0.75;
	got = genotype_get_var_base_proportion(geno, base, 4);
	mu_assert(got==expected,"Incorrect variant base proportion calculated.");
	free(geno);
	return NULL;
}

char *test_genotype_get_genotype_t_as_string(){
	genotype_t *geno = genotype_init_genotype();
	char base = 'A';
	int base_count = 3;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->a_count == 3, "Incorrectly set A count.");
	base = 'T';
	base_count = 1;
	genotype_set_base_count(geno, base, base_count);
	mu_assert(geno->t_count == 1, "Incorrectly set T count.");
	char *expected = "AAAT";
	char *got = genotype_get_genotype_t_as_string(geno);
	mu_assert(strcmp(expected,got)==0,"Wrong string returned");
	free(geno);
	free(got);
	return NULL;
}

char *test_genotype_generate_genotype_list_for_cn_and_ref_base(){
	//int norm_cn, int tum_cn, char *ref_base
	int norm_cn = 2;
	int tum_cn = 2;
	//Make sure we aren't using the cache...
	char *ref_base = "A";
	genotype_store_t *store = genotype_generate_genotype_list_for_cn_and_ref_base(norm_cn, tum_cn, ref_base);
	mu_assert(store != NULL, "Genotype store not generated properly");
	mu_assert(List_count(store->normal_genos)==7,"Wrong number of normal genotypes.");
	mu_assert(List_count(store->tumour_genos)==7,"Wrong number of tumour genotypes.");

	genotype_destroy_genotype_store(store);
	genotype_clear_genotype_cache();
	
	tum_cn = 3;
	store = genotype_generate_genotype_list_for_cn_and_ref_base(norm_cn, tum_cn, ref_base);
	mu_assert(store != NULL, "Genotype store not generated properly");
	mu_assert(store->somatic_count==9,"Wrong number of somatic genotypes.");
	mu_assert(store->hom_count==3,"Wrong number of hom snp genotypes.");
	mu_assert(store->het_count==12,"Wrong number of het snp genotypes.");
	mu_assert(List_count(store->tumour_genos)==10,"Wrong number of tumour genotypes.");

	genotype_destroy_genotype_store(store);
	genotype_clear_genotype_cache();
	
	tum_cn = 4;
	store = genotype_generate_genotype_list_for_cn_and_ref_base(norm_cn, tum_cn, ref_base);
	mu_assert(store != NULL, "Genotype store not generated properly");
	mu_assert(List_count(store->normal_genos)==7,"Wrong number of normal genotypes.");
	mu_assert(List_count(store->tumour_genos)==13,"Wrong number of tumour genotypes.");
	mu_assert(store->somatic_count==12,"Wrong number of somatic genotypes.");
	
	genotype_destroy_genotype_store(store);
	genotype_clear_genotype_cache();
	
	return NULL;
}

char *all_tests() {
   mu_suite_start();
   mu_run_test(test_genotype_init_genotype); 
   mu_run_test(test_genotype_calculate_genotypes);
   mu_run_test(test_genotype_add_base_to_count);
   mu_run_test(test_genotype_hard_copy_genotype_t_list);
   mu_run_test(test_genotype_set_base_count);
   mu_run_test(test_genotype_get_base_count);
   mu_run_test(test_genotype_copy_genotype);
   mu_run_test(test_genotype_get_var_base_proportion);
   mu_run_test(test_genotype_get_genotype_t_as_string);
   mu_run_test(test_genotype_generate_genotype_list_for_cn_and_ref_base);
   return NULL;
}

RUN_TESTS(all_tests);
