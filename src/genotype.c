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

#include <genotype.h>
#include <dbg.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define ELEMENT_TYPE genotype_t_ptr
#define ELEMENTS_PER_NODE 4
#include <List.c>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

#define ELEMENT_TYPE combined_genotype_t_ptr
#define ELEMENTS_PER_NODE 8
#include <List.c>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

#define ELEMENT_TYPE combined_genotype_t_ptr_List_ptr
#define ELEMENTS_PER_NODE 8
#include <List.c>
#undef ELEMENT_TYPE
#undef ELEMENTS_PER_NODE

char bases[4] = {'A','C','G','T'};
char *bases_str[4] = {"A","C","G","T"};

//Making a huge assumption that copy number will not go above 250 here....
genotype_t_ptr_List *geno_cache[250][4] = {};

genotype_t_ptr_List *genotype_hard_copy_genotype_t_list(genotype_t_ptr_List *new_list, genotype_t_ptr_List *old) {
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, old, first, next, cur) {
		genotype_t *tmp = genotype_copy_genotype(cur);
		genotype_t_ptr_List_push(new_list,tmp);
	}
	return new_list;
}

void genotype_add_base_to_count(genotype_t *geno, const char base){
	assert(geno != NULL);
	switch(base){
		case 'A':
			geno->a_count++;
			break;

		case 'C':
			geno->c_count++;
			break;

		case 'G':
			geno->g_count++;
			break;

		case 'T':
			geno->t_count++;
			break;

		default:
			sentinel("Incorrect base passed to add: %c",base);
			break;
	};
	return;
error:
	return;
}

void genotype_set_base_count(genotype_t *geno, const char base, int count){
	assert(geno != NULL);
	switch(base){
		case 'A':
			geno->a_count = count;
			break;

		case 'C':
			geno->c_count = count;
			break;

		case 'G':
			geno->g_count = count;
			break;

		case 'T':
			geno->t_count = count;
			break;

		default:
			break;
	};
	return;
}

int genotype_get_base_count(genotype_t *geno, const char base){
	assert(geno != NULL);
	switch(base){
		case 'A':
			return geno->a_count;

		case 'C':
			return geno->c_count;

		case 'G':
			return geno->g_count;

		case 'T':
			return geno->t_count;

		default:
			return -1;
	};
}

genotype_t *genotype_copy_genotype(genotype_t *geno){
	assert(geno != NULL);
	genotype_t *new_gen = genotype_init_genotype();
	check(new_gen != NULL,"Error trying generating genotype_t");
   check_mem(new_gen);
   int i=0;
   for(i=0; i<4; i++){
   	genotype_set_base_count(new_gen,bases[i],genotype_get_base_count(geno,bases[i]));
   }
   return new_gen;
error:
	if(new_gen) free(new_gen);
	return NULL;
}

void genotype_init_genotype_t(genotype_t *gen){
	int i=0;
   for(i=0; i<4; i++){
   	genotype_set_base_count(gen, bases[i], 0);
   }
   gen->var_base = '\0';
   gen->var_base_prop = 0.0;
   return;
}

genotype_t *genotype_init_genotype(){
	genotype_t *geno = malloc(sizeof(struct genotype_t));
	check_mem(geno);
	genotype_init_genotype_t(geno);
	return geno;
error:
	if (geno) free(geno);
	return NULL;
}

genotype_t_ptr_List *genotype_generate_unique_genotype_list(genotype_t_ptr_List *li) {
	int geno_size = List_count(li);
	genotype_t **genos;
	genos = malloc(sizeof(genotype_t *) * geno_size);
	check_mem(genos);
	memset(genos,0,(sizeof(genotype_t *) * geno_size));
	genotype_t_ptr_List *uniq = genotype_t_ptr_List_create();
	int i=0;
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, li,first,next,cur) {
		int j=0;
		int match = 0;
		for(j=0;j<geno_size;j++){
			if(genos[j] > 0 && genotype_equals(genos[j],cur)==1){
				match = 1;
			}
		}
		if(match == 0){
			//genos[j] = malloc(sizeof(struct genotype_t));
			//check_mem(genos[i]);
			genos[i] = cur;
			genotype_t_ptr_List_push(uniq,genos[i]);
		}else{
			free(cur);
		}
		i++;
	}

	free(genos);
	genotype_t_ptr_List_destroy(li);
	return uniq;
 error:
	List_clear_destroy(genotype_t_ptr, li);
	if(genos) free(genos);
	return NULL;
}

genotype_t_ptr_List *genotype_calculate_genotypes(int copy_num, char *ref_base) {
	//Use a store of copy number then reference base to get the list (if it doesn't exist we create it).
	//Use the copy number to store the location in the cn-cache array. if it's not there or NULL we default to the creation.
	int bs=0;
	genotype_t_ptr_List *unique = NULL;

	if(geno_cache[copy_num] != NULL){
		//This size will always be 4 and ordered by ref_base;
		for (bs=0; bs<4; bs++){
			if(strcmp(ref_base,bases_str[bs]) == 0){
				if(geno_cache[copy_num][bs] != NULL){
					return geno_cache[copy_num][bs];
				}
				//Break so bs is carried down and we don't need to work it out again!.
				break;
			}
		}
	}

	//If this entry doesn't exist yet we create it.
	genotype_t_ptr_List *genos = genotype_t_ptr_List_create();
	int i=0;
	for(i=0; i<copy_num; i++){
		// If genotypes have already been made, we need to append to them.
		if(List_count(genos)>0){
			genotype_t_ptr_List *tmp = genotype_t_ptr_List_create();
			//We need four copies of each genotype to append one of each bean to.
			LIST_FOR_EACH_ELEMENT(genotype_t_ptr, genos, first, next, cur) {
				genotype_t *gen = cur;
				//Create 2 copies, one will be reference, the next will be mutant base.
				genotype_t *gen1 = genotype_copy_genotype(gen);
				genotype_t *gen2 = genotype_copy_genotype(gen1);
				genotype_add_base_to_count(gen2,ref_base[0]);
				genotype_t_ptr_List_push(tmp,gen2);
				int j=0;
				for(j=0; j<4; j++){
					int base_count = genotype_get_base_count(gen1,bases[j]);
					check(base_count >= 0, "Error fetching base count for base: %c with genotype: %s.",bases[j],genotype_get_genotype_t_as_string(gen1));
					if(base_count>=1 && ref_base[0] != bases[j]){
						genotype_add_base_to_count(gen1,bases[j]);
						genotype_t_ptr_List_push(tmp,gen1);
					}
				}
				if(genotype_equals(List_last(tmp),gen1)!=1) free(gen1);
			}
			List_clear_destroy(genotype_t_ptr, genos);
			genos = tmp;
		}else{
			//Create a genotype for each possible base.
			int j=0;
			for(j=0; j<4; j++){
				genotype_t *gen = genotype_init_genotype();
				check(gen != NULL, "Error generating genotype.");
				genotype_add_base_to_count(gen,bases[j]);
				genotype_t_ptr_List_push(genos,gen);
			}
		}
	}

	unique = genotype_generate_unique_genotype_list(genos);

	//Now we have a list of genotypes we calculate variant base and variant base proportion.
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, unique, first, next, gen) {
		//Variant base
		char var_b = genotype_get_var_base(gen,ref_base[0]);
		check(var_b,"Error retrieving variant base for genotype: %s given ref base: %s.",genotype_get_genotype_t_as_string(gen),ref_base);
		gen->var_base = var_b;
		//variant base proportion
		gen->var_base_prop = genotype_get_var_base_proportion(gen, ref_base[0], copy_num);
		check(gen->var_base_prop >= 0, "Error calculating variant base proportion for genotype: %s given ref base: %s",genotype_get_genotype_t_as_string(gen),ref_base);
		int i=0;
		for(i=0;i<4;i++){
			if(bases[i] == var_b){
				gen->var_base_idx = i;
			}
		}
	}

	geno_cache[copy_num][bs] = unique;
	return unique;

 error:
	if (unique) List_clear_destroy(genotype_t_ptr, unique);
	if (genos) List_clear_destroy(genotype_t_ptr, genos);
	return NULL;
}

long double genotype_get_var_base_proportion(genotype_t *gen, const char ref_base, int copy_num) {
	if(genotype_get_base_count(gen,ref_base) == copy_num){
		return (long double)0;
	}
	int i=0;
	for(i=0; i<4; i++){
		if(ref_base != bases[i] && genotype_get_base_count(gen,bases[i])){
			long double prop = (long double)genotype_get_base_count(gen,bases[i]) / (long double)copy_num;
			return prop;
		}
	}
	return (long double)-1;
}

char *genotype_get_genotype_t_as_string(genotype_t *geno) {
	char *gen_str = malloc(
			       (sizeof(char) * genotype_get_base_count(geno,'A')) +
			       (sizeof(char) * genotype_get_base_count(geno,'C')) +
			       (sizeof(char) * genotype_get_base_count(geno,'G')) +
			       (sizeof(char) * genotype_get_base_count(geno,'T')) +
			       + sizeof(char) * 5);
	check_mem(gen_str);
	gen_str = strcpy(gen_str,"");
	int i=0;
	for(i=0; i<4; i++){
		int count = genotype_get_base_count(geno,bases[i]);
		int app = 0;
		while(app<count){
			strcat(gen_str,bases_str[i]);
			app++;
		}
	}
	return gen_str;
 error:
	if(gen_str) free(gen_str);
	return NULL;
}

char genotype_get_var_base(genotype_t *geno, const char ref_base){
	int i=0;
	for(i=0; i<4; i++){
		if(ref_base != bases[i] && genotype_get_base_count(geno,bases[i]) > 0){
			return bases[i];
		}
	}
	return ref_base;
}

combined_genotype_t_ptr_List *genotype_create_combined_List(genotype_t_ptr_List *li, genotype_t *norm) {
	combined_genotype_t_ptr_List *out = combined_genotype_t_ptr_List_create();
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, li,first,next,cur) {
		combined_genotype_t *combo = malloc(sizeof(combined_genotype_t));
		check_mem(combo);
		combo->norm_geno = norm;
		combo->tum_geno = cur;
		combo->prob = 0.0;
		combined_genotype_t_ptr_List_push(out,combo);
	}
	return out;
 error:
	if (out) List_clear_destroy(combined_genotype_t_ptr, out);
	return NULL;
}

genotype_t *genotype_find_ref_genotype(genotype_t_ptr_List *tum_genos,int tum_cn, const char ref_base) {
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, tum_genos, first, next, tum) {
		if(genotype_get_base_count(tum, ref_base) == tum_cn){
			return tum;
		}
	}
	return NULL;
}

combined_genotype_t_ptr_List *genotype_calculate_related_genotypes(genotype_t *norm,genotype_t_ptr_List *tum_genos,int *count, int tum_cn, const char ref_base) {
	assert(norm != NULL);
	assert(tum_genos != NULL);
	//Iterate through each tumour genotype and see if it can fit with the normal genotype passed. If it can, add it to the list.
	genotype_t_ptr_List *store = genotype_t_ptr_List_create();
	//Lastly make the list into an array.
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, tum_genos, first, next, tum) {
		//The variant bases should match (unless hom ref in tumour) and normal genotype is not reference - SNPs.
		if(norm->var_base != ref_base){
			//Hom SNPs
			if(norm->var_base == tum->var_base && norm->var_base_prop == 1.0 && tum->var_base_prop == 1.0){
				genotype_t_ptr_List_push(store,tum);
				//HET SNPs
			}else if((tum->var_base == norm->var_base || tum->var_base == ref_base) && norm->var_base_prop != 1.0){
				genotype_t_ptr_List_push(store,tum);
			}
		}else if(norm->var_base == ref_base && tum->var_base != ref_base){
			//Reference normal genotype, so all these are somatics so ignore the reference tumour genotype.
			genotype_t_ptr_List_push(store,tum);
		}
	}
	combined_genotype_t_ptr_List *list = genotype_create_combined_List(store, norm);
	*count += List_count(list);
	genotype_t_ptr_List_destroy(store);
	return list;
}

int genotype_get_size_of_list_of_lists(combined_genotype_t_ptr_List_ptr_List *li) {
	int count = 0;
	LIST_FOR_EACH_ELEMENT(combined_genotype_t_ptr_List_ptr, li,first,next,this) {
		count += List_count(this);
	}
	return count;
}

void genotype_put_genotype_combos_into_array(combined_genotype_t **combos,combined_genotype_t_ptr_List_ptr_List *list) {
	//int size = genotype_get_size_of_list_of_lists(list);
	int i=0;
	LIST_FOR_EACH_ELEMENT(combined_genotype_t_ptr_List_ptr, list,first,next,this) {
		LIST_FOR_EACH_ELEMENT(combined_genotype_t_ptr, this, first, next, here){
			combos[i] = here;
			i++;
		}
		combined_genotype_t_ptr_List_destroy(this);
	}
	return;
}

void genotype_fill_het_snp_norms_list(genotype_t_ptr_List *het_snps_norm,combined_genotype_t **het_snp_norm_genotypes) {
	int i=0;
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, het_snps_norm, first, next, this) {
		combined_genotype_t *g = malloc(sizeof(struct combined_genotype_t));
		check_mem(g);
		g->norm_geno = this;
		g->prob = 0.0;
		het_snp_norm_genotypes[i] = g;
		i++;
	}
	return;
 error:
	return;
}

void genotype_set_snp_and_mut_genotypes(genotype_t_ptr_List *norm_genos, genotype_t_ptr_List *tum_genos, int norm_cn, int tum_cn, char *ref_base, genotype_store_t *store) {
	int somatic_count = 0;
	int het_count = 0;
	int hom_count = 0;

	combined_genotype_t **somatic_combos = NULL;
	combined_genotype_t **het_snp_combos = NULL;
	combined_genotype_t **hom_snp_combos = NULL;

	combined_genotype_t_ptr_List_ptr_List *hom_snps = combined_genotype_t_ptr_List_ptr_List_create();
	combined_genotype_t_ptr_List_ptr_List *het_snps = combined_genotype_t_ptr_List_ptr_List_create();
	combined_genotype_t_ptr_List_ptr_List *somatics = combined_genotype_t_ptr_List_ptr_List_create();
	genotype_t_ptr_List *het_snps_norm = genotype_t_ptr_List_create();
	combined_genotype_t *ref_genotype = malloc(sizeof(combined_genotype_t));
	//Iterate through normal genotypes.
	LIST_FOR_EACH_ELEMENT(genotype_t_ptr, norm_genos, first, next, norm) {
		int ref_count = genotype_get_base_count(norm, ref_base[0]);
		if(ref_count == norm_cn){
			ref_genotype->norm_geno = norm;
			ref_genotype->tum_geno = genotype_find_ref_genotype(tum_genos,tum_cn,ref_base[0]);
			ref_genotype->prob = 0.0;
			//somatic genotypes
			combined_genotype_t_ptr_List *list = genotype_calculate_related_genotypes(norm,tum_genos,&somatic_count,tum_cn,ref_base[0]);
			combined_genotype_t_ptr_List_ptr_List_push(somatics,list);
		}else if(ref_count > 0){
			//HET SNP normal genotype
			genotype_t_ptr_List_push(het_snps_norm,norm);
			combined_genotype_t_ptr_List *list = genotype_calculate_related_genotypes(norm,tum_genos,&het_count,tum_cn,ref_base[0]);
			combined_genotype_t_ptr_List_ptr_List_push(het_snps,list);
		}else{
			//HOM snp genotype
			combined_genotype_t_ptr_List *list = genotype_calculate_related_genotypes(norm,tum_genos,&hom_count,tum_cn,ref_base[0]);
			combined_genotype_t_ptr_List_ptr_List_push(hom_snps,list);
		}
	}
	combined_genotype_t **het_snp_norm_genotypes = malloc(sizeof(combined_genotype_t *) * List_count(het_snps_norm));
	check_mem(het_snp_norm_genotypes);
	int het_norm_count = List_count(het_snps_norm);
	genotype_fill_het_snp_norms_list(het_snps_norm,het_snp_norm_genotypes);
	somatic_combos = malloc(sizeof(combined_genotype_t *) * somatic_count);
	check_mem(somatic_combos);
	memset(somatic_combos,0,(sizeof(combined_genotype_t *) * somatic_count));
	het_snp_combos = malloc(sizeof(combined_genotype_t *) * het_count);
	check_mem(het_snp_combos);
	memset(het_snp_combos,0,(sizeof(combined_genotype_t *) * het_count));
	hom_snp_combos = malloc(sizeof(combined_genotype_t *) * hom_count);
	check_mem(hom_snp_combos);
	memset(hom_snp_combos,0,(sizeof(combined_genotype_t *) * hom_count));

	//Now assign everything in the lists to the appropriate pointer arrays.
	genotype_put_genotype_combos_into_array(somatic_combos,somatics);
	genotype_put_genotype_combos_into_array(het_snp_combos,het_snps);
	genotype_put_genotype_combos_into_array(hom_snp_combos,hom_snps);

	//and put them in the store
	store->ref_genotype = ref_genotype;
	store->het_snp_genotypes = het_snp_combos;
	store->het_count = het_count;
	store->hom_snp_genotypes = hom_snp_combos;
	store->hom_count = hom_count;
	store->somatic_genotypes = somatic_combos;
	store->somatic_count = somatic_count;
	store->het_snp_norm_genotypes = het_snp_norm_genotypes;
	store->het_norm_count = het_norm_count;
	store->ref_geno_norm_prob = 0;
	store->ref_geno_tum_prob = 0;
	store->tum_max = somatic_count;
	store->norm_max = het_norm_count;
	if(het_count>store->tum_max){
		store->tum_max = het_count;

	}
	if(hom_count>store->tum_max){
		store->tum_max = hom_count;
	}
	store->total_max = store->tum_max;
	if(hom_count>store->norm_max){
		store->norm_max = hom_count;
	}

	if(store->norm_max>store->total_max){
		store->total_max = store->norm_max;
	}

	//Free up the lists.
	combined_genotype_t_ptr_List_ptr_List_destroy(hom_snps);
	combined_genotype_t_ptr_List_ptr_List_destroy(het_snps);
	combined_genotype_t_ptr_List_ptr_List_destroy(somatics);
	genotype_t_ptr_List_destroy(het_snps_norm);
	return;
 error:
	combined_genotype_t_ptr_List_ptr_List_destroy(hom_snps);
	combined_genotype_t_ptr_List_ptr_List_destroy(het_snps);
	combined_genotype_t_ptr_List_ptr_List_destroy(somatics);
	genotype_t_ptr_List_destroy(het_snps_norm);
	if(somatic_combos) free(somatic_combos);
	if(het_snp_combos) free(het_snp_combos);
	if(hom_snp_combos) free(hom_snp_combos);
	if(het_snp_norm_genotypes) free(het_snp_norm_genotypes);
	if(ref_genotype) free(ref_genotype);
	return;
}

int genotype_equals(genotype_t *g_a,genotype_t *g_b) {
	if(g_a->a_count != g_b->a_count){
		return 0;
	}
	if(g_a->c_count != g_b->c_count){
		return 0;
	}
	if(g_a->g_count != g_b->g_count){
		return 0;
	}
	if(g_a->t_count != g_b->t_count){
		return 0;
	}
	return 1;
}

int genotype_get_total_base_count(genotype_t *geno){
	int count = 0;
	count += geno->a_count;
	count += geno->c_count;
	count += geno->g_count;
	count += geno->t_count;
	return count;
}

genotype_store_t *genotype_generate_genotype_list_for_cn_and_ref_base(int norm_cn, int tum_cn, char *ref_base){
	genotype_store_t *store = malloc(sizeof(struct genotype_store_t));
	check_mem(store);
	store->normal_genos = genotype_calculate_genotypes(norm_cn, ref_base);
	check(store->normal_genos != NULL,"Error calculating normal genotypes by copy number.");
	store->tumour_genos = genotype_calculate_genotypes(tum_cn, ref_base);
	check(store->tumour_genos != NULL,"Error calculating tumour genotypes by copy number.");
	//Means to store probabilities of genotypes, linking them together via the combined genotype;
	//During this calculation we set the number of SNP and number of mutant genotypes fields as well
	genotype_set_snp_and_mut_genotypes(store->normal_genos, store->tumour_genos, norm_cn, tum_cn, ref_base, store);

	return store;
error:
	genotype_destroy_genotype_store(store);
	return NULL;
}

void genotype_destroy_genotype_store(genotype_store_t *store){
	if(store) {
		if(store->ref_genotype) free(store->ref_genotype);
		if(store->het_snp_genotypes){
			int i=0;
			for(i=0;i<store->het_count;i++){
				if(store->het_snp_genotypes[i]) free(store->het_snp_genotypes[i]);
			}
			free(store->het_snp_genotypes);
		}
		if(store->hom_snp_genotypes){
			int i=0;
			for(i=0;i<store->hom_count;i++){
				if(store->hom_snp_genotypes[i]) free(store->hom_snp_genotypes[i]);
			}
			free(store->hom_snp_genotypes);
		}
		if(store->somatic_genotypes){
			int i=0;
			for(i=0;i<store->somatic_count;i++){
				if(store->somatic_genotypes[i]) free(store->somatic_genotypes[i]);
			}
			free(store->somatic_genotypes);
		}
		if(store->het_snp_norm_genotypes){
			int i=0;
			for(i=0;i<store->het_norm_count;i++){
				if(store->het_snp_norm_genotypes[i]) free(store->het_snp_norm_genotypes[i]);
			}
			free(store->het_snp_norm_genotypes);
		}
		free(store);
	}
	return;
}

void genotype_clear_genotype_cache(){
	int i=0;
	for(i=0;i<250;i++){
		int j=0;
		for(j=0;j<4;j++){
			if(geno_cache[i][j] != NULL && List_count(geno_cache[i][j]) > 0 ){
			  List_clear_destroy(genotype_t_ptr, geno_cache[i][j]);
			  geno_cache[i][j] = NULL;
			}
		}
	}
}
