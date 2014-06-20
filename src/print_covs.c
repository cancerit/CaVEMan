/**   LICENSE
* Copyright (c) 2014 Genome Research Ltd. 
* 
* Author: Cancer Genome Project cgpit@sanger.ac.uk 
* 
* This file is part of caveman_c. 
* 
* caveman_c is free software: you can redistribute it and/or modify it under 
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

#include "alg_bean.h"
#include "covs_access.h"
#include <stdio.h>
#include <assert.h>
#include <List.h>

int main(int argc, char *argv[]){
	assert(argc==3);
	char *alg_bean_loc = argv[1];
	printf("%s\n",alg_bean_loc);
	char *cov_array_loc = argv[2];
	printf("%s\n",cov_array_loc);
	FILE *alg_bean_file = fopen(alg_bean_loc,"r");
	alg_bean_t *alg = alg_bean_read_file(alg_bean_file);
	fclose(alg_bean_file);
	
	int ********arr = covs_access_read_covs_from_file(cov_array_loc,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
									
	
	cov_access_print_cov_array(arr,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
	
	
	covs_access_free_cov_array_given_dimensions(List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
				List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base),arr);
									
	alg_bean_destroy(alg);	
	return 0;
	
}