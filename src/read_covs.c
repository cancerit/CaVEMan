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

#include <covs_access.h>
#include <alg_bean.h>
#include <stdio.h>


int main(int argc, char *argv[]){

FILE *alg_file = fopen(argv[1],"r");
alg_bean_t *alg = alg_bean_read_file(alg_file);

//Create an empty covariate array.
int ********covs = covs_access_read_covs_from_file(argv[2],List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));

long double ********probs = covs_access_read_probs_from_file(argv[3],List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));
							
							
cov_access_print_cov_and_prob_array(covs,probs,List_count(alg->read_order),List_count(alg->strand),List_count(alg->lane),
							List_count(alg->rd_pos),List_count(alg->map_qual),List_count(alg->base_qual),List_count(alg->ref_base),List_count(alg->call_base));

return 1;
}