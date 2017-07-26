/**   LICENSE
* Copyright (c) 2014-2015 Genome Research Ltd.
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

#ifndef _covs_access_h
#define _covs_access_h

#include <stdint.h>
#include "zlib.h"


uint64_t ********covs_access_generate_cov_array_given_dimensions(int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
void covs_access_free_cov_array_given_dimensions(int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size, uint64_t ********arr);
int covs_access_write_covs_to_file(char *file_loc,uint64_t ********arr,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
uint64_t ********covs_access_read_covs_from_file(char *file_loc,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
int cov_access_compare_two_cov_arrays(uint64_t ********first,uint64_t ********second,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
int cov_access_compare_two_prob_arrays(long double ********first,long double ********second,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
void cov_access_print_cov_array(uint64_t ********arr,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
void covs_access_merge_count_arrays(uint64_t ********arr_1, uint64_t ********ar_2, int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
long double ********covs_access_generate_probability_array(uint64_t ********arr,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
void cov_access_print_prob_array(long double ********arr,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
void covs_access_free_prob_array_given_dimensions(int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size, long double ********arr);
int covs_access_write_probs_to_file(char *file_loc,long double ********arr,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
long double ********covs_access_read_probs_from_file(char *file_loc,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);
void cov_access_print_cov_and_prob_array(uint64_t ********arr,long double ********arr_2,int rd_order_size,int strand_size, int lane_size, int rpos_size, int mq_size, int bq_size, int ref_b_size, int called_b_size);

#endif
