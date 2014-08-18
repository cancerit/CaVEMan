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

#ifndef _config_file_access_h
#define _config_file_access_h

#include <stdio.h>

int config_file_access_read_config_file(FILE *file, char *tum_bam_file, char *norm_bam_file, char *ref_idx,
			char *ignore_regions_file, char *alg_bean_loc, char *results, char *list_loc, int *includeSW,
			int *includeSingleEnd, int *includeDups, char *version);

int config_file_access_write_config_file(FILE *file, char *tum_bam_file, char *norm_bam_file, char *ref_idx,
			char *ignore_regions_file, char *alg_bean_loc, char *results, char *list_loc, int includeSW,
			int includeSingleEnd, int includeDups);

#endif
