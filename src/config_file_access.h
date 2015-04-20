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

#ifndef _config_file_access_h
#define _config_file_access_h

#include <stdio.h>

int config_file_access_read_config_file(FILE *file, char *tum_bam_file, char *norm_bam_file, char *ref_idx,
			char *ignore_regions_file, char *alg_bean_loc, char *results, char *list_loc, int *includeSW,
			int *includeSingleEnd, int *includeDups, char *version, char *norm_cn, char *tum_cn);

int config_file_access_write_config_file(FILE *file, char *tum_bam_file, char *norm_bam_file, char *ref_idx,
			char *ignore_regions_file, char *alg_bean_loc, char *results, char *list_loc, int includeSW,
			int includeSingleEnd, int includeDups, char *norm_cn, char *tum_cn);

#endif
