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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_alg_bean_tests.h"
#include "check_algos_tests.h"
#include "check_bam_access_tests.h"
#include "check_cn_access_tests.h"
#include "check_config_file_access_tests.h"
#include "check_covs_access_tests.h"
#include "check_fai_access_tests.h"
#include "check_genotype_tests.h"
#include "check_ign_region_access_tests.h"
#include "check_output_tests.h"
#include "check_split_access_tests.h"

int main (void)
{
    int number_failed;
    Suite *s = check_alg_bean_tests_suite ();
    SRunner *sr = srunner_create (s);
    srunner_add_suite (sr, check_algos_tests_suite());
    srunner_add_suite (sr, check_bam_access_tests_suite());
    srunner_add_suite (sr, check_cn_access_tests_suite());
    srunner_add_suite (sr, check_config_file_access_tests_suite());
    srunner_add_suite (sr, check_covs_access_tests_suite());
    srunner_add_suite (sr, check_fai_access_tests_suite());
    srunner_add_suite (sr, check_genotype_tests_suite());
    srunner_add_suite (sr, check_ign_region_access_tests_suite());
    srunner_add_suite (sr, check_output_tests_suite());
    srunner_add_suite (sr, check_split_access_tests_suite());

    srunner_run_all (sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed (sr);
    srunner_free (sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
