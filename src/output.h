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
*/

#ifndef _output_h
#define _output_h

#include <algos.h>

int output_vcf_variant_position(estep_position_t *pos, FILE *out, char *chrom);
int output_vcf_header(FILE *out, char *tum_bam, char *norm_bam, char *ref_seq_loc,
													char *assembly, char *species, char *norm_prot, char *tum_prot,
													char *norm_plat, char *tum_plat);
char *output_generate_info_lines();
char *output_generate_format_lines();
char *output_generate_reference_contig_lines(char *bam_file, char *assembly, char *species);
int output_append_position_to_no_analysis(char *chr_name, int start_one_base, int stop);
int output_flush_no_analysis(char *chr_name);
void output_set_no_analysis_file(FILE *file);
void output_set_no_analysis_section_list(List *sections);

#endif
