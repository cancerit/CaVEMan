CaVEMan
======

A C implementation of the CaVEMan program. Uses an expectation maximisation
approach to calling single base substitutions in paired data.
Designed for use with a compute farm/cluster most steps in the program make
use of an index parameter. The split step is designed to divide the genome into
chunks of adjustable size to optimise for runtime/memory usage requirements.
For simple execution of CaVEMan please see [cgpCaVEManWrapper](https://github.com/cancerit/cgpCaVEManWrapper)

| Master                                                                                                              | Dev                                                                                                              |
| ------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| [![Build Status](https://travis-ci.org/cancerit/CaVEMan.svg?branch=master)](https://travis-ci.org/cancerit/CaVEMan) | [![Build Status](https://travis-ci.org/cancerit/CaVEMan.svg?branch=dev)](https://travis-ci.org/cancerit/CaVEMan) |

Installation
============

See INSTALL.TXT

## Prerequisites

* As of version 1.6.0 CaVEMan supports cram files (with index).
* BWA Mapped, indexed, duplicate marked/removed bam files, for both a normal and tumour sample
* Reference.fasta and index
* A one based bed style format file of regions to ignore during analysis (see specified format).
* [zlib](https://zlib.net/) >= 1.2.3.5
* [linasm](http://linasm.sourceforge.net/index.php)

## Optional inputs (will result in more accurate calls)

* Normal and tumour copy number files (see specified format).
* A normal contamination of tumour value

## Processing flow

CaVEMan is executed in several distinct steps in the order listed below.

### Setup

Generates a config file for use with the remaining CaVEMan steps (it'll save you a lot of time typing commandline args).
Also generates a file named 'alg_bean' in the run directory.

	bin/caveman setup || ./bin/setupCaveman
	Usage: caveman setup -t tum.bam -n norm.bam -r reference.fa.fai -g ignore_regions.tab -e tum_cn.bed -j norm_cn.bed [-f path] [-l path] [-a path] [-wzu]

	-t  --tumour-bam [file]             Location of tumour bam
	-n  --normal-bam [file]             Location of normal bam
	-r  --reference-index [file]        Location of reference fasta index
	-g  --ignore-regions-file [file]    Location of tsv ignore regions file

	Optional
	-c  --config-file [file]            File to write caveman run config file [default:'./caveman.cfg.ini']
	-f  --results-folder [file]         Folder to write results [default:'./results']
	-l  --split-file [file]             File to write list of split sections [default:'./splitList']
	-a  --alg-bean-file [file]          Location to write alg-bean [default:'./alg_bean']
	-w  --include-smith-waterman        Include SW mapped reads in the analysis
	-z  --include-single-end            Use single end reads for this analysis
	-u  --include-duplicates            Include reads marked as duplicates in the analysis
	-e  --tumour-copy-no-file [file]    Location of tumour copy number bed file (if the extension is not .bed the file will 		be treated as 1 based start). If no copy number file is supplied then the default cn of 2 will be used
	-j  --normal-copy-no-file [file]    Location of normal copy number bed file (if the extension is not .bed the file will 		be treated as 1 based start). If no copy number file is supplied then the default cn of 2 will be used
	-h	--help                          Display this usage information.


### Split

Requires one job per entry in your reference.fasta.fai file.
Each job creates a list of segments to be analysed, these are determined by total read count
and do not include reads from ignored regions. The size of sections can be tuned using
the -m, -c and -e parameters. Once all jobs complete successfully you will need to concatenate
all split files into a single file with the name passed to the setup step in -f parameter
(or splitList if you used the default).

	bin/caveman split || ./bin/splitCaveman
	Usage: caveman split -i jobindex -f path [-c int] [-m int] [-e int]

	-f --config-file file           Path to the config file produced by setup [default: 'caveman.cfg.ini'].
	-i  --index int                 Job index (e.g. from $LSB_JOBINDEX)

	Optional
	-c  --increment int             Increment to use when deciding split sizes
	-m  --max-read-count double     Proportion of read-count to allow as a max in a split section
	-e  --read-count int            Guide for maximum read count in a section
	-h	help                        Display this usage information.

### Mstep

Requires one job per entry in the merged split file. The parameter -i referes to a line in the merged split file.
-a can be used to tune the size of section downloaded from bam at a time, this allows tuning
of the memory footprint. The mstep will create a file for each job under the results folder
(specified as a parameter in the setup step).

	bin/caveman mstep || ./bin/mstepCaveman
	Usage: caveman mstep -f config.file -i jobindex [-m int] [-a int]

	-f --config-file file                  Path to the config file produced by setup [default: 'caveman.cfg.ini'].
	-i --index int                         Job index.

	Optional
	-a --split_size int                    Size of section to retrieve at a time from bam file. Allows memory footprint tuning [default:50000].
	-m  --min-base-qual int                Minimum base quality to include in analysis [default:11]
	-h	help                                Display this usage information.

### Merge

Runs as a single job, merging all the cov_array files generated by the mstep into a file
representing the profile of the whole genome. The resulting file is named cov_array and stored
in the root run folder. Another file named prob_array is also created.

	bin/caveman merge || ./bin/mergeCaveman
	Usage: caveman merge -f config_file [-c path] [-p path]

	-f --config-file file                Path to the config file produced by setup. [default: 'caveman.cfg.ini']

	Optional
	-c  --covariate-file filename        Location to write merged covariate array [default: covs_arr]
	-p  --probabilities-file filename    Location to write probability array [default: probs_arr]
	-h	help                              Display this usage information.

### Estep

The final step in calling variants using CaVEMan. As was the case with the mstep, a job per
entry in the merged split list is required. Copy number (-e, -j), and normal contamination (-k)
are required (See default settings for advice on obtaining results if you don't have these).
Each job will create several files named according to the corresponding line in the merged split file.
*_muts.vcf and *_snps.vcf are the calls above the given cutoffs for somatic and SNP probabilities alike.
The *.no_analysis.bed file lists sections not analysed by caveman, due to being ignored or lacking coverage.
Using the debugoption will output another file named *.dbg.vcf containing  a line for every position
in the genome that was analysed with read counts (as seen by CaVEMan) and top two probabilities
calculated.

	bin/caveman estep || ./bin/estepCaveman
	Usage: caveman estep -i jobindex [-f file] [-m int] [-k float] [-b float] [-p float] [-q float] [-x int] [-y int] [-c 	float] [-d float] [-a int]

	-i  --index [int]                                Job index (e.g. from $LSB_JOBINDEX)

	Optional
	-f  --config-file [file]                         Path to the config file produced by setup. [default:'caveman.cfg.ini']
	-m  --min-base-qual [int]                        Minimum base quality for inclusion of a read position [default:11]
	-c  --prior-mut-probability [float]              Prior somatic probability [default:0.000006]
	-d  --prior-snp-probability [float]              Prior germline mutant probability [default:0.000100]
	-k  --normal-contamination [float]               Normal contamination of tumour [default:0.100000]
	-b  --reference-bias [float]                     Reference bias [default:0.950000]
	-p  --mut-probability-cutoff [float]             Minimum probability call for a somatic mutant position to be output 		[default:0.800000]
	-q  --snp-probability-cutoff [float]             Minimum probability call for a germline mutant position to be output 		[default:0.950000]
	-x  --min-tum-coverage [int]                     Minimum tumour coverage for analysis of a position [default:1]
	-y  --min-norm-coverage [int]                    Minimum normal coverage for analysis of a position [default:1]
	-a  --split-size [int]                           Size of section to retrieve at a time from bam file. Allows memory 								 footprint tuning [default:50000].
	-s  --debug                                      Adds an extra output to a debug file. Every base analysed has an 			output
	-g  --cov-file [file]                            File location of the covariate array. [default:'covs_arr']
	-o  --prob-file [file]                           File location of the prob array. [default:'probs_arr']
	-v  --species-assembly [string]                  Species assembly (eg 37/GRCh37), required if bam header SQ lines do 								 not contain AS and SP information.
	-w  --species [string]                           Species name (eg Human), required if bam header SQ lines do not 								 contain AS and SP information.
	-n  --normal-copy-number [int]                   Copy number to use when filling gaps in the normal copy number file 								 [default:2].
	-t  --tumour-copy-number [int]                   Copy number to use when filling gaps in the tumour copy number file 								 [default:2].
	-l  --normal-protocol [string]                   Normal protocol. Ideally this should match -r but not checked 									 (WGS|WGX|RNA) [default:WGS].
	-r  --tumour-protocol [string]                   Tumour protocol. Ideally this should match -l but not checked 									 (WGS|WGX|RNA) [default:WGS].
	-P  --normal-platform [string]                   Normal platform. Overrides the values retrieved from bam header.
	-T  --tumour-platform [string]                   Tumour platform. Overrides the values retrieved from bam header.
	-M  --max-copy-number [int]                      Maximum copy number permitted. If exceeded the copy number for the 								 offending region will be set to this value. [default:10].
	-h  --help                                     Display this usage information.


## File Formats


### Copy Number

Copy number files are taken in a 1-based bed style tab separated format, or BED format if the suffix is .bed.
Where the columns are chromosome,start,stop,copynumber(integer). A separate file is required
for normal and tumour. Each file should have a copy number assigned for every region requested
to be analysed (NB, CaVEMan set CN to 2 in regions where copy number is 0).

	Example:
		1	0	20000	2
		2	0	2500	4
		2	2501	500000	6

### Ignored regions file

A 1-based bed style tab separated format file of regions to be ignored during analysis.
An example might be regions known to have extreme depth of mapped reads through mismapping.


## Default Settings

There may be cases where copy number is not available. In order to extract the best results
it is advised to use copy number 2 in the normal and 5 in the tumour in combination with
a normal contamination of 0.1 . This gives CaVEMan a broad range over which variants will
be called compared to a copy number of 2 in the tumour.

LICENCE
=======

Copyright (c) 2014-2018 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

This file is part of CaVEMan.

CaVEMan is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
