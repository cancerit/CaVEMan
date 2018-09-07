# CHANGES

## 1.13.3

* Add gzbuffer call after gzopen to ensure we don't hit the limit where many contigs are printed
* Added unit test for additional methods
* Fixes #77
* Update license dates

## 1.13.2

* Futher fixes for `sentinel` calls.
* Correct version of REQUIRED_MIN_LIBZ as we are looking for greater than.

## 1.13.1

* Ensure all `sentinel` calls return non-zero exit code.

## 1.13.0

* Overlapping reads now handled #78
* Fixes for #64, #65, #61

## 1.12.1

* Resolves #74

## 1.12.0

* Now outputs snp and mut vcf files as gzipped.
* Requires [zlib](https://zlib.net/) >= 1.2.3.5

## 1.11.3

* Makefile update to cope with central/prefix install of htslib

## 1.11.2

* Update htslib version (1.11.1 updated an unused reference to this archive)

## 1.11.0

* Rearrangement to htslib pileup code gives large speed increases

## 1.10.1

* Correction to fclose checking in estep.c

## 1.10.0

* Added checks to all fflush and fclose calls

## 1.9.5

* Removed dependency on ENA during compilation

## 1.9.4

* Correct main method to ensure failure is passed through from running the by section main methods.

## 1.9.3

* Fix bug in ignore regions file reading where bed file resulted in incorrect coords

## 1.9.2

* Resolves #45 - New [samtools/htslib 1.3](https://github.com/samtools/htslib/releases/tag/1.3) to remove need for patch.

## 1.9.1

* Removed unnecessary dependancy on `rsync` when `cp` will do.

## 1.9.0

* Corrections to install methods to ensure all relevant files were installed.

## 1.8.0

* Commmandline params checked
* Readlength taken into account in split sections
* setup.sh installs perl scripts

## 1.7.3

* Fixed generateCavemanUMNormVCF - command line help has errors #28
* Fixed check BAM headers in setup for readgroups #34

## 1.7.2

* RG lane id search no longer requires ID to be at the beginning of the RG line.

## 1.7.1

* Added header read to bam_access_get_by_position_counts method. Fixes #41

## 1.7.0

* Updated merge script to check number of files against splitList.

## 1.6.4

* Fixed #37, Fixed bug where error thrown when no CN file, and should be using default.
* Fixed #38, removed need for stack memory usage in array reading and writing.

## 1.6.2

* Fixed bug in generateCavemanVCFUnatchedNormalPanel.c where q was not read at commandline.

## 1.6.0

* Cram support added using htslib
* Added setup script to download and compile htslib with patch (required for CRAM support).
* Added 2015 and date additional info to licenses in all files
