# CHANGES

## 1.15.3

* Remove requirement for CWD in config ini file

## 1.15.2

* Correct genotype bug since inception. Now uses / rather than | in sample genotype section (GT) of VCF.

## 1.15.1

* Fixed bug in ignore region calculations where ignored regions ends at the end of a chr.

## 1.15.0

* Fixed bug in fai_access test
* Updated htslib to 1.10.2
* Added missing libgen.h header

## 1.14.1

* Adding more useful error message to cn_access file reading method

## 1.14.0

* Making species and assembly required commandline arguments in the estep
* Change method by which read position index is resolved. Results in time saving.
* linasm is now required [linasm](http://linasm.sourceforge.net/index.php)
* New log and exp functions from [linasm](http://linasm.sourceforge.net/index.php)

## 1.13.16

* Correct occasional memory blowout caused by split step not excluding an ignore region.

## 1.13.15

* Add checks to sam_iter_next to ensure result is checked for errors

## 1.13.14

* Fix array length index, check loop length and error messages related to 1.13.13.

## 1.13.13

* Add AMPLICON and TARGETED to accepted protocol list. Change RNA-seq, but also keep RNA.
* https://www.ebi.ac.uk/ena/submit/reads-library-strategy

## 1.13.12

* Correct bad malloc on sscanf.

## 1.13.11

* Modify `fai_access_get_count_length_all_contigs` to use `strtok` rather than sscanf.

## 1.13.10

### Behaviour change

**Where the proper pair filter flag is used, this code now checks that the paired-end orientation is also used.**
**This will mean that mate-pair orientation (F/F or R/R) will be rejected**
**If you wish to use mate-pair data, please use previous version**

* Where a proper pair filter is used, now check for the correct paired-end orientation of F/R.
* If this is not met the read is ignored.

## 1.13.9

* Fix version check to allow equal to or greater than comparisons

## 1.13.8

* Use gzputs instead of gzprintf when writing reference contig lines to avoid buffer overflow

## 1.13.7

* Modify means by which zlib version is detected.

## 1.13.5

* Remove printf that causes all contigs to be printed to stderr.

## 1.13.4

* Fix incorrect file path introduced in 1.13.3
* Fixes #82

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
