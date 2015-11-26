### 1.7.3
* Fixed generateCavemanUMNormVCF - command line help has errors #28
* Fixed check BAM headers in setup for readgroups #34

### 1.7.2
* RG lane id search no longer requires ID to be at the beginning of the RG line.

### 1.7.1
* Added header read to bam_access_get_by_position_counts method. Fixes #41

### 1.7.0
* Updated merge script to check number of files against splitList.

### 1.6.4
* Fixed #37, Fixed bug where error thrown when no CN file, and should be using default.
* Fixed #38, removed need for stack memory usage in array reading and writing.

### 1.6.2
* Fixed bug in generateCavemanVCFUnatchedNormalPanel.c where q was not read at commandline.

### 1.6.0
* Cram support added using htslib
* Added setup script to download and compile htslib with patch (required for CRAM support).
* Added 2015 and date additional info to licenses in all files


