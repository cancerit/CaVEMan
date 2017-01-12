cwlVersion: v1.0
class: CommandLineTool
baseCommand: caveman setup
inputs:
  tumloc:
    type: File?
    inputBinding:
      prefix: --tumour-bam
      separate: true
      position: 1
  normloc:
    type: File?
    inputBinding:
      prefix: --normal-bam
      separate: true
      position: 2
  referenceIndex:
    type: File?
    inputBinding:
      prefix: --reference-index
      separate: true
      position: 3
  ignoreRegions:
    type: File?
    inputBinding:
      prefix: --ignore-regions-file
      separate: true
      position: 4
  configFile:
    type: File?
    inputBinding:
      prefix: --config-file
      separate: true
      position: 5
  resultsFolder:
    type: File?
    inputBinding:
      prefix: --results-folder
      separate: true
      position: 6
  splitFile:
    type: File?
    inputBinding:
      prefix: --split-file
      separate: true
      position: 7
  algBeanFile:
    type: File?
    inputBinding:
      prefix: --alg-bean-file
      separate: true
      position: 8
  include-smith-waterman:


-w  --include-smith-waterman        Include SW mapped reads in the analysis
-z  --include-single-end            Use single end reads for this analysis
-u  --include-duplicates            Include reads marked as duplicates in the analysis
-e  --tumour-copy-no-file [file]    Location of tumour copy number bed file (if the extension is not .bed the file will be treated as 1 based start). If no copy number file is supplied then the default cn of 2 will be used
-j  --normal-copy-no-file [file]    Location of normal copy number bed file (if the extension is not .bed the file will be treated as 1 based start). If no copy number file is supplied then the default cn of 2 will be used