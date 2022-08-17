def get_chr_count (genome, ignorecontigs) {
  echo "$genome\t$ignorecontigs"
}


process setup {

  input:
    //setup
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    file ('ign.file')

    //Optional setup - or empty file and blanket CN
    val normcn
    val tumcn
    //tumcn
    // file ('caveman.cfg.ini')
    // file ('alg_bean')
    // path ('results')
    // file ('splitList')
    val includeSW
    val includeSE
    val includeDups

  publishDir "$baseDir", mode: 'copy', pattern: 'caveman.cfg.ini', overwrite: true
  publishDir "$baseDir", mode: 'copy', pattern: 'alg_bean', overwrite: true

  output:
    //setup outputs
    path "caveman.cfg.ini", emit: configFile
    path "alg_bean", emit: algBeanFile

  script:
    def applySW = includeSW != "NO_SW" ? "-w $includeSW" : ''
    def applySE = includeSE != "NO_SE" ? "-z $includeSE" : ''
    def applyDup = includeDups != "NO_DUPS" ? "-u $includeDups" : ''
    def applyNCN = normcn != "NO_NCN" ? "-j $normcn" : ''
    def applyTCN = tumcn != "NO_TCN" ? "-e $tumcn" : ''
    """ 
    caveman setup \
    ${applySW} \
    ${applySE} \
    ${applyDup} \
    ${applyNCN} \
    ${applyTCN} \
    -l splitList \
    -a alg_bean \
    -c caveman.cfg.ini \
    -t mt.bam \
    -n wt.bam \
    -r genome.fa.fai \
    -g ign.file \
    -f results \
    """
}

process split {

  input:
    path 'caveman.cfg.ini'
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    path ('ign.file')
    each index
    val maxSplitReadNum

  publishDir "$baseDir/splits/", mode: 'symlink', pattern: 'splitList.*', overwrite: true

  output:
    path "splitList.*", emit: splitFiles
    path "caveman.cfg.ini", emit: configFile
    path "readpos.*", optional: true, emit: rposFiles

  script:
    """
    caveman split \
    -f caveman.cfg.ini \
    -i $index \
    -e $maxSplitReadNum \
    """
}

process split_concat {
  input:
    path splitFileList

  publishDir "$baseDir", mode: 'copy', pattern: 'splitList', overwrite: true

  output:
    path 'splitList', emit: splitList

  script:
  """
  for splitFile in ${splitFileList}
  do
    cat \$splitFile >> splitList
  done
  """
}

process mstep {
  input:
    path 'caveman.cfg.ini'
    path 'alg_bean'
    path 'splitList'
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    path ('ign.file')
    each index
    val mstepSplitSize
    val mstepMinBaseQual
    path rposFiles

  publishDir "$baseDir/results/", mode: 'symlink', pattern: 'results/*/*.covs', overwrite: true

  output:
    path "results/*/*.covs", emit: mstepResults

  script:
    """
    caveman mstep \
    -i $index \
    -f caveman.cfg.ini \
    -a $mstepSplitSize \
    -m $mstepMinBaseQual
    """

}

// process merge{
//   input:
//     path 'caveman.cfg.ini'
//     path 'alg_bean'
//     path 'splitList'
//     path 'results'

//   output:
//     path 'covs_arr', emit: covArrFile
//     path 'probs_arr', emit: probsArrFile

//   publishDir "$baseDir/", mode: 'copy', pattern: 'covs_arr', overwrite: true
//   publishDir "$baseDir/", mode: 'copy', pattern: 'probs_arr', overwrite: true

//   script:
//   """
//   caveman merge \
//   -c caveman.cfg.ini \
//   -p probs_arr \
//   -f covs_arr
//   """
// }

/*
process estep {
// estep -i %d -k %f -g %s -o %s -v %s -w %s -f %s -l %s -r %s
// const my $CAVEMAN_ESTEP_MUT_PRIOR_EXT => q{ -c %s};
// const my $CAVEMAN_ESTEP_SNP_PRIOR_EXT => q{ -d %s};
// const my $CAVEMAN_ESTEP_NPLATFORM_EXT => q{ -P %s};
// const my $CAVEMAN_ESTEP_TPLATFORM_EXT => q{ -T %s};
// caveman estep -i jobindex [-f file] [-m int] [-k float] [-b float] [-p float] [-q float] [-x int] [-y int] [-c float] [-d float] [-a int]

// -i  --index [int]                                Job index (e.g. from $LSB_JOBINDEX)

// Optional
// -f  --config-file [file]                         Path to the config file produced by setup. [default:'caveman.cfg.ini']
// -m  --min-base-qual [int]                        Minimum base quality for inclusion of a read position [default:11]
// -c  --prior-mut-probability [float]              Prior somatic probability [default:0.000006]
// -d  --prior-snp-probability [float]              Prior germline mutant probability [default:0.000100]
// -k  --normal-contamination [float]               Normal contamination of tumour [default:0.100000]
// -b  --reference-bias [float]                     Reference bias [default:0.950000]
// -p  --mut-probability-cutoff [float]             Minimum probability call for a somatic mutant position to be output [default:0.800000]
// -q  --snp-probability-cutoff [float]             Minimum probability call for a germline mutant position to be output [default:0.950000]
// -x  --min-tum-coverage [int]                     Minimum tumour coverage for analysis of a position [default:1]
// -y  --min-norm-coverage [int]                    Minimum normal coverage for analysis of a position [default:1]
// -a  --split-size [int]                           Size of section to retrieve at a time from bam file. Allows memory footprint tuning [default:50000].
// -s  --debug                                      Adds an extra output to a debug file. Every base analysed has an output
// -g  --cov-file [file]                            File location of the covariate array. [default:'covs_arr']
// -o  --prob-file [file]                           File location of the prob array. [default:'probs_arr']
// -v  --species-assembly [string]                  Species assembly (eg 37/GRCh37), required if bam header SQ lines do not contain AS and SP information.
// -w  --species [string]                           Species name (eg Human), required if bam header SQ lines do not contain AS and SP information.
// -n  --normal-copy-number [int]                   Copy number to use when filling gaps in the normal copy number file [default:2].
// -t  --tumour-copy-number [int]                   Copy number to use when filling gaps in the tumour copy number file [default:2].
// -l  --normal-protocol [string]                   Normal protocol. Ideally this should match -r but not checked (WGS|WXS|RNA|RNA-Seq|AMPLICON|TARGETED) [default:WGS].
// -r  --tumour-protocol [string]                   Tumour protocol. Ideally this should match -l but not checked (WGS|WXS|RNA|RNA-Seq|AMPLICON|TARGETED) [default:WGS].
// -P  --normal-platform [string]                   Normal platform. Overrides the values retrieved from bam header.
// -T  --tumour-platform [string]                   Tumour platform. Overrides the values retrieved from bam header.
// -M  --max-copy-number [int]                      Maximum copy number permitted. If exceeded the copy number for the offending region will be set to this value. [default:10].
// -h	help                                         Display this usage information.

}

process merge_results {
// mergeCavemanResults -s %s -o %s -f %s
// Usage:
//     mergeCavemanResults [options] [file(s)...]

//     Required parameters: --output -o File to output result to. --splitlist
//     -s File containing list of split section --file-match -f ls style
//     pattern match of files to be merged. e.g. 
//     --help -h Brief help message.
}

process add_ids {
//   cgpAppendIdsToVcf.pl -i %s -o %s 
//   cgpAppendIdsToVcf.pl --help
// Usage

// Usage:
//     cgpAppendIdsToVcf.pl [-h] -i this.vcf -o this_with_ids.vcf

//       General Options:

//       --help      (-h)       Brief documentation

//             --file      (-i)       The file to append IDs to.

//             --outFile   (-o)       The output filename

//             Optional parameters:

//             --idstart   (-g)       Will set a sequential id generator to the given integer value. If not present will assign UUIDs to each variant.

//             --version   (-v)       Prints version information.

//       Examples:

//         cgpAppendIdsToVcf.pl -f this.vcf -o this_with_ids.vcf -g 1
}

process flag {
  cgpFlagCaVEMan.pl -i %s -o %s -s %s -m %s -n %s -b %s -g %s -umv %s -ref %s -t %s -sa %s
}
*/
