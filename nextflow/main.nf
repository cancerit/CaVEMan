/*
 * Copyright (c) 2014-2022 Genome Research Ltd
 *
 * Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
 *
 * This file is part of CaVEMan.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * 1. The usage of a range of years within a copyright statement contained within
 * this distribution should be interpreted as being equivalent to a list of years
 * including the first and last year specified and all consecutive years between
 * them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
 * 2009, 2011-2012’ should be interpreted as being identical to a statement that
 * reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
 * statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
 * identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
 * 2009, 2010, 2011, 2012’.
 */

nextflow.enable.dsl=2

def helpMessage() {
    log.info """
    Usage:
      nextflow run https://github.com/cancerit/CaVEMan -entry <entry-point> --help
    Available entry points:
      - caveman
      - caveman_np
    """.stripIndent()
}

def helpCavemanMessage() {
    log.info """
    Usage:
      nextflow run https://github.com/cancerit/CaVEMan -entry caveman --help
  Required parameters:
    --outdir        Folder to output result to.
    --genomefa      Path to reference index genome file *.fa.fai[.gz]
    --mutBam        Tumour BAM/CRAM file (co-located index and bas files)
    --controlBam    Normal BAM/CRAM file (co-located index and bas files)
    --ignoreContigs Location of tsv ignore regions file



    --simrep      Full path to tabix indexed simple/satellite repeats.
    --filter      VCF filter rules file (see FlagVcf.pl for details)
    --genes       Full path to tabix indexed coding gene footprints.
    --unmatched   Full path to tabix indexed gff3 of unmatched normal panel
                   - see pindel_np_from_vcf.pl
  Optional
    --seqtype    Sequencing protocol, expect all input to match [WGS]
    --assembly   Name of assembly in use
                  -  when not available in BAM header SQ line.
    --species    Species
                  -  when not available in BAM header SQ line.
    --exclude    Exclude this list of ref sequences from processing, wildcard '%'
                  - comma separated, e.g. NC_007605,hs37d5,GL%
    --badloci    Tabix indexed BED file of locations to not accept as anchors
                  - e.g. hi-seq depth from UCSC
    --skipgerm   Don't output events with more evidence in normal BAM.
    --softfil    VCF filter rules to be indicated in INFO field as soft flags
    --limit      When defined with '-cpus' internally thread concurrent processes.
                  - requires '-p', specifically for pindel/pin2vcf steps
    --debug      Don't cleanup workarea on completion.
    --apid       Analysis process ID (numeric) - for cgpAnalysisProc header info
                  - not necessary for external use
    """.stripIndent()
}

def get_chr_count (genome, ignorecontigs) {
  echo "$genome\t$ignorecontigs"
}


process setup {

  input:
    //setup
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    path ('ign.file')

    //Optional setup - or empty file and blanket CN
    val normcn
    val tumcn
    //tumcn
    val configF
    val algBeanF
    path ('results')
    val splitList
    val includeSW
    val includeSE
    val includeDups

  output:
    //setup outputs
    path "${configF}", emit: configFile
    path "${algBeanF}", emit: algBeanFile

  log.info """\
  CaVEMan:setup
  """.stripIndent()

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
    -l $splitList \
    -a $algBeanF \
    -c $configF \
    -t mt.bam \
    -n wt.bam \
    -r genome.fa.fai \
    -g ign.file \
    -f results \
    """
}

process split {

  input:
    path configFile
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    path ('ign.file')
    each index
    val maxSplitReadNum
    val splitList

  output:
    path "${splitList}*": emit splitFiles

  log.info """\
  CaVEMan:split
  """.stripIndent()

  script:
    """
    caveman split \
    -f $configFile \
    -i $index \
    -e $maxSplitReadNum \
    """
}

process split_concat {
  input:
    path 'splitList.*'
    path splitList

  output:
    path $splitList, emit: mergedSplitList

  log.info """\
  CaVEMan:split_concat
  """.stripIndent()

  script:
  """
  cat splitList.* > $splitList
  """
}

process mstep {
  input:
    path configFile
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    path ('ign.file')
    each index
    val mstepSplitSize
    val mstepMinBaseQual
    val splitList
    path ('results')

  output:

  log.info """\
  CaVEMan:mstep
  """.stripIndent()

  script:
    """
    caveman mstep \
    -i $index \
    -f $configFile \
    -a $mstepSplitSize \
    -m $mstepMinBaseQual \
    """

}
//    -i jobindex [-f path] [-c int] [-m int] [-e int] 
//   split -i %d -f %s -e %d
//   caveman split -i jobindex [-f path] [-m int] [-e int] 

// -i  --index [int]                 Job index (e.g. from $LSB_JOBINDEX)

// Optional
// -f  --config-file [file]          Path to the config file produced by setup [default:'caveman.cfg.ini'].
// -e  --read-count [int]            Guide for maximum read count in a section [default:350000]
// -h    help                          Display this usage information.
/*}

process split_concat {
	// # uncoverable subroutine
	// my $options = shift;
	// my $tmp = $options->{'tmp'};
	// my $out = $options->{'out_file'};
	// my $target = $options->{'target_files'};
	// return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	// my $command = sprintf('cat %s > %s',$target,$out);
	// PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	// return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

process mstep {
// mstep -i %d -f %s
// caveman mstep -i jobindex [-f file] [-m int] [-a int]

// -i  --index [int]                 Job index.

// Optional
// -f  --config-file [file]          Path to the config file produced by setup [default: 'caveman.cfg.ini'].
// -a  --split_size [int]            Size of section to retrieve at a time from bam file. Allows memory footprint tuning [default:50000].
// -m  --min-base-qual [int]         Minimum base quality for inclusion of a read position in analysis [default:11]
// -h	help                          Display this usage information.

}

process merge {
// merge -c %s -p %s -f %s
// caveman merge [-f file] [-c file] [-p file] 

// Optional
// -f  --config-file file               Path to the config file produced by setup. [default: 'caveman.cfg.ini']

// -c  --covariate-file filename        Location to write merged covariate array [default: covs_arr]
// -p  --probabilities-file filename    Location to write probability array [default: probs_arr]
// -h	help                             Display this usage information.
}

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

// Print help if no workflow specified
workflow { 
    if ( params.help ) {
        helpMessage()
        exit 0
    }
}

workflow caveman {
  //wf for running caveman
  //Help message
  if ( params.help ) {
    helpCavemanMessage()
    exit 0
  }
  log.info """\
    CaVEMan:caveman - NF Pipeline
    ------------------------------
    //setup
    genomefa: ${params.genomefa}
    mutBam: ${params.mutBam}
    controlBam: ${params.controlBam}
    ignoreContigs: ${params.ignoreContigs}
  """.stripIndent()

  main:
    genome = tuple file(params.genomefa), file("${params.genomefa}.fai")
    mt = tuple file(params.mutBam), file("${params.mutBam}.bai")
    wt = tuple file(params.controlBam), file("${params.controlBam}.bai")
    ign_file = file(params.ignoreContigs)
    results = file(params.resultsDir)
    
    subCaVEMan(
      //setup
      genome,
      mt,
      wt,
      ign_file,

      //Optional setup - or empty file and blanket CN
      params.normcn,
      params.tumcn,
      params.configF,
      params.algBeanF,
      results,
      params.splitList,
      params.includeSW,
      params.includeSE,
      params.includeDups,

      params.maxSplitReadNum
    )
}

workflow subCaVEMan {
  //sub workflow, params and help in entrypoint wf
  take:
    //setup
    genome
    mt
    wt
    ign_file

    //Optional setup - or empty file and blanket CN
    normcn
    tumcn
    configF
    algBeanF
    results
    splitList
    includeSW
    includeSE
    includeDups

    //Optional spliut 
    maxSplitReadNum

  main:
    setup(
      genome,
      mt,
      wt,
      ign_file,
      normcn,
      tumcn,
      configF,
      algBeanF,
      results,
      splitList,
      includeSW,
      includeSE,
      includeDups
    )
    // contigs(
    //   genome
    // )
    Channel.of(1..file("${params.genomefa}.fai").countLines()).set({contig_index})
    split(
      setup.out.configFile,
      genome,
      mt,
      wt,
      ign_file,
      contig_index,
      maxSplitReadNum,
      splitList
    ) | collect | 
    split_concat(

    )
    mstep(
      setup.out.configFile,
      genome,
      mt,
      wt,
      ign_file,
      split.out.splitFiles,
      mstepSplitSize,
      mstepMinBaseQual,
      splitList,
      results
    )
}
