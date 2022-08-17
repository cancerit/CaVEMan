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

def getUsableContigIndices(genomeFai, ignoreContigs) {
  fai_lines = file(genomeFai).readLines()
  def contig_fai = []
  for (fai_line in fai_lines){
    contig_fai.push(fai_line.split("\\s+")[0])
  }
  contig_fai = contig_fai.reverse()
  def ignore_contigs_list  = file(ignoreContigs).readLines()
  def unique_fai = (contig_fai + ignore_contigs_list ) - contig_fai.intersect( ignore_contigs_list )
  def indices = [] 
  for(contig in unique_fai){
    indices.add((contig_fai.findIndexOf { "$it" == "$contig" })+1)
  }
  return indices
}

include { setup; split; split_concat; mstep } from './modules/caveman-core.nf'

//mstep; merge; estep; combine;
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
    ignoreRegions: ${params.ignoreRegions}
  """.stripIndent()

  main:
    genome = tuple file(params.genomefa), file("${params.genomefa}.fai")
    mt = tuple file(params.mutBam), file("${params.mutBam}.bai")
    wt = tuple file(params.controlBam), file("${params.controlBam}.bai")
    ign_file = file(params.ignoreRegions)
    results = file(params.resultsDir)
    splitF = file(params.splitList)
    ignoreContigFile = file(params.ignoreContigs)
    
    subCaVEMan(
      //setup
      genome,
      mt,
      wt,
      ign_file,
      ignoreContigFile,
      //Optional setup - or empty file and blanket CN
      params.normcn,
      params.tumcn,
      params.configF,
      params.algBeanF,
      results,
      splitF,
      params.includeSW,
      params.includeSE,
      params.includeDups,

      params.maxSplitReadNum,

      params.mstepSplitSize,
      params.mstepMinBaseQual
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
    ignoreContigFile
    //Optional setup - or empty file and blanket CN
    normcn
    tumcn
    configF
    algBeanF
    results
    splitF
    includeSW
    includeSE
    includeDups

    //Optional spliut 
    maxSplitReadNum

    //Optional mstep
    mstepSplitSize
    mstepMinBaseQual

  main:
    setup(
      genome,
      mt,
      wt,
      ign_file,
      normcn,
      tumcn,
      // configF,
      // algBeanF,
      // results,
      // splitF,
      includeSW,
      includeSE,
      includeDups
    )
    //Only want the indices of usable contigs - subtract contents of the ignoreContigsFile from the genome.fa.fa file
    Channel.fromList(getUsableContigIndices("${params.genomefa}.fai", ignoreContigFile)).set({contig_index})
    split(
      setup.out.configFile,
      genome,
      mt,
      wt,
      ign_file,
      contig_index,
      maxSplitReadNum,
    )
    //Concatenate split files
    split_concat(
      split.out.splitFiles.collect(),
    )
    Channel.of(1..file("$baseDir/splitList").countLines()).set({mstep_index})
    mstep(
      setup.out.configFile,
      setup.out.algBeanFile, 
      split_concat.out.splitList,
      genome,
      mt,
      wt,
      ign_file,
      mstep_index,
      mstepSplitSize,
      mstepMinBaseQual,
      split.out.rposFiles.collect()
    )
    // merge(
    //   setup.out.configFile,
    //   setup.out.algBeanFile, 
    //   split_concat.out.splitList,
    //   mstep.mstepResults.collect()
    // )
}
