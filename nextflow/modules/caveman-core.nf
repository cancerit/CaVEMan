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
    val includeSW
    val includeSE
    val includeDups

  publishDir "$launchDir", mode: 'copy', pattern: 'caveman.cfg.ini', overwrite: true
  publishDir "$launchDir", mode: 'copy', pattern: 'alg_bean', overwrite: true

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

  publishDir "$launchDir/splits/", mode: 'symlink', pattern: 'splitList.*', overwrite: true

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

  publishDir "$launchDir/", mode: 'copy', pattern: 'splitList.tmp', overwrite: true

  output:
    path "splitList", emit: splitList

  script:
    """
    for splitFile in ${splitFileList}
    do
      cat \$splitFile >> splitList
    done
    """
}

process generate_mstep_indices {
  input:
    path 'fileForLineCount'

  output:
    stdout emit: mstep_index

  script:
    """
    sed -n '=' fileForLineCount
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
    val mstepSplitSize
    val minBaseQual
    path rposFiles
    each index

  publishDir "$launchDir/", mode: 'symlink', pattern: 'results/*/*.covs', overwrite: true

  output:
    path "results/*/*.covs", emit: mstepResults

  script:
    """
    caveman mstep \
    -f caveman.cfg.ini \
    -a $mstepSplitSize \
    -m $minBaseQual \
    -i $index
    """

}

process merge {
  input:
    path 'caveman.cfg.ini'
    path 'alg_bean'
    path 'splitList'
    path 'mstepResults'
    path 'results'

  output:
    path 'covs_arr', emit: covArrFile
    path 'probs_arr', emit: probsArrFile

  publishDir "$launchDir/", mode: 'copy', pattern: 'covs_arr', overwrite: true
  publishDir "$launchDir/", mode: 'copy', pattern: 'probs_arr', overwrite: true

  script:
    """
    caveman merge \
    -c covs_arr \
    -p probs_arr \
    -f caveman.cfg.ini
    """
}

process estep {

  input:
    path 'caveman.cfg.ini'
    path 'covs_arr'
    path 'probs_arr'
    path 'alg_bean'
    path 'splitList'
    path 'results'
    tuple path('genome.fa'), path('genome.fa.fai')
    tuple path('mt.bam'), path('mt.bam.bai')
    tuple path('wt.bam'), path('wt.bam.bai')
    path ('ign.file')
    path rposFiles
    val minBQ
    val spec
    val ass
    val priorMutProb
    val priorSnpProb
    val normalContamination
    val refBias
    val mutProbCutoff
    val snpProbCutoff
    val minTumCvg
    val minNormCvg
    val normProtocol
    val tumProtocol
    val normCNFill
    val tumCNFill
    val normSeqPlatform
    val tumSeqPlatform
    each index

  output:
    path "results/*/*.vcf.gz", emit: estepVCFs
    path "results/*/*.bed", emit: estepNoAnalysis
    val index, emit: estep_idx

  publishDir "$launchDir/", mode: 'symlink', pattern: 'results/*/*.vcf.gz', overwrite: true
  publishDir "$launchDir/", mode: 'symlink', pattern: 'results/*/*.bed', overwrite: true

  script:
    def applyNormPlat = normSeqPlatform != "NO_PLAT" ? "-P $normSeqPlatform" : ''
    def applyTumPlat = tumSeqPlatform != "NO_PLAT" ? "-T $tumSeqPlatform" : ''
    """
    caveman estep \
    -f caveman.cfg.ini \
    -g covs_arr \
    -o probs_arr \
    -m $minBQ \
    -v $ass \
    -w $spec \
    -c $priorMutProb \
    -d $priorSnpProb \
    -k $normalContamination \
    -b $refBias \
    -p $mutProbCutoff \
    -q $snpProbCutoff \
    -x $minTumCvg \
    -y $minNormCvg \
    -l $normProtocol \
    -r $tumProtocol \
    -n $normCNFill \
    -t $tumCNFill \
    ${applyNormPlat} \
    ${applyTumPlat} \
    -i $index
    """

}

process merge_results {
  input:
    path 'results'
    path 'splitList'
    val indexes
    val analysisName

  output:
    path "${analysisName}.muts.raw.vcf", emit: mergedMutsVCF
    path "${analysisName}.snps.raw.vcf", emit: mergedSnpVCF
    path "${analysisName}.no_analysis.bed.gz", emit: mergedNoAnalysis

  publishDir "$launchDir/", mode: 'copy', pattern: '*.no_analysis.bed.gz*', overwrite: true

  script:
    """
    mergeCavemanResults -s splitList -o ${analysisName}.muts.raw.vcf -f results/%/%.muts.vcf.gz && 
    mergeCavemanResults -s splitList -o ${analysisName}.snps.raw.vcf -f results/%/%.snps.vcf.gz && \
    mergeCavemanResults -s splitList -o ${analysisName}.no_analysis.bed -f results/%/%.no_analysis.bed && \
    bgzip ${analysisName}.no_analysis.bed && tabix -p bed ${analysisName}.no_analysis.bed.gz\
    """

}

process add_ids {
  input:
    path 'muts.raw.vcf'
    path 'snp.raw.vcf'
    val analysisName

  output:
    path "${analysisName}.muts.vcf.gz", emit: mutsVCFIIDs
    path "${analysisName}.muts.vcf.gz.tbi", emit: mutsVCFIIDsTbx
    path "${analysisName}.snps.vcf.gz", emit: snpsVCFIIDs
    path "${analysisName}.snps.vcf.gz.tbi", emit: snpsVCFIIDsTbx

  publishDir "$launchDir/", mode: 'copy', pattern: '*.snps.vcf.gz*', overwrite: true
  publishDir "$launchDir/", mode: 'copy', pattern: '*.muts.vcf.gz*', overwrite: true

  script:
    """
    cgpAppendIdsToVcf.pl -i muts.raw.vcf -o ${analysisName}.muts.vcf && \
    cgpAppendIdsToVcf.pl -i snp.raw.vcf -o ${analysisName}.snps.vcf && \
    bgzip ${analysisName}.muts.vcf && bgzip ${analysisName}.snps.vcf &&
    tabix -p vcf ${analysisName}.muts.vcf.gz && tabix -p vcf ${analysisName}.snps.vcf.gz
    """
}

/*
process flag {
  cgpFlagCaVEMan.pl -i %s -o %s -s %s -m %s -n %s -b %s -g %s -umv %s -ref %s -t %s -sa %s
}
*/
