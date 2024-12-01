#! /usr/bin/env nextflow
nextflow.enable.dsl=2

process FastqToSam {
  tag "$readname"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // [readgroup, [left.fq.gz, right.fq.gz], increment_readgroup]
  tuple val(readname), path(readpairs), val(i_readname)

  output: // increment_readgroup.bam since one readgroup can have multiple lanes
  path("*.bam")

  script:
  """
  #! /usr/bin/env bash
  ${gatk_app} --java-options "${java_options}" FastqToSam \
    --FASTQ ${readpairs.getAt(0)} \
    --FASTQ2 ${readpairs.getAt(1)} \
    --OUTPUT ${i_readname}.bam \
    --READ_GROUP_NAME ${readname} \
    --SAMPLE_NAME ${readname}_name \
    --LIBRARY_NAME ${readname}_lib \
    --PLATFORM ILLUMINA \
    --SEQUENCING_CENTER ISU \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${i_readname}.bam
  """
}

process MarkIlluminaAdapters {
  tag "${bam.fileName}"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // reads.bam
  path(bam)

  output: // reads_marked.bam
  path "${bam.simpleName}_marked.bam"

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" MarkIlluminaAdapters \
    --INPUT $bam \
    --OUTPUT ${bam.simpleName}_marked.bam \
    --METRICS ${bam.simpleName}_marked_metrics.txt \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_marked.bam
  touch ${bam.simpleName}_marked_metrics.txt
  """
}

process CreateSequenceDictionary {
  tag "${genome_fasta.baseName}"
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"

  input:
  path(genome_fasta)

  output:
  path("${genome_fasta.baseName}.dict")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" CreateSequenceDictionary \
    -R ${genome_fasta} \
    -O ${genome_fasta.baseName}.dict
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta.baseName}.dict
  """
}

process MarkDuplicates {
  tag "$i_readname"
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"

  input:  // [readgroup, unmapped reads, mapped reads]
  tuple val(i_readname), path(merge_bam), path(merge_bai)

  output: // merged bam and bai files
  tuple val("${i_readname}"), path("${i_readname}_markduplicates.bam"), path("${i_readname}_markduplicates.bai")

  script:
  """
  ${gatk_app} \
       MarkDuplicates \
       --INPUT ${merge_bam} \
       --OUTPUT ${i_readname}_markduplicates.bam  \
       --CREATE_INDEX true \
       --VALIDATION_STRINGENCY SILENT \
       --METRICS_FILE ${i_readname}_markduplicates.metrics
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${i_readname}_markduplicates.bam
  touch ${i_readname}_markduplicates.bai
  """
}

process samtools_faidx {
  tag "${genome_fasta.simpleName}"
  label 'samtools'
  publishDir "${params.outdir}/03_PrepGATK"

  input:
  path(genome_fasta)

  output:
  path("${genome_fasta}.fai")

  script:
  """
  #! /usr/bin/env bash
  $samtools_app faidx $genome_fasta
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta}.fai
  """
}

process bedtools_makewindows {
  tag "${genome_fai.simpleName}"
  label 'bedtools'
  publishDir "${params.outdir}/03_PrepGATK"

  input:  // genome.fai
  path(genome_fai)

  output: // genome_coords.bed
  path("*_coords.bed")

  script:
  """
  #! /usr/bin/env bash
  cat ${genome_fai} | awk -F'\t' '{print \$1"\t"\$2}' > genome_length.txt
  $bedtools_app makewindows -w "$params.window" -g genome_length.txt |\
    awk '{print \$1"\t"\$2+1"\t"\$3}' |\
    sed \$'s/\t/:/1' |\
    sed \$'s/\t/-/g' > ${genome_fai.simpleName}_coords.bed
  """

  stub:
  """
  #! /usr/bin/env bash
  echo "0-1000" > ${genome_fai.simpleName}_coords.bed
  echo "1000-2000" >> ${genome_fai.simpleName}_coords.bed
  echo "2000-3000" >> ${genome_fai.simpleName}_coords.bed
  """
}
//  ${fai.simpleName}_coords.bed


// Consolidate GVCFs
process CombineGVCFs {
  tag "ConsolidateGVCFs"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple path(gvcf), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  GVCFFILES=`echo $gvcf | sed 's/ / --variant /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" CombineGVCFs \
    -R $genome_fasta \
    --variant \$GVCFFILES \
    --output all_combined.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch all_combined.vcf
  """
}

// Joint-Call Cohorts
process GenotypeGVCFs {
  tag "JointCallCohorts"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple path(all_combined_gvcf), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" GenotypeGVCFs \
    -R $genome_fasta \
    -V $all_combined_gvcf \
    --output output.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch output.vcf
  """
}

process merge_vcf {
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // multiple SNP vcf files
  path(vcfs)

  output: // merged into one vcf file
  path("first_round_merged.vcf")

  script:
  """
  #! /usr/bin/env bash
  set -euv
  
  cat ${vcfs.getAt(0)} | grep "^#" > first_round_merged.vcf
  cat ${vcfs} | grep -v "^#" >> first_round_merged.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch first_round_merged.vcf
  """
}

process vcftools_snp_only {
  tag "${merged_vcf.fileName}"
  label 'vcftools'
  publishDir "${params.outdir}/05_FilterSNPs", mode: 'copy'

  input:  // merged SNP vcf file
  path(merged_vcf)

  output: // vcf file only containing SNPs
  path("${merged_vcf.simpleName}_snps-only.*")

  script:
  """
  #! /usr/bin/env bash
  $vcftools_app \
    --vcf $merged_vcf \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --out ${merged_vcf.simpleName}_snps-only
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${merged_vcf.simpleName}_snps-only.outputfiles
  """
}

process SortVcf {
  tag "$vcf.fileName"
  label 'gatk'
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:  // [SNP.vcf, genome.dict]
  tuple path(vcf), path(dict)

  output: // sorted SNP.vcf
  path ("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SortVcf \
  --INPUT $vcf \
  --SEQUENCE_DICTIONARY $dict \
  --CREATE_INDEX true \
  --OUTPUT ${vcf.simpleName}_sorted.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${vcf.simpleName}_sorted.vcf
  """
}
// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process calc_DPvalue {
  tag "$sorted_vcf.fileName"
  label 'datamash'

  input:  // sorted SNP vcf
  path(sorted_vcf)

  output: // DP value (number) to filter vcf in downstream process
  stdout()

  script:
  """
  #! /usr/bin/env bash
  grep -v "^#" $sorted_vcf | cut -f 8 | grep -oe ";DP=.*" | cut -f 2 -d ';' | cut -f 2 -d "=" > dp.txt
  cat dp.txt | $datamash_app mean 1 sstdev 1 > dp.stats
  cat dp.stats | awk '{print \$1+5*\$2}'
  """

  stub:
  """
  #! /usr/bin/env bash
  # idk what default number should be here
  echo "75"
  """
}

process VariantFiltration {
  tag "$sorted_snp_vcf.fileName"
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:  // [sorted snp vcf, DP filter, genome files ... ]
  tuple path(sorted_snp_vcf), val(dp), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // filtered to identified SNP variants
  path("${sorted_snp_vcf.simpleName}.marked.vcf")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" VariantFiltration \
    --reference $genome_fasta \
    --sequence-dictionary $genome_dict \
    --variant $sorted_snp_vcf \
    --filter-expression \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > $dp\" \
    --filter-name "FAIL" \
    --output ${sorted_snp_vcf.simpleName}.marked.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${sorted_snp_vcf.simpleName}.marked.vcf
  """
}
// java -Xmx100g -Djava.io.tmpdir=$TMPDIR -jar

process keep_only_pass {
  tag "${snp_marked_vcf.fileName}"
  publishDir "$params.outdir/05_FilterSNPs", mode: 'copy'

  input:
  path(snp_marked_vcf)

  output:
  path("${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf")

  script:
  """
  #! /usr/bin/env bash
  cat $snp_marked_vcf | grep "^#" > ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  cat $snp_marked_vcf | grep -v "^#" | awk '\$7=="PASS"' >> ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${snp_marked_vcf.simpleName}_snp-only.pass-only.vcf
  """
}
