#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process SamToFastq {
  tag "${bam.fileName}"
  label 'gatk'
  publishDir "${params.outdir}/01_MarkAdapters/"

  input:  // reads.bam
  path(bam)

  output: // reads_interleaved.fq
  tuple val("${bam.simpleName}"), path("${bam.simpleName}_newR1.fq"), path("${bam.simpleName}_newR2.fq")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SamToFastq \
    --INPUT $bam \
    --FASTQ ${bam.simpleName}_newR1.fq \
    --SECOND_END_FASTQ ${bam.simpleName}_newR2.fq \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INCLUDE_NON_PF_READS true \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_newR1.fq
  touch ${bam.simpleName}_newR2.fq
  """
}

// INTERLEAVE=true
// USE_JDK_DEFLATER=true USE_JDK_INFLATER=true

process bwamem2_index {
  tag "${genome_fasta.simpleName}"
  label 'bwamem'
  publishDir "${params.outdir}/02_MapReads"

  input:
  path(genome_fasta)

  output: // [genome.fasta, [genome_index files]]
  tuple path("$genome_fasta"), path("${genome_fasta}*")

  script:
  """
  #! /usr/bin/env bash
  $bwamem2_app index $genome_fasta
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${genome_fasta}.variousindexfiles
  """
}

process bwamem2_mem {
  tag "$readname"
  label 'bwamem'
  publishDir "${params.outdir}/02_MapReads"

  input:
  tuple path(genome_fasta), path(genome_index), val(readname), path(readpairs)

  output: // reads_mapped_2_genome.bam
  tuple val("$readname"), path("${readname}_mapped.bam")

  script:
  """
  #! /usr/bin/env bash
  PROC1=\$((`nproc` * 3/4))
  $bwamem2_app mem -t \${PROC1} ${genome_fasta} ${readpairs} |\
     $samtools_app view --threads 1 -bS - > ${readname}_mapped.bam
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${readname}_mapped.bam
  """
}

process MergeBamAlignment {
  tag "$i_readname"
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"

  input:  // [readgroup, unmapped reads, mapped reads]
  tuple val(i_readname), path(read_unmapped), path(read_mapped), path(genome_fasta), path(genome_dict)

  output: // merged bam and bai files
  tuple path("${i_readname}_merged.bam"), path("${i_readname}_merged.bai")

  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" MergeBamAlignment \
  --REFERENCE_SEQUENCE $genome_fasta \
  --UNMAPPED_BAM ${read_unmapped} \
  --ALIGNED_BAM ${read_mapped} \
  --OUTPUT ${i_readname}_merged.bam \
  --CREATE_INDEX true \
  --ADD_MATE_CIGAR true \
  --CLIP_ADAPTERS false \
  --CLIP_OVERLAPPING_READS true \
  --INCLUDE_SECONDARY_ALIGNMENTS true \
  --MAX_INSERTIONS_OR_DELETIONS -1 \
  --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
  --ATTRIBUTES_TO_RETAIN XS \
  --USE_JDK_DEFLATER true \
  --USE_JDK_INFLATER true
  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${i_readname}_merged.bam
  touch ${i_readname}_merged.bai
  """
}

process MarkDuplicates {
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"
  input:
  tuple path(bam), path(bai), path(genome_fasta), path(genome_dict)
  output:
  path("${bam.simpleName}_marked.bam")
  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" MarkDuplicates \
    --INPUT $bam \
    --OUTPUT ${bam.simpleName}_marked.bam \
    --METRICS_FILE ${bam.simpleName}_marked.metrics \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    --ASSUME_SORT_ORDER "queryname" \
    --CREATE_MD5_FILE true
  """
}

process SortAndFixTags {
  label 'gatk'
  publishDir "${params.outdir}/03_PrepGATK"
  input:
  tuple path(bam), path(genome_fasta), path(genome_dict)
  output:
  tuple path("${bam.simpleName}_sorted.bam"), path("${bam.simpleName}_sorted.bai")
  script:
  """
  #! /usr/bin/env bash
  $gatk_app --java-options "${java_options}" SortSam \
    --INPUT $bam \
    --OUTPUT /dev/stdout \
    --SORT_ORDER "coordinate" \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false \
  | $gatk_app --java-options "${java_options}" SetNmMdAndUqTags \
    --INPUT /dev/stdin \
    --OUTPUT ${bam.simpleName}_sorted.bam \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE $genome_fasta
  """
}

process gatk_HaplotypeCaller {
  tag "$window"
  label 'gatk'
  clusterOptions = "-N 1 -n 36 -t 01:00:00"
  publishDir "${params.outdir}/04_GATK", mode: 'copy'
  errorStrategy 'retry'
  maxRetries 3

  input:  // [window, reads files ..., genome files ...]
  tuple val(window), path(bam), path(bai), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  BAMFILES=`echo $bam | sed 's/ / -I /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" HaplotypeCaller \
    -R $genome_fasta \
    -I \$BAMFILES \
    -L $window \
    --output ${window.replace(':','_')}.vcf \
    ${gatk_HaplotypeCaller_params}

  """

  stub:
  """
  touch ${window.replace(':','_')}.vcf
  """
}

process gatk_HaplotypeCaller_invariant {
  tag "$window:${bam.simpleName}"
  label 'gatk'
  publishDir "${params.outdir}/04_GATK", mode: 'copy'

  input:  // [window, reads files ..., genome files ...]
  tuple val(window), path(bam), path(bai), path(genome_fasta), path(genome_dict), path(genome_fai)

  output: // identified SNPs as a vcf file
  path("*.vcf")

  script:
  """
  #! /usr/bin/env bash
  BAMFILES=`echo $bam | sed 's/ / -I /g' | tr '[' ' ' | tr ']' ' '`
  $gatk_app --java-options "${java_options}" HaplotypeCaller \
    -ERC BP_RESOLUTION \
    -R $genome_fasta \
    -I \$BAMFILES \
    -L $window \
    --output ${bam.simpleName}_${window.replace(':','_')}.vcf \
    ${gatk_HaplotypeCaller_params}

  """

  stub:
  """
  #! /usr/bin/env bash
  touch ${bam.simpleName}_${window.replace(':','_')}.vcf
  """
}
// --include-invariant -ERC GVCF

