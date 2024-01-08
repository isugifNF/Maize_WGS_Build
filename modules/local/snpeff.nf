#! /usr/bin/env nextflow

nextflow.enable.dsl=2

process snpEff_Build {
  tag "${genome}"
  label 'snpeff'
  publishDir "${params.outdir}/snpEff_Build"
  
  input:
  tuple val(genome), path(genome_fasta), path(gff_file)

  output:
  tuple path("${genome}/snpEff.config"), path("${genome}/data/${genome}")

  script:
  """
  mkdir -p data/${genome}
  cp $gff_file data/${genome}/genes.gff
  cp $fasta_file data/${genome}/sequences.fa
  echo "${genome}.genome : ${genome}" >> snpEff.config
  snpEff build -gff3 -v ${genome}
  """
}

process snpEff_Annotate {
  tag "${vcf_file.simpleName}"
  label 'snpeff'
  publishDir "${params.outdir}/snpEff_Annotate"

  input:
  tuple val(genome), path(vcf_file)

  output:
  path("${vcf_file.simpleName}.eff.vcf")
  # other files
  
  script:
  """
  snpEff eff -v ${genome} -dataDir data/ $vcf_file > ${vcf_file.simpleName}.eff.vcf
  """
}
