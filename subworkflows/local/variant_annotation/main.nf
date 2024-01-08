#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { snpEff_Build; snpEff_Annotate; } from '../../../modules/local/snpEff.nf'

workflow SNP_EFF_ANNOTATION {
  take:
  genome_ch
  gff_ch
  fasta_ch
  vcf_ch

  main:
  // Building snpEff database
  genome_ch
  | combine(gff_ch)
  | combine(fasta_ch)
  | snpEff_Build

  // Annotating variants using snpEff
  annotated_vcf_ch = vcf_ch
    | combine(genome_ch)
    | snpEff_Annotate

  emit:
  annotated_vcf = annotated_vcf_ch
}
