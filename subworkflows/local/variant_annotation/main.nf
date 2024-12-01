#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { snpEff_Build; snpEff_Annotate; } from '../../../modules/local/snpEff.nf'

workflow SNP_EFF_ANNOTATION {
  take:
  genome_ch
  gff_ch
  vcf_ch

  main:

  annotated_vcf_ch = genome_ch
  | combine(gff_ch)
  | snpEff_Build    // Building snpEff database
  | combine(genome_ch)
  | combine(vcf_ch)
  | snpEff_Annotate // Annotating variants using snpEff

  emit:
  annotated_vcf = annotated_vcf_ch
}
