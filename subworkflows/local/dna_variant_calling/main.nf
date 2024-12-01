#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { FastqToSam;
          MarkIlluminaAdapters;
          CreateSequenceDictionary;
          samtools_faidx;
          bedtools_makewindows;
          CombineGVCFs;
          GenotypeGVCFs;
          merge_vcf;
          vcftools_snp_only;
          SortVcf;
          calc_DPvalue;
          VariantFiltration;
          keep_only_pass; } from '../../../modules/local/GATK.nf'

include { SamToFastq as SamToFastq_DNA
          bwamem2_index;
          bwamem2_mem; 
          MergeBamAlignment as MergeBamAlignment_DNA;
          MarkDuplicates as MarkDuplicates_DNA;
          SortAndFixTags as SortAndFixTags_DNA;
          gatk_HaplotypeCaller as gatk_HaplotypeCaller_DNA;
          gatk_HaplotypeCaller_invariant; } from '../../../modules/local/DNAseq.nf'

workflow DNA_VARIANT_CALLING {
  take:
  genome_ch
  reads_ch


  main:
  // == Since one sample may be run on multiple lanes
  i = 1

  // == Prepare mapped and unmapped read files
  cleanreads_ch = reads_ch
    | map { n -> [n.getAt(0), n.getAt(1), "${i++}_"+n.getAt(0)] }
    | FastqToSam
    | MarkIlluminaAdapters
    | SamToFastq_DNA
    | map { n -> [ n.getAt(0).replaceFirst("_marked",""), [ n.getAt(1), n.getAt(2)] ] }

  mapped_ch = genome_ch
    | bwamem2_index
    | combine(cleanreads_ch)
    | bwamem2_mem

  unmapped_ch = MarkIlluminaAdapters.out
    | map { n -> [n.simpleName.replaceFirst("_marked",""), n] }

  genome_ch 
    | (CreateSequenceDictionary & samtools_faidx )

  MergeBamAlignment_ch = unmapped_ch 
    | join(mapped_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | MergeBamAlignment_DNA
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | MarkDuplicates_DNA
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | SortAndFixTags_DNA

  if(params.invariant) {
    allbambai_ch = MergeBamAlignment.out // do these need to be merged by read?
  } else {
    allbai_ch = MergeBamAlignment_ch
      | map { n -> n.getAt(1)}
      | collect 
      | map { n -> [n]}

      allbambai_ch = MergeBamAlignment_ch
      | map { n -> n.getAt(0)}
      | collect
      | map { n -> [n]}
      | combine(allbai_ch) 
  }

  // == Run Gatk Haplotype by interval window
  reads_part1_ch = samtools_faidx.out
    | bedtools_makewindows
    | splitText(){it.trim()}
    | combine(allbambai_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
  

  // If invariant, stop at output.vcf (all sites). If not, only keep SNPs with SNP filtering.
  if(params.invariant){
    vcf_out_ch = reads_part1_ch
      | gatk_HaplotypeCaller_invariant
      | collect
      | map { n -> [n] }
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | CombineGVCFs
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | GenotypeGVCFs

  } else {
    reads_part1_ch
      | gatk_HaplotypeCaller_DNA
      | collect
      | merge_vcf
      | vcftools_snp_only
      | combine(CreateSequenceDictionary.out)
      | SortVcf
      | calc_DPvalue

    // == Filter resulting SNPs
    vcf_out_ch = SortVcf.out
      | combine(calc_DPvalue.out.map{n-> n.replaceAll("\n","")})
      | combine(genome_ch)
      | combine(CreateSequenceDictionary.out)
      | combine(samtools_faidx.out)
      | VariantFiltration
      | keep_only_pass
  }

  emit:
  vcf = vcf_out_ch
}
