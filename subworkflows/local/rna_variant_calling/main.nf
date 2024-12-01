#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { FastqToSam;
          MarkIlluminaAdapters;
          CreateSequenceDictionary;
          MarkDuplicates;
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

include { SamToFastq as SamToFastq_RNA;
          STAR_index;
          STAR_align;
          MergeBamAlignment as MergeBamAlignment_RNA; 
          SplitNCigarReads;
          gatk_HaplotypeCaller as gatk_HaplotypeCaller_RNA; } from '../../../modules/local/RNAseq.nf'

workflow RNA_VARIANT_CALLING {
  take:
  genome_ch
  reads_ch
  gtf_ch

  main:
  // == Since one sample may be run on multiple lanes
  i = 1

  // == Prepare mapped and unmapped read files
  cleanreads_ch = reads_ch
    | map { n -> [n.getAt(0), n.getAt(1), "${i++}_"+n.getAt(0)] }
    | FastqToSam
    | MarkIlluminaAdapters
    | SamToFastq_RNA
    | map { n -> [ n.getAt(0).replaceFirst("_marked",""), [ n.getAt(1), n.getAt(2)] ] }

  if(params.star_index_file){
    star_index_ch = channel.fromFile(params.star_index_file, checkIfExists:true)
  } else {
    star_index_ch = genome_ch
     | combine(gtf_ch)
     | STAR_index
  }

  mapped_ch = star_index_ch
    | combine(cleanreads_ch)
    | STAR_align
    | map { n -> [ n.getAt(0), n.getAt(1)]}

  unmapped_ch = MarkIlluminaAdapters.out
    | map { n -> [n.simpleName.replaceFirst("_marked",""), n] }

  genome_ch 
    | (CreateSequenceDictionary & samtools_faidx )

  MergeBamAlignment_ch = unmapped_ch 
    | join(mapped_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | MergeBamAlignment_RNA

  MergeBamAlignment_ch
    | MarkDuplicates
    | combine(genome_ch)
    | combine(samtools_faidx.out)
    | combine(CreateSequenceDictionary.out)
    | SplitNCigarRead
  allbai_ch = SplitNCigarReads.out
    | map { n -> n.getAt(2)}
    | collect
    | map { n -> [n]}
  allbambai_ch = SplitNCigarReads.out
    | map { n -> n.getAt(1)}
    | collect
    | map { n -> [n]}
    | combine(allbai_ch)

  // == Run Gatk Haplotype by interval window
  def split_outputs_ch = SplitNCigarReads.out.map { n -> [n[1], n[2]] }

  reads_part1_ch = samtools_faidx.out
    | bedtools_makewindows
    | splitText(){it.trim()}
    | combine(split_outputs_ch)
    | combine(genome_ch)
    | combine(CreateSequenceDictionary.out)
    | combine(samtools_faidx.out)
  
  reads_part1_ch
    | gatk_HaplotypeCaller_RNA
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

  emit:
  vcf = vcf_out_ch
}