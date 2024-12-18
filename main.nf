#! /usr/bin/env nextflow

nextflow.enable.dsl=2

include { DNA_VARIANT_CALLING } from './subworkflows/local/dna_variant_calling/main.nf'

include { RNA_VARIANT_CALLING } from './subworkflows/local/rna_variant_calling/main.nf'

include { LONGREAD_VARIANT_CALLING } from './subworkflows/local/long_read_variant_calling/main.nf'

include { SNP_EFF_ANNOTATION } from './subworkflows/local/variant_annotation/main.nf'

def helpMsg() {
  log.info """
   Usage:
   The typical command for running the pipeline are as follows:

   DNAseq:
     nextflow run main.nf --genome GENOME.fasta --reads "*_{R1,R2}.fastq.gz" --seq "dna" -profile singularity
     nextflow run main.nf --genome GENOME.fasta --reads_file READ_PATHS.txt --seq "dna" -profile singularity

   RNAseq:
     nextflow run main.nf --genome GENOME.fasta --gtf "genes.gtf" --reads "*_{R1,R2}.fastq.gz" --seq "rna" -profile singularity
  
   PacBio Long Reads:
     nextflow run main.nf --genome GENOME.fasta --long_reads "*.fastq.gz" --seq "longread" -profile singularity

   Mandatory arguments:
    --seq                   Specify input sequence type as 'dna', 'rna', or 'longread' [default:'${params.seq}'].
    --genome                Reference genome fasta file, against which reads will be mapped to find Variant sites

   Read input arguments:
    --reads                 Paired-end reads in fastq.gz format, will need to specify glob (e.g. "*_{R1,R2}.fastq.gz")
    --reads_file            Text file (tab delimited) with three columns [readname left_fastq.gz right_fastq.gz]. Will need full path for files.
    --long_reads            Long read file in fastq.gz format, will need to specify glob (e.g. "*.fastq.gz")

   Optional analysis arguments:
    --invariant             Output invariant sites [default:false]
    --gtf                   Gene Transfer Format file, only required for RNAseq input [default:false]
    --window                Window size passed to bedtools for parallel GATK Haplotype calls [default:${params.window}]

  Optional annotation of variants arguments:
    --annotate              Determine if annotation steps will be executed [default: ${params.annotate}]
    --gff                   Provide path to Gene Feature Format file as it's required for the annotation pipeline [default: ${params.gff}]
    --vcf                   Provide path to VCF file, to skip variant calling to speedup annotation pipeline [default: ${params.vcf}]

   Optional configuration arguments:
    -profile                Configuration profile to use. Can use multiple (comma separated)
                            Available: local, slurm, singularity, docker [default:local]
    --container_img         Container image used for singularity and docker [default:'${params.container_img}']
    
   GATK:
    --gatk_app              Link to gatk executable [default: '$gatk_app']
    --java_options          Java options for gatk [default:'${java_options}']
    --gatk_cluster_options  GATK cluster options [default:'${params.gatk_cluster_options}']
    --gatk_haplotype_caller_params  Additional parameters to pass to GATK HaplotypeCaller [default:'${params.gatk_haplotype_caller_params}']
    
   Aligners:
    --bwamem2_app           Link to bwamem2 executable [default: '$bwamem2_app']
    --star_app              Link to star executable [default: '$star_app']
    --star_index_params     Parameters for star index [default: '$star_index_params']
    --star_index_params     Parameters to pass to STAR index [default:'${params.star_index_params}']
    --star_index_file       Optional: speedup by providing a prebuilt STAR indexed genome [default: '${params.star_index_file}']
    --pbmm2_app             Link to pbmm2 executable [default: '$pbmm2_app']

   Other:
    --samtools_app          Link to samtools executable [default: '$samtools_app']
    --bedtools_app          Link to bedtools executable [default: '$bedtools_app']
    --datamash_app          Link to datamash executable [default: '$datamash_app']
    --vcftools_app          Link to vcftools executable [default: '$vcftools_app']

   Optional other arguments:
    --outdir                Output directory [default:'${params.outdir}']
    --threads               Threads per process [default:4 for local, 16 for slurm] 
    --queueSize             Maximum jobs to submit to slurm [default:${params.queueSize}]
    --account               HPC account name for slurm sbatch, atlas and ceres requires this
    --help                  Print this help message
 
"""
}

if(params.help){
  helpMsg()
  exit 0
}

def parameters_valid = ['help','outdir',
  'genome','gtf','reads','reads_file','long_reads','invariant','seq',
  'annotate','gff','vcf',
  'singularity_img','docker_img','container_img',
  'gatk_app','gatk_haplotype_caller_params',
  'star_app','star_index_params','star_index_file','bwamem2_app','samtools_app','bedtools_app','datamash_app','vcftools_app',
  'pbmm2_app','snpeff_params',
  'java_options','window','queueSize','queue-size','account', 'threads', 'gatk_cluster_options'] as Set

def parameter_diff = params.keySet() - parameters_valid
if (parameter_diff.size() != 0){
   exit 1, "[Pipeline error] Parameter(s) $parameter_diff is(are) not valid in the pipeline!\n"
}

workflow {
  // == Read in genome and reads channels
  if(params.genome) {
    genome_ch = channel.fromPath(params.genome, checkIfExists:true)
      | view {file -> "Genome file : $file "}
  } else {
    exit 1, "[Missing File(s) Error] This pipeline requires a reference '--genome [GENOME.fasta]' \n"
  }

  if (params.annotate) {
    if (params.gff && params.vcf) {
      gff_ch = channel.fromPath(params.gff, checkIfExists:true)
        | view {file -> "GFF file : $file "}
      vcf_ch = channel.fromPath(params.vcf, checkIfExists:true)
        | view {file -> "VCF file : $file "}
    } else {
      exit 1, "[Missing File(s) Error] Annotation mode requires both '--gff [annotation.gff]' and '--vcf [variants.vcf]' files\n"
    }
  } else if (params.seq != "longread") {
    reads_ch = channel.fromFilePairs(params.reads, checkIfExists:true)
      | view {files -> "Read files : $files "}
  } else if (params.reads_file) {
    reads_ch = channel.fromPath(params.reads_file, checkIfExists:true)
      | splitCsv(sep:'\t')
      | map { n -> [ n.getAt(0), [n.getAt(1), n.getAt(2)]] }
      | view {files -> "Read files : $files "}
  } else if (params.seq == "longread") {
    reads_ch = channel.fromPath(params.reads, checkIfExists:true)
      | view { files -> "Long read file : $files " }
  } else {
    exit 1, "[Missing File(s) Error] This pipeline requires either paired-end read files as a glob '--reads [*_{r1,r2}.fq.gz]' or as a tab-delimited text file '--reads_file [READS_FILE.txt]'\n"
  }

  if(!params.annotate){
    if (params.seq == "dna"){
      DNA_VARIANT_CALLING(genome_ch, reads_ch)
    } else if (params.seq == "rna") {
      if(params.gtf) {
        gtf_ch = channel.fromPath(params.gtf, checkIfExists:true)
      }else{
        exit 1, "[Missing File(s) Error] This pipeline requires a gtf file '--gtf [GENOME.gtf]' \n"
      }
      RNA_VARIANT_CALLING(genome_ch, reads_ch, gtf_ch)
    } else if (params.seq == "longread") {
      LONGREAD_VARIANT_CALLING(genome_ch, reads_ch)
    }
  } else {
    SNP_EFF_ANNOTATION(genome_ch, gff_ch, vcf_ch)
  }
}
