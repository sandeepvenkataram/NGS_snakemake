#####
#
#   Sandeep Venkataram
#   Jan 2023
#   
#   Snakemake pipeline for NGS analysis
#   Pipeline starts with raw FASTQ files
#   Does QC analysis of FASTQ files
#   Maps reads to a reference genome
#   Calls and filters variants
#   Aggregates all variants into a single VCF
#
#   Example FASTQ data should be downloaded from NCBI SRA and stored in FASTQ_Input/
#   SRR15034727
#   SRR15034729
#
#   Pipeline is designed to be run in a user-specified conda environment containing all necessary tools:
#
#   Bowtie2
#   GATK4
#   FASTQC
#   Samtools
#   Picardtools
#   VCFTools
#   BGZip
#   Tabix
#
#####

import pandas as pd
from snakemake.utils import validate
from os import path
import glob

configfile: "config.yaml"


validate(config, "config.schema.yaml")


samples = pd.read_table(config["samples"], sep=",")
samples.set_index("sample", drop=False, inplace=True)
sample_list = samples["sample"].tolist()

validate(samples, "samples.schema.yaml")


ref_prefix = config["reference_genome_dir"] + config["reference_genome_base"]


def num_files(wildcards):
    return (len(get_fastq_files(wildcards)))

def get_fastq_files(wildcards):
    return (glob.glob(config["fastq_dir"]+wildcards.sample+"*"))

def get_trimmomatic_mode(wildcards):
    my_num_files = num_files(wildcards)
    if my_num_files == 1:
        return "SE"
    return "PE"

def get_trimmomatic_input(wildcards):
    return(" ".join(get_fastq_files(wildcards)))

def get_trimmomatic_output(wildcards):
    my_files = get_fastq_files(wildcards)
    if len(my_files)==1:
        return (config["trimmed_fastq_dir"] + wildcards.sample + ".fastq.gz")
    else:
        return(" ".join([
            config["trimmed_fastq_dir"] + wildcards.sample + ".fastq.gz",
            config["trimmed_fastq_dir"] + wildcards.sample + "_unpaired.fastq.gz",
            config["trimmed_fastq_dir"] + wildcards.sample + "_2.fastq.gz",
            config["trimmed_fastq_dir"] + wildcards.sample + "_unpaired_2.fastq.gz"
        ]))

def get_bowtie_input(wildcards):
    my_files = get_fastq_files(wildcards)
    if len(my_files)==1:
        return ("-U "+config["trimmed_fastq_dir"] + wildcards.sample + ".fastq.gz")
    else:
        return(" ".join([
            "-1 "+config["trimmed_fastq_dir"] + wildcards.sample + ".fastq.gz",
            "-2 " + config["trimmed_fastq_dir"] + wildcards.sample + "_2.fastq.gz"
        ]))


rule all:
    input:
        expand("FASTQC/{sample}_1_fastqc.html", sample=sample_list),  #fastqc output
        config["variant_dir"] + "all_clones_merged_forVEP.vcf"  #aggregate all VCF files
    default_target: True


rule index_reference:
    input:
        fasta = ref_prefix + ".fasta",
        bed_mask = config["reference_genome_dir"] + config["repetitive_elements_bed"]
    output:
        ref_prefix + ".dict",
        ref_prefix + ".fasta.fai",
        ref_prefix + ".1.bt2",
        config["reference_genome_dir"] + config["repetitive_elements_bed"] + ".idx"
    message:
        "Indexing Reference Genome"
    conda:
        config["conda_env"]
    params:
        ref_prefix=ref_prefix,
    log:
        "logs/reference_indexing",
    script:
        "scripts/reference_indexing.sh"


rule fastqc:
    input:
        get_fastq_files
    output:
        "FASTQC/{sample}_1_fastqc.html",
        "FASTQC/{sample}_2_fastqc.html"
    message:
        "Generating FASTQC for {wildcards.sample}"
    conda:
        config["conda_env"]
    log:
        "logs/{sample}_fastqc.log",
    shell:
        "fastqc -o FASTQC/ -f fastq {input} &>{log}"


rule trimmomatic:
    input:
        get_fastq_files
    output:
        config["trimmed_fastq_dir"] + "{sample}.fastq.gz"
    params:
        input_string = get_trimmomatic_input,
        output_string = get_trimmomatic_output,
        trimmomatic_mode = get_trimmomatic_mode,
        trimmomatic_clip = config["trimmomatic_clip"]
    message:
        "Trimming FASTQ files with Trimmomatic for {wildcards.sample}"
    conda:
        config["conda_env"]
    threads: 4
    log:
        "logs/{sample}_trimmomatic.log",
    shell:
        "trimmomatic {params.trimmomatic_mode} -threads {threads} {params.input_string} {params.output_string} {params.trimmomatic_clip} &> {log}"


rule bowtie2_map:
    input:
        config["trimmed_fastq_dir"] + "{sample}.fastq.gz",
        ref_prefix + ".1.bt2",
        ref_prefix + ".fasta.fai"
    output:
        sam_file = config["mapped_genome_dir"] + "{sample}_bowtie2.sam",
        rg_out = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.bam",
        md_bam = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam",
        md_stats = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.stats",
        bam_idex = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam.bai",
    message:
        "Mapping FASTQ files with Bowtie2 for {wildcards.sample}"
    params:
        ref_prefix = ref_prefix,
        bowtie2_input = get_bowtie_input
    conda:
        config["conda_env"]
    threads: 4
    log:
        "logs/{sample}",
    script:
        "scripts/mapping.sh"


rule call_and_filter_variants:
    input:
        bam = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam",
        bam_bai = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam.bai",
        picard_index = ref_prefix + ".dict",
        bed_mask = config["reference_genome_dir"] + config["repetitive_elements_bed"],
        bed_index = config["reference_genome_dir"] + config["repetitive_elements_bed"] + ".idx",
        ref_fasta = ref_prefix + "." + config["reference_genome_fasta_suffix"],
    output:
        bam_out = config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.realigned.bam",
        vcf_recode_gz = config["variant_dir"] + "{sample}.recode.vcf.gz",
        vcf_recode_gz_tbi = config["variant_dir"] + "{sample}.recode.vcf.gz.tbi",
    message:
        "GATK HaplotypeCaller and Filtration for {wildcards.sample}"
    conda:
        config["conda_env"]
    params:
        prefix = config["variant_dir"] + "{sample}",
        ref_fasta = ref_prefix + "." + config["reference_genome_fasta_suffix"],
        filter_expression = config["variant_filter_expression"],
    log:
        "logs/{sample}",
    script:
        "scripts/variant_calling.sh"


rule variant_merge:
    input:
        vcfs = expand(
            config["variant_dir"] + "{sample}.recode.vcf.gz", sample=sample_list
        ),
    output:
        vcf = config["variant_dir"] + "all_clones_merged.vcf",
        vep = config["variant_dir"] + "all_clones_merged_forVEP.vcf",
    message:
        "Merge all vcf"
    conda:
        config["conda_env"]
    log:
        "logs/variant",
    script:
        "scripts/variant_merge.sh"