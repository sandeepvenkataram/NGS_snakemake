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


configfile: "config.yaml"


validate(config, "config.schema.yaml")


samples = pd.read_table(config["samples"], sep=",")
samples.set_index("sample", drop=False, inplace=True)
sample_list = samples["sample"].tolist()

validate(samples, "samples.schema.yaml")


ref_prefix = config["reference_genome_dir"] + config["reference_genome_base"]
ref_suffixes = ["dict", "fasta.fai", "fasta.amb", "1.bt2"]


trimmomatic_mode = "PE"
if not "fastq_suffix_2" in samples.columns:
    trimmomatic_mode = "SE"

fastq_file_dict = dict()
trimmed_file_dict = dict()
bowtie2_input_dict = dict()
for sample in samples["sample"].tolist():
    result = config["fastq_dir"] + sample + "." + samples.at[sample, "fastq_suffix_1"]
    result_trimmed = config["trimmed_fastq_dir"] + sample + ".fastq.gz"
    result_bowtie = "-1 " + result_trimmed
    if trimmomatic_mode == "PE":
        result += (
            " "
            + config["fastq_dir"]
            + sample
            + "."
            + samples.at[sample, "fastq_suffix_2"]
        )

        result_trimmed += (
            " " + config["trimmed_fastq_dir"] + sample + "_unpaired.fastq.gz"
        )
        result_trimmed += " " + config["trimmed_fastq_dir"] + sample + "_2.fastq.gz"
        result_trimmed += (
            " " + config["trimmed_fastq_dir"] + sample + "_unpaired_2.fastq.gz"
        )

        result_bowtie += " -2 " + config["trimmed_fastq_dir"] + sample + "_2.fastq.gz"
    else:
        result_bowtie = "-U " + result_trimmed
    fastq_file_dict[sample] = result
    trimmed_file_dict[sample] = result_trimmed
    bowtie2_input_dict[sample] = result_bowtie


rule all:
    input:
        expand("FASTQC/{sample}_fastqc.html", sample=sample_list),  #fastqc output
        config["variant_dir"] + "all_clones_merged_forVEP.vcf"  #aggregate all VCF files
    default_target: True


rule index_reference:
    input:
        ref_prefix + ".fasta",
    output:
        ref_prefix + ".dict",
        ref_prefix + ".fasta.fai",
        ref_prefix + ".1.bt2"
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
        lambda wildcards: fastq_file_dict[wildcards.sample],
    output:
        "FASTQC/{sample}_fastqc.html",
    message:
        "Generating FASTQC for {sample}"
    conda:
        config["conda_env"]
    log:
        "logs/{sample}_fastqc.log",
    shell:
        "fastqc -o FASTQC/ -f fastq {input} &>{log}"


rule trimmomatic:
    input:
        config["fastq_dir"] + "{sample}." + samples.at[sample, "fastq_suffix_1"],
    output:
        config["trimmed_fastq_dir"] + "{sample}.fastq.gz",
    params:
        input_string=lambda wildcards: fastq_file_dict[wildcards.sample],
        output_string=lambda wildcards: trimmed_file_dict[wildcards.sample],
        trimmomatic_mode=trimmomatic_mode,
        trimmomatic_clip=config["trimmomatic_clip"],
    message:
        "Trimming FASTQ files with Trimmomatic for {sample}"
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
    output:
        samfile=config["mapped_genome_dir"] + "{sample}_bowtie2.sam",
        rg_out=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.bam",
        mdbam=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam",
        mdstats=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.stats",
        bamidex=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam.bai",
    message:
        "Mapping FASTQ files with Bowtie2 for {sample}"
    params:
        ref_prefix=ref_prefix,
        bowtie2_input=bowtie2_input_dict[sample],
        rg_out=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.bam",
        mdbam=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam",
        mdstats=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.stats",
    conda:
        config["conda_env"]
    threads: 4
    log:
        "logs/{sample}",
    script:
        "scripts/mapping.sh"


rule call_and_filter_variants:
    input:
        bam=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam",
        bambai=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.bam.bai",
        picardindex=ref_prefix + ".dict",
        bed_mask=config["reference_genome_dir"] + config["repetitive_elements_bed"],
        bed_index=config["reference_genome_dir"] + config["repetitive_elements_bed"] + ".idx",
        ref_fasta=ref_prefix + "." + config["reference_genome_fasta_suffix"],
    output:
        bamout=config["mapped_genome_dir"] + "{sample}_bowtie2.sorted.dups.realigned.bam",
        vcfrecodegz=config["variant_dir"] + "{sample}.recode.vcf.gz",
        vcfrecodegztbi=config["variant_dir"] + "{sample}.recode.vcf.gz.tbi",
    message:
        "GATK HaplotypeCaller and Filtration for {sample}"
    conda:
        config["conda_env"]
    params:
        prefix=config["variant_dir"] + "{sample}",
        ref_fasta=ref_prefix + "." + config["reference_genome_fasta_suffix"],
        filter_expression=config["variant_filter_expression"],
    log:
        "logs/{sample}",
    script:
        "scripts/variant_calling.sh"


rule variant_merge:
    input:
        vcfs=expand(
            config["variant_dir"] + "{sample}.recode.vcf.gz", sample=sample_list
        ),
    output:
        vcf=config["variant_dir"] + "all_clones_merged.vcf",
        vep=config["variant_dir"] + "all_clones_merged_forVEP.vcf",
    message:
        "Merge all vcf"
    conda:
        config["conda_env"]
    log:
        "logs/variant",
    script:
        "scripts/variant_merge.sh"