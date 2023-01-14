#!/usr/bin/env bash


bowtie2 -x ${snakemake_params[ref_prefix]} ${snakemake_params[bowtie2_input]} --sensitive -S ${snakemake_output[sam_file]}  &> ${snakemake_log[0]}_bowtie2.log
gatk AddOrReplaceReadGroups -I ${snakemake_output[sam_file]} -O ${snakemake_output[rg_out]} -LB lib1 -PL illumina -PU unit1 -SM ${snakemake_wildcards[sample]} --SORT_ORDER coordinate  &> ${snakemake_log[0]}_add_RG.log
gatk MarkDuplicates -I ${snakemake_output[rg_out]} -O ${snakemake_output[md_bam]} -M ${snakemake_output[md_stats]}  &> ${snakemake_log[0]}_mark_dups.log
samtools index ${snakemake_output[md_bam]}  &> ${snakemake_log[0]}_samtools_index.log