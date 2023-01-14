#!/usr/bin/env bash

gatk HaplotypeCaller -I ${snakemake_input[bam]} -O ${snakemake_params[prefix]}.vcf -R ${snakemake_params[ref_fasta]} -bamout ${snakemake_output[bam_out]}  &> ${snakemake_log[0]}_haplotype_caller.log
gatk VariantFiltration --variant ${snakemake_params[prefix]}.vcf --output ${snakemake_params[prefix]}.filtered.vcf -R ${snakemake_input[ref_fasta]} --filter-name low_quality --filter-expression "${snakemake_params[filter_expression]}" --mask-name repetitive_region --mask ${snakemake_input[bed_mask]}  &> ${snakemake_log[0]}_var_filter.log
vcftools --vcf ${snakemake_params[prefix]}.filtered.vcf --out ${snakemake_params[prefix]} --remove-filtered-all --recode-INFO-all --recode  &> ${snakemake_log[0]}_vcf_recode.log
bgzip -f -i ${snakemake_params[prefix]}.recode.vcf  &> ${snakemake_log[0]}_vcf_bgzip.log
tabix -f -p vcf ${snakemake_params[prefix]}.recode.vcf.gz  &> ${snakemake_log[0]}_vcf_tabix.log