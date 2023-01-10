#!/usr/bin/env bash

vcf-merge -c any ${snakemake_input[vcfs]} > ${snakemake_output[vcf]} 2> ${snakemake_log[0]}_vcf_merge.log
cut -f 1-8 ${snakemake_output[vcf]} > ${snakemake_output[vep]} 2> ${snakemake_log[0]}_create_vep.log