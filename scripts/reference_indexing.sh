#!/usr/bin/env bash

picard CreateSequenceDictionary -R ${snakemake_input[0]} -O ${snakemake_output} &> ${snakemake_log[0]}_picard_dict.log
samtools faidx ${snakemake_input[0]}  &> ${snakemake_log[0]}_samtools_faidx.log
bowtie2-build ${snakemake_input[0]} ${snakemake_params[ref_prefix]} &> ${snakemake_log[0]}_bowtie2_index.log