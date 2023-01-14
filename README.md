 # NGS_snakemake

Sandeep Venkataram

Jan 2023
  
## Snakemake pipeline for NGS analysis

This pipeline starts with raw FASTQ files and generates a single VCF containing variants that pass quality filters.

Example FASTQ data should be downloaded from NCBI SRA and stored in FASTQ_Input/. I recommend using the fasterq-dump tool from the sra-toolkit available from NCBI.

    SRR15034727
    
    SRR15034729

This pipeline is designed to be run in a user-specified conda environment (specified in config.yaml "conda_env") containing all necessary tools:

    Bowtie2

    GATK4
    
    FASTQC
    
    Samtools
    
    Picardtools
    
    VCFTools
    
    BGZip
    
    Tabix
