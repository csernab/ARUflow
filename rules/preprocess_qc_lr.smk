"""
Author: Carlos Serna
Affiliation: Antimicrobial Resistance Unit (ARU) - UCM
Aim: Snakemake workflow to process short paired-end and long reads
Date: 03/08/2020
Latest modification:
"""

########################################
# Preprocess and QC rules - Long reads #
########################################
configfile: "config.yaml"

import io
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError


#---- VARIABLES ----#
# Input files (short reads)
INPUTDIR_SR = config["project_dir"] + "/" + config["short_raw_data"]
INPUTDIR_LR = config["project_dir"] + "/" + config["long_raw_data"]

# Output files
# fitlong
OUTDIR_FILTLONG = config["project_dir"] + "/output_data/" + config["project"] + "/filtlong"
# nanoQC
OUTDIR_NANOQC = config["project_dir"] + "/output_data/" + config["project"] + "/nanoQC"

# Get wildcards from raw short data files (sample ID and forward/reverse reads)
SAMPLE_LIST,READS = glob_wildcards(INPUTDIR_SR + "/{sample}_{read}.fastq.gz")
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
READS_SET = set(READS)


#---- RULES----#

rule all_preprocess_lr:
    input:
        fastq_nanopore_filt = expand("{directory}/{sample}.filt.fastq.gz", directory = OUTDIR_FILTLONG, sample=SAMPLE_SET),
        html_nanoqc = expand("{directory}/{sample}_nanoqc.html", directory = OUTDIR_NANOQC, sample=SAMPLE_SET)



rule filtlong:
    """
    Filtering long reads by quality
    """
    input:
        r1_illumina = INPUTDIR_SR + "/{sample}_R1.fastq.gz",
        r2_illumina = INPUTDIR_SR + "/{sample}_R2.fastq.gz",
        fastq_nanopore = INPUTDIR_LR + "/{sample}.fastq.gz"
    output:
        fastq_nanopore_filt = OUTDIR_FILTLONG + "/{sample}.filt.fastq.gz"
    conda:
        "../envs/qc_nanopore_env.yaml"
    params:
    	minlenght = config["params_filtlong"]["min_length"],
    	keeppercent = config["params_filtlong"]["keep_percent"],
    	targetbases = config["params_filtlong"]["target_bases"]
    shell:
        """
        export LC_ALL=C
        filtlong -1 {input.r1_illumina} -2 {input.r2_illumina} \
        --min_length {params.minlenght} --keep_percent {params.keeppercent} --target_bases {params.targetbases} \
        {input.fastq_nanopore} | gzip > {output.fastq_nanopore_filt}
        """


rule nanoqc:
    """
    """
    input:
        fastq_nanopore_filt = OUTDIR_FILTLONG + "/{sample}.filt.fastq.gz"
    output:
        outdir_nanoqc = temp(directory(OUTDIR_NANOQC + "/{sample}_nanoqc")),
        html_nanoqc = OUTDIR_NANOQC + "/{sample}_nanoqc.html"
    conda:
        "../envs/qc_nanopore_env.yaml"
    shell:
        """
        nanoQC -o {output.outdir_nanoqc} {input.fastq_nanopore_filt}
        mv {output.outdir_nanoqc}/nanoQC.html {output.html_nanoqc}
        mv {output.outdir_nanoqc}/NanoQC.log logs/{wildcards.sample}.nanoQC.log
        """
