"""
Author: Carlos Serna
Affiliation: Antimicrobial Resistance Unit (ARU) - UCM
Aim: Snakemake workflow to process short paired-end and long reads
Date: 01/07/2020
Latest modification:
"""

#########################################
# Preprocess and QC rules - Short reads #
#########################################
configfile: "config.yaml"

import io
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError


#---- VARIABLES ----#
# Input files (short reads)
INPUTDIR = config["project_dir"] + "/" + config["short_raw_data"]

# Output files
# FastQC
OUTDIR_FASTQC = config["project_dir"] + "/output_data/" + config["project"] + "/fastqc"
# Trimmomatic
OUTDIR_TRIMMOMATIC = config["project_dir"] + "/output_data/" + config["project"] + "/trimmomatic"
# MultiQC
OUTDIR_MULTIQC = config["project_dir"] + "/output_data/" + config["project"] + "/multiqc"

# Path trimmomatic - adapters
TRIMMOMATIC = config["project_dir"] + "/" + config["path_trimmomatic"]
ADAPTERS = config["project_dir"] + "/" + config["path_adapters"]

# Get wildcards from raw short data files (sample ID and forward/reverse reads)
SAMPLE_LIST,READS = glob_wildcards(INPUTDIR + "/{sample}_{read}.fastq.gz")
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
READS_SET = set(READS)


#---- RULES----#

rule all_preprocess:
    input:
        html = expand("{directory}/{sample}_{read}_fastqc.html", directory = OUTDIR_FASTQC, sample=SAMPLE_SET, read=READS_SET),
        zip = expand("{directory}/{sample}_{read}_fastqc.zip", directory = OUTDIR_FASTQC, sample=SAMPLE_SET, read=READS_SET),
        trimmedData = expand("{directory}/{sample}_{read}.trim.fastq.gz", directory = OUTDIR_TRIMMOMATIC, sample=SAMPLE_SET, read=READS_SET),
        html_trim = expand("{directory}/{sample}_{read}_trim_fastqc.html", directory = OUTDIR_FASTQC, sample=SAMPLE_SET, read=READS_SET),
        zip_trim = expand("{directory}/{sample}_{read}_trim_fastqc.zip", directory = OUTDIR_FASTQC, sample=SAMPLE_SET, read=READS_SET),
        notrim_html = OUTDIR_MULTIQC + "/report_raw_data_multiqc.html",
        trim_html = OUTDIR_MULTIQC + "/report_clean_data_multiqc.html"


rule fastqc_notrim:
    """
    Raw short reads quality control
    """
    input:
        INPUTDIR + "/{sample}_{read}.fastq.gz"
    output:
        html = OUTDIR_FASTQC + "/{sample}_{read}_fastqc.html",
        zip = OUTDIR_FASTQC + "/{sample}_{read}_fastqc.zip"
    params: ""
    threads: config["params_fastqc"]["threads"]
    wrapper:
        "0.63.0/bio/fastqc"


rule trimmomatic_pe:
    """
    Raw short reads trimming
    """
    input:
        r1 = INPUTDIR + "/{sample}_R1.fastq.gz",
        r2 = INPUTDIR + "/{sample}_R2.fastq.gz",
        trim_jar = os.path.join(TRIMMOMATIC, "trimmomatic-0.39.jar"),
        adapter_path = os.path.join(ADAPTERS, config["params"]["adapters"])
    output:
        # paired reads
        r1_paired = OUTDIR_TRIMMOMATIC + "/{sample}_R1.trim.fastq.gz",
        r2_paired = OUTDIR_TRIMMOMATIC + "/{sample}_R2.trim.fastq.gz",
        # unpaired redas
        r1_unpaired = OUTDIR_TRIMMOMATIC + "/{sample}_R1.unpaired.fastq.gz",
        r2_unpaired = OUTDIR_TRIMMOMATIC + "/{sample}_R2.unpaired.fastq.gz"
    log:
        trim_log = "logs/{sample}_trimmomatic.log"
    threads: config["params"]["threads"]
    resources: cpus=config["params"]["threads"]
    params:
    	leading = config["params"]["leading"],
    	trailing = config["params"]["trailing"],
    	window = config["params"]["window"],
    	minlen = config["params"]["minlen"]
    shell:
        """java -jar {input.trim_jar} PE \
        -phred33 -threads {threads} -trimlog {log.trim_log} {input.r1} {input.r2} \
        {output.r1_paired} \
        {output.r1_unpaired} \
        {output.r2_paired} \
        {output.r2_unpaired} \
        ILLUMINACLIP:{input.adapter_path}:2:30:10:1:TRUE\
        LEADING:{params.leading} TRAILING:{params.trailing} \
        SLIDINGWINDOW:{params.window} MINLEN:{params.minlen}"""


rule fastqc_trim:
    """
    Clean short reads quality control
    """
    input:
        OUTDIR_TRIMMOMATIC + "/{sample}_{read}.trim.fastq.gz"
    output:
        html = OUTDIR_FASTQC + "/{sample}_{read}_trim_fastqc.html",
        zip = OUTDIR_FASTQC + "/{sample}_{read}_trim_fastqc.zip"
    threads: config["params_fastqc"]["threads"]
    wrapper:
        "0.63.0/bio/fastqc"


rule multiqc:
    """
    Summary FastQC results (notrim and trim)
    """
    input:
        notrim = expand("{directory}/{sample}_{read}_fastqc.zip", directory= OUTDIR_FASTQC, sample=SAMPLE_SET, read=READS_SET),
        trim = expand("{directory}/{sample}_{read}_trim_fastqc.zip", directory= OUTDIR_FASTQC, sample=SAMPLE_SET, read=READS_SET)
    output:
        notrim_html = OUTDIR_MULTIQC + "/report_raw_data_multiqc.html",
        trim_html = OUTDIR_MULTIQC + "/report_clean_data_multiqc.html"
    conda:
        "../envs/qc_env.yaml"
    params:
        project = config["project"],
        project_dir = config["project_dir"],
        out_multiqc = OUTDIR_MULTIQC
    shell:
        """
        #notrim
        multiqc -n {output.notrim_html} {input.notrim} #run multiqc
        #repeat for trimmed data
        multiqc -n {output.trim_html} {input.trim} #run multiqc
        mkdir -p {params.project_dir}/results
        mkdir -p {params.project_dir}/results/{params.project}
        cp {output.notrim_html} {params.project_dir}/results/{params.project}
        cp {output.trim_html} {params.project_dir}/results/{params.project}
        cp {params.out_multiqc}/report_raw_data_multiqc_data/multiqc.log {params.project_dir}/logs/multiqc_raw.log
        cp {params.out_multiqc}/report_clean_data_multiqc_data/multiqc.log {params.project_dir}/logs/multiqc_clean.log
        """
