"""
Author: Carlos Serna
Affiliation: Antimicrobial Resistance Unit (ARU) - UCM
Aim: Snakemake workflow to process short paired-end and long reads
Date: 01/07/2020
Latest modification:
"""

#################################################################################
# Extensive - short reads: Antimicrobial Resistance and Plasmids Identification #
#################################################################################
configfile: "config.yaml"


import io
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#---- VARIABLES ----#
# Short reads
INPUTDIR = config["short_raw_data"]

# Output quality control
# FastQC
OUTDIR_FASTQC = config["project_dir"] + "/output_data/" + config["project"] + "/fastqc"
# MultiQC
OUTDIR_MULTIQC = config["project_dir"] + "/output_data/" + config["project"] + "/multiqc"

# Input clean reads
OUTDIR_TRIMMOMATIC = config["project_dir"] + "/output_data/" + config["project"] + "/trimmomatic"

# Output assembly
# Shovill
OUTDIR_SHOVILL = config["project_dir"] + "/output_data/" + config["project"] + "/shovill"
# Quast
OUTDIR_QUAST = config["project_dir"] + "/output_data/" + config["project"] + "/quast"

# Abricate
OUTDIR_ABRICATE_RESFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/abricate/resfinder"
OUTDIR_ABRICATE_PLASMIDFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/abricate/plasmidfinder"

# Get wildcards from clean short data files (sample ID and forward/reverse reads)
SAMPLE_LIST,READS = glob_wildcards(INPUTDIR + "/{sample}_{read}.fastq.gz")
# Unique the output variables from glob_wildcards
SAMPLE_SET = set(SAMPLE_LIST)
READS_SET = set(READS)

#---- Sub-WORKFLOWS----#

# subworkflow preprocess_qc:
#     workdir:
#         "./rules"
#     snakefile:
#         "./rules/preprocess_qc_sr.smk"
#     configfile:
#         "config.yaml"

include: "rules/preprocess_qc_sr.smk"
include: "rules/chord_diagram_sr_include.smk"

#---- RULES----#

rule all_extensive_sr:
    input:
        trimmedData = expand("{directory}/{sample}_{read}.trim.fastq.gz", directory = OUTDIR_TRIMMOMATIC, sample=SAMPLE_SET, read=READS_SET),
        notrim_html = OUTDIR_MULTIQC + "/report_raw_data_multiqc.html",
        trim_html = OUTDIR_MULTIQC + "/report_clean_data_multiqc.html",
        out_shovill = expand("{directory}/{sample}.shovill", directory = OUTDIR_SHOVILL, sample=SAMPLE_SET),
        assembly_shovill = expand(OUTDIR_SHOVILL + "/{sample}.fa", sample=SAMPLE_SET),
        quast_report = config["project_dir"] + "/results/" + config["project"] + "/quast_report.html",
        report_mlst = config["project_dir"] + "/results/" + config["project"] + "/MLST_report.tsv",
        out_abricate_resfinder = expand(OUTDIR_ABRICATE_RESFINDER + "/{sample}.tsv", sample=SAMPLE_SET),
        out_abricate_plasmidfinder = expand(OUTDIR_ABRICATE_PLASMIDFINDER + "/{sample}.tsv", sample=SAMPLE_SET),
        summary_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder.tsv",
        summary_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder.tsv",
        plot_resfinder = "results/" + config["project"] + "/extensive_sr_resfinder_analysis.png",
        plot_plasmidfinder = "results/" + config["project"] + "/extensive_sr_plasmidfinder_analysis.png",
        plot_count = "results/" + config["project"] + "/extensive_sr_count_analysis.png",
        chord_folder = config["project_dir"] + "/results/" + config["project"] + "/chord_sr/"


rule shovill:
    """
    Shovill assembler con los short reads sin TRIMMING
    """
    input:
        r1_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R1.trim.fastq.gz",
        r2_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R2.trim.fastq.gz"
        # r1 = INPUTDIR + "/{sample}_R1.fastq.gz",
        # r2 = INPUTDIR + "/{sample}_R2.fastq.gz"
    output:
        out_shovill = directory(OUTDIR_SHOVILL + "/{sample}.shovill"),
        assembly_shovill = OUTDIR_SHOVILL + "/{sample}.fa"
    conda:
        "envs/assembly_sr_env.yaml"
    shell:
        """
        shovill --outdir {output.out_shovill} --R1 {input.r1_clean} --R2 {input.r2_clean} --noreadcorr
        mv {output.out_shovill}/contigs.fa {output.assembly_shovill}
        cp {output.out_shovill}/shovill.log logs/{wildcards.sample}.shovill.log
        """

rule quast:
    """
    """
    input:
        assembly_shovill = expand(OUTDIR_SHOVILL + "/{sample}.fa", sample=SAMPLE_SET)
    output:
        quast_folder = temp(directory(OUTDIR_QUAST)),
        quast_report = config["project_dir"] + "/results/" + config["project"] + "/quast_report.html"
    conda:
        "envs/assembly_sr_env.yaml"
    shell:
        """
        quast.py -o {output.quast_folder} {input.assembly_shovill}
        mv {output.quast_folder}/report.html {output.quast_report}
        """

rule mlst:
    """
    """
    input:
        assembly_shovill = expand(OUTDIR_SHOVILL + "/{sample}.fa", sample=SAMPLE_SET)
    output:
        report_mlst = config["project_dir"] + "/results/" + config["project"] + "/MLST_report.tsv"
    conda:
        "envs/mlst_env.yaml"
    shell:
        """
        mlst {input.assembly_shovill} > {output.report_mlst}
        """

rule abricate_resfinder:
    """
    """
    input:
        assembly_shovill = OUTDIR_SHOVILL + "/{sample}.fa"
    output:
        out_abricate_resfinder = OUTDIR_ABRICATE_RESFINDER + "/{sample}.tsv"
    conda:
        "envs/abricate_env.yaml"
    params:
    	minid = config["abricate_params"]["minid"],
    	mincov = config["abricate_params"]["mincov"]
    shell:
        """
        abricate --db resfinder --minid {params.minid} --mincov {params.mincov} {input.assembly_shovill} > {output.out_abricate_resfinder}
        """

rule abricate_plasmidfinder:
    """
    """
    input:
        assembly_shovill = OUTDIR_SHOVILL + "/{sample}.fa"
    output:
        out_abricate_plasmidfinder = OUTDIR_ABRICATE_PLASMIDFINDER + "/{sample}.tsv"
    conda:
        "envs/abricate_env.yaml"
    params:
    	minid = config["abricate_params"]["minid"],
    	mincov = config["abricate_params"]["mincov"]
    shell:
        """
        abricate --db plasmidfinder --minid {params.minid} --mincov {params.mincov} {input.assembly_shovill} > {output.out_abricate_plasmidfinder}
        """


rule abricate_summary_resfinder:
    """
    """
    input:
        out_abricate = expand(OUTDIR_ABRICATE_RESFINDER + "/{sample}.tsv", sample=SAMPLE_SET)
    output:
        summary_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder.tsv"
    conda:
        "envs/abricate_env.yaml"
    shell:
        """
        abricate --summary {input.out_abricate} > {output.summary_resfinder}
        """

rule abricate_summary_plasmidfinder:
    """
    """
    input:
        out_abricate = expand(OUTDIR_ABRICATE_PLASMIDFINDER + "/{sample}.tsv", sample=SAMPLE_SET)
    output:
        summary_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder.tsv"
    conda:
        "envs/abricate_env.yaml"
    shell:
        """
        abricate --summary {input.out_abricate} > {output.summary_plasmidfinder}
        """


rule analysis_plot:
    """
    """
    input:
        out_report_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder.tsv",
        out_report_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder.tsv"
    output:
        plot_resfinder = "results/" + config["project"] + "/extensive_sr_resfinder_analysis.png",
        plot_plasmidfinder = "results/" + config["project"] + "/extensive_sr_plasmidfinder_analysis.png",
        plot_count = "results/" + config["project"] + "/extensive_sr_count_analysis.png",
        Rplots = temp("Rplots.pdf")
    conda:
        "envs/r.yaml"
    script:
        "scripts/extensive_sr_analysis.R"
