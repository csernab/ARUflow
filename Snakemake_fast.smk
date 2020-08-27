"""
Author: Carlos Serna
Affiliation: Antimicrobial Resistance Unit (ARU) - UCM
Aim: Snakemake workflow to process short paired-end and long reads data
Date: 01/07/2020
Latest modification: 12/08/2020
"""

#################################################################################
# Fast: Antimicrobial Resistance, Virulence factors and Plasmids Identification #
#################################################################################
configfile: "config.yaml"


import io
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError


#---- VARIABLES ----#
# Databases
INPUTDIR = config["short_raw_data"]
RESFINDER = config["resfinder_db"]
PLASMIDFINDER = config["plasmidfinder_db"]

# Output quality control
# FastQC
OUTDIR_FASTQC = config["project_dir"] + "/output_data/" + config["project"] + "/fastqc"
# MultiQC
OUTDIR_MULTIQC = config["project_dir"] + "/output_data/" + config["project"] + "/multiqc"

# Input clean reads
OUTDIR_TRIMMOMATIC = config["project_dir"] + "/output_data/" + config["project"] + "/trimmomatic"

# Output ariba
OUTDIR_RESFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/ariba/resfinder"
OUTDIR_PLASMIDFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/ariba/plasmidfinder"


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

#---- RULES----#

rule all_fast:
    input:
        trimmedData = expand("{directory}/{sample}_{read}.trim.fastq.gz", directory = OUTDIR_TRIMMOMATIC, sample=SAMPLE_SET, read=READS_SET),
        notrim_html = OUTDIR_MULTIQC + "/report_raw_data_multiqc.html",
        trim_html = OUTDIR_MULTIQC + "/report_clean_data_multiqc.html",
        reports_resfinder = expand(OUTDIR_RESFINDER + "/{sample}.tsv", sample=SAMPLE_SET),
        reports_plasmidfinder = expand(OUTDIR_PLASMIDFINDER + "/{sample}.tsv", sample=SAMPLE_SET),
        out_report_resfinder = "results/" + config["project"] + "/resfinder_final_report.csv",
        out_report_plasmidfinder = "results/" + config["project"] + "/plasmidfinder_final_report.csv",
        plot_resfinder = "results/" + config["project"] + "/fast_resfinder_analysis.png",
        plot_plasmidfinder = "results/" + config["project"] + "/fast_plasmidfinder_analysis.png",



rule ariba_resfinder:
    """

    """
    input:
        database = RESFINDER,
        r1_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R1.trim.fastq.gz",
        r2_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R2.trim.fastq.gz"
    output:
        out_resfinder = temp(directory(OUTDIR_RESFINDER + "/{sample}.out")),
        reports_resfinder = OUTDIR_RESFINDER + "/{sample}.tsv"
    conda:
        "envs/fast_env.yaml"
    shell:
        """
        ariba run {input.database} {input.r1_clean} {input.r2_clean} {output.out_resfinder}
        mv {output.out_resfinder}/report.tsv {output.reports_resfinder}
        """


rule ariba_plasmidfinder:
    """

    """
    input:
        database = PLASMIDFINDER,
        r1_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R1.trim.fastq.gz",
        r2_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R2.trim.fastq.gz"
    output:
        out_plasmidfinder = temp(directory(OUTDIR_PLASMIDFINDER + "/{sample}.out")),
        reports_plasmidfinder = OUTDIR_PLASMIDFINDER + "/{sample}.tsv"
    conda:
        "envs/fast_env.yaml"
    shell:
        """
        ariba run {input.database} {input.r1_clean} {input.r2_clean} {output.out_plasmidfinder}
        mv {output.out_plasmidfinder}/report.tsv {output.reports_plasmidfinder}
        """


rule ariba_summary_resfinder:
   """
   """
    input:
        reports_final = expand(OUTDIR_RESFINDER + "/{sample}.tsv", sample=SAMPLE_SET)
    output:
        out_report_resfinder = "results/" + config["project"] + "/resfinder_final_report.csv",
        out_phandango = temp("results/" + config["project"] + "/resfinder_final_report.csv.phandango.csv")
    conda:
        "envs/fast_env.yaml"
    shell:
        """
        ariba summary --no_tree --cluster_cols assembled,ref_seq {output.out_report_resfinder} {input.reports_final}
        mv {output.out_report_resfinder}.csv {output.out_report_resfinder}
        """


rule ariba_summary_plasmidfinder:
   """
   """
    input:
        reports_final = expand(OUTDIR_PLASMIDFINDER + "/{sample}.tsv", sample=SAMPLE_SET)
    output:
        out_report_plasmidfinder = "results/" + config["project"] + "/plasmidfinder_final_report.csv",
        out_phandango = temp("results/" + config["project"] + "/plasmidfinder_final_report.csv.phandango.csv")
    conda:
        "envs/fast_env.yaml"
    shell:
        """
        ariba summary --no_tree --cluster_cols assembled,ref_seq {output.out_report_plasmidfinder} {input.reports_final}
        mv {output.out_report_plasmidfinder}.csv {output.out_report_plasmidfinder}
        """


rule analysis_plot:
    input:
        out_report_resfinder = "results/" + config["project"] + "/resfinder_final_report.csv",
        out_report_plasmidfinder = "results/" + config["project"] + "/plasmidfinder_final_report.csv"
    output:
        plot_resfinder = "results/" + config["project"] + "/fast_resfinder_analysis.png",
        plot_plasmidfinder = "results/" + config["project"] + "/fast_plasmidfinder_analysis.png",
        Rplots = temp("Rplots.pdf")
    conda:
        "envs/r.yaml"
    script:
        "scripts/fast_analysis.R"
