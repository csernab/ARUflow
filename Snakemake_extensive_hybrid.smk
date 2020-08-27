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
INPUTDIR_SR = config["short_raw_data"]
# Long reads
INPUTDIR_LR = config["long_raw_data"]

# Output quality control
# FastQC
OUTDIR_FASTQC = config["project_dir"] + "/output_data/" + config["project"] + "/fastqc"
# MultiQC
OUTDIR_MULTIQC = config["project_dir"] + "/output_data/" + config["project"] + "/multiqc"
# Input clean reads
OUTDIR_TRIMMOMATIC = config["project_dir"] + "/output_data/" + config["project"] + "/trimmomatic"
# fitlong
OUTDIR_FILTLONG = config["project_dir"] + "/output_data/" + config["project"] + "/filtlong"
# nanoQC
OUTDIR_NANOQC = config["project_dir"] + "/output_data/" + config["project"] + "/nanoQC"

# Output assembly
# Unicycler
OUTDIR_UNICYCLER = config["project_dir"] + "/output_data/" + config["project"] + "/unicycler"
# Quast
OUTDIR_QUAST = config["project_dir"] + "/output_data/" + config["project"] + "/quast"

# Abricate
OUTDIR_ABRICATE_RESFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/abricate/resfinder_hybrid"
OUTDIR_ABRICATE_PLASMIDFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/abricate/plasmidfinder_hybrid"

# Get wildcards from clean short data files (sample ID and forward/reverse reads)
SAMPLE_LIST,READS = glob_wildcards(INPUTDIR_SR + "/{sample}_{read}.fastq.gz")
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
#
# subworkflow preprocess_qc_lr:
#     workdir:
#         "./rules"
#     snakefile:
#         "./rules/preprocess_qc_lr.smk"
#     configfile:
#         "config.yaml"

include: "rules/preprocess_qc_sr.smk"
include: "rules/preprocess_qc_lr.smk"
include: "rules/chord_diagram_hybrid_include.smk"

#---- RULES----#

rule all_fast:
    input:
        trimmedData = expand("{directory}/{sample}_{read}.trim.fastq.gz", directory = OUTDIR_TRIMMOMATIC, sample=SAMPLE_SET, read=READS_SET),
        notrim_html = OUTDIR_MULTIQC + "/report_raw_data_multiqc.html",
        trim_html = OUTDIR_MULTIQC + "/report_clean_data_multiqc.html",
        fastq_nanopore_filt = expand("{directory}/{sample}.filt.fastq.gz", directory = OUTDIR_FILTLONG, sample=SAMPLE_SET),
        html_nanoqc = expand("{directory}/{sample}_nanoqc.html", directory = OUTDIR_NANOQC, sample=SAMPLE_SET),
        out_unicycler = expand("{directory}/{sample}.unicycler", directory = OUTDIR_UNICYCLER, sample=SAMPLE_SET),
        assembly_unicycler = expand("{directory}/{sample}.fasta", directory = OUTDIR_UNICYCLER, sample=SAMPLE_SET),
        quast_report = config["project_dir"] + "/results/" + config["project"] + "/quast_report_hybrid.html",
        report_mlst = config["project_dir"] + "/results/" + config["project"] + "/MLST_report_hybrid.tsv",
        out_abricate_resfinder = expand(OUTDIR_ABRICATE_RESFINDER + "/{sample}.tsv", sample=SAMPLE_SET),
        out_abricate_plasmidfinder = expand(OUTDIR_ABRICATE_PLASMIDFINDER + "/{sample}.tsv", sample=SAMPLE_SET),
        summary_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder_hybrid.tsv",
        summary_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder_hybrid.tsv",
        plot_resfinder = "results/" + config["project"] + "/extensive_hybrid_resfinder_analysis.png",
        plot_plasmidfinder = "results/" + config["project"] + "/extensive_hybrid_plasmidfinder_analysis.png",
        plot_count = "results/" + config["project"] + "/extensive_hybrid_count_analysis.png",
        chord_folder = config["project_dir"] + "/results/" + config["project"] + "/chord_hybrid/"


rule unicycler:
    input:
        r1_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R1.trim.fastq.gz",
        r2_clean = OUTDIR_TRIMMOMATIC + "/{sample}_R2.trim.fastq.gz",
        long_reads_filt = OUTDIR_FILTLONG + "/{sample}.filt.fastq.gz"
    output:
        out_unicycler = directory(OUTDIR_UNICYCLER + "/{sample}.unicycler"),
        assembly_unicycler = OUTDIR_UNICYCLER + "/{sample}.fasta"
    conda:
        "envs/assembly_hybrid_env.yaml"
    shell:
        """
        unicycler -1 {input.r1_clean} -2 {input.r2_clean} -l {input.long_reads_filt} -o {output.out_unicycler} --no_correct
        mv {output.out_unicycler}/assembly.fasta {output.assembly_unicycler}
        cp {output.out_unicycler}/unicycler.log logs/{wildcards.sample}.unicycler.log
        """


rule quast:
    """
    """
    input:
        assembly_unicycler = expand(OUTDIR_UNICYCLER + "/{sample}.fasta", sample=SAMPLE_SET)
    output:
        quast_folder = temp(directory(OUTDIR_QUAST)),
        quast_report = config["project_dir"] + "/results/" + config["project"] + "/quast_report_hybrid.html"
    conda:
        "envs/assembly_hybrid_env.yaml"
    shell:
        """
        quast.py -o {output.quast_folder} {input.assembly_unicycler}
        mv {output.quast_folder}/report.html {output.quast_report}
        """


rule mlst:
    """
    """
    input:
        assembly_unicycler = expand(OUTDIR_UNICYCLER + "/{sample}.fasta", sample=SAMPLE_SET)
    output:
        report_mlst = config["project_dir"] + "/results/" + config["project"] + "/MLST_report_hybrid.tsv"
    conda:
        "envs/mlst_env.yaml"
    shell:
        """
        mlst {input.assembly_unicycler} > {output.report_mlst}
        """


rule abricate_resfinder:
    """
    """
    input:
        assembly_unicycler = OUTDIR_UNICYCLER + "/{sample}.fasta"
    output:
        out_abricate_resfinder = OUTDIR_ABRICATE_RESFINDER + "/{sample}.tsv"
    conda:
        "envs/abricate_env.yaml"
    params:
    	minid = config["abricate_params_hybrid"]["minid"],
    	mincov = config["abricate_params_hybrid"]["mincov"]
    shell:
        """
        abricate --db resfinder --minid {params.minid} --mincov {params.mincov} {input.assembly_unicycler} > {output.out_abricate_resfinder}
        """


rule abricate_plasmidfinder:
    """
    """
    input:
        assembly_unicycler = OUTDIR_UNICYCLER + "/{sample}.fasta"
    output:
        out_abricate_plasmidfinder = OUTDIR_ABRICATE_PLASMIDFINDER + "/{sample}.tsv"
    conda:
        "envs/abricate_env.yaml"
    params:
    	minid = config["abricate_params_hybrid"]["minid"],
    	mincov = config["abricate_params_hybrid"]["mincov"]
    shell:
        """
        abricate --db plasmidfinder --minid {params.minid} --mincov {params.mincov} {input.assembly_unicycler} > {output.out_abricate_plasmidfinder}
        """


rule abricate_summary_resfinder:
    """
    """
    input:
        out_abricate = expand(OUTDIR_ABRICATE_RESFINDER + "/{sample}.tsv", sample=SAMPLE_SET)
    output:
        summary_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder_hybrid.tsv"
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
        summary_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder_hybrid.tsv"
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
        out_report_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder_hybrid.tsv",
        out_report_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder_hybrid.tsv"
    output:
        plot_resfinder = "results/" + config["project"] + "/extensive_hybrid_resfinder_analysis.png",
        plot_plasmidfinder = "results/" + config["project"] + "/extensive_hybrid_plasmidfinder_analysis.png",
        plot_count = "results/" + config["project"] + "/extensive_hybrid_count_analysis.png",
        Rplots = temp("Rplots.pdf")
    conda:
        "envs/r.yaml"
    script:
        "scripts/extensive_sr_analysis.R"
