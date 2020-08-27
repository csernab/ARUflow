"""
Author: Carlos Serna
Affiliation: Antimicrobial Resistance Unit (ARU) - UCM
Aim: Snakemake workflow to process short paired-end and long reads
Date: 01/07/2020
Latest modification:
"""

configfile: "config.yaml"


import io
import os
import pandas as pd
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

# Abricate
OUTDIR_ABRICATE_RESFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/abricate/resfinder"
OUTDIR_ABRICATE_PLASMIDFINDER = config["project_dir"] + "/output_data/" + config["project"] + "/abricate/plasmidfinder"

rule all_dir:
    input:
        chord_folder = config["project_dir"] + "/results/" + config["project"] + "/chord_sr/"

rule plasmids_chord:
    """
    """
    input:
        summary_resfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_resfinder.tsv",
        summary_plasmidfinder = config["project_dir"] + "/results/" + config["project"] + "/summary_abricate_plasmidfinder.tsv"
    #     outdir_input = config["project_dir"] + "/output_data/" + config["project"]
    output:
        chord_folder = directory(config["project_dir"] + "/results/" + config["project"] + "/chord_sr/")
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plasmids_chord.R"
