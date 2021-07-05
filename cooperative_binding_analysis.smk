import os
import sys 
import pandas as pd
samples = pd.read_table(config["sample_metadata"])
samples.index = samples["sample"]
def get_raw_fasta (wildcards):
    read1 = samples.loc[wildcards.sample, "paired_read_1_path"]
    read2 = samples.loc[wildcards.sample, "paired_read_2_path"]
    return [read1, read2]
    

include: "snakemakes/trim_galore_pe.smk"
include: "snakemakes/bismark_align_pe.smk"
include: "snakemakes/suppress_context.smk"
include: "snakemakes/process_bam.smk"
include: "snakemakes/generate_metadata.smk"
include: "snakemakes/map_fragments_to_regions.smk"
include: "snakemakes/define_footprints_and_fix_gaps.smk"
include: "snakemakes/footprint_matix.smk" 
include: "snakemakes/assign_binding_states.smk"
include: "snakemakes/plot_binding_states.smk"

rule fastqc:
    input:
        "trimmed/{sample}.fq.gz"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    wrapper:
        "0.42.0/bio/fastqc"
