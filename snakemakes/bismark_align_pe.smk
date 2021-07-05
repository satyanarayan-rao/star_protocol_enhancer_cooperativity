rule bismark_align_pe:
    input:
        read1 = "trimmed/{sample}_R1_val_1.fq.gz", 
        read2 = "trimmed/{sample}_R2_val_2.fq.gz"
    params:
        genome = "/beevol/home/satyanarr/workplace/data/ucsc/dm/dm3/", 
        out_dir = "bismark_mapped",
        extra = "--gzip --no_dovetail -p 4",
        basename = lambda wildcards: "--basename {b}".format(b = wildcards.sample)
    output:
        "bismark_mapped/{sample}_pe.bam", 
	"bismark_mapped/{sample}_PE_report.txt"
    script:
        "wrapper/bismark/align_pe.py"

rule bismark2report: 
    input: 
        "bismark_mapped/{sample}_PE_report.txt"
    output: 
        "bismark_align_report/{sample}.html"
    shell:
        "bismark2report --output {output} --alignment_report {input}" 

rule sort_bismark_bam:
    input:
        bismark_mapped_bam = "bismark_mapped/{sample}_pe.bam", 
    params:
    output:
        sorted_bam = "bismark_mapped/{sample}_pe_sorted.bam"
    shell:
        "samtools sort -n {input.bismark_mapped_bam} -o {output.sorted_bam}" 
rule merge_sorted_bam: 
    input:
        bam_file_list = lambda wildcards: config["bam_merge_config"][wildcards.cell_line]
    params:
    output:
        merged_bam = "bismark_mapped/merged_{cell_line}.bam"
    shell:
        "sh scripts/merge_bam.sh \"{input.bam_file_list}\" {output.merged_bam} {wildcards.cell_line}"
