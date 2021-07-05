rule suppress_context:
    input:
        bam_file = "bismark_mapped/merged_{cell_line}.bam"
    output:
        out_bam = "bismark_mapped/suppressed_merged_{cell_line}.bam"
    shell:
        "sh scripts/suppress_context.sh {input.bam_file} {output.out_bam}"
