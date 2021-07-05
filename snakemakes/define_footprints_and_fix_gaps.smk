rule define_footprints:
    input:
        overlapping_or_adjacent = "overlapping_or_adjacent_fragments/{sample}_overlapping_or_adjacent.bed.gz"
    params:
         wobble_gap = lambda wildcards: config["process_footprints"][wildcards.setting]["wobble_gap"],
         min_fp_len = lambda wildcards: config["process_footprints"][wildcards.setting]["min_fp_len"]
    output:
        footprints_on_fragments = "footprints_on_fragments/{sample}_with_{setting}.bed.gz",
    shell:
        "zcat {input.overlapping_or_adjacent} | python scripts/call_footprints.py"
        " | python scripts/fix_wobble.py {params.wobble_gap}"
        " | python scripts/keep_all_footprint_ge_th.py {params.min_fp_len}" 
        " | gzip - > {output.footprints_on_fragments}" 
