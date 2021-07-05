rule bam2fragment_level_methylation_bedgz:
    input:
        bam_file = "bismark_mapped/{sample}.bam"
    params:
    output:
        read_level_methylation_status = "fragment_level_methylation/{sample}_fragment_methylation_vec.bed.gz"
    shell:
        "samtools view {input.bam_file} | python scripts/prepare_methylation_string.py | gzip - > {output.read_level_methylation_status}" 
rule overlapping_or_adjacent_reads:
    input:
        read_level_methylation_status = "fragment_level_methylation/{sample}_fragment_methylation_vec.bed.gz" 
    output:
        overlapping_or_adjacent = "overlapping_or_adjacent_fragments/{sample}_overlapping_or_adjacent.bed.gz"
    shell:
        "zcat {input.read_level_methylation_status}"
        " | egrep \"_overlapping|_adjacent\" "
        " | sort --parallel=8 -k1,1 -k2,2n -k3,3n |  gzip - > {output.overlapping_or_adjacent}"

