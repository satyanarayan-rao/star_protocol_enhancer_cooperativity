rule map_fragments_to_regions:
    input:
        overlapping_or_adjacent = "footprints_on_fragments/{sample}_with_wobble_1_min_fp_10.bed.gz",
        regions_bed = "input_bed/{bed}.bed"
    params:
    output:
        mapped_to_regions = "fragments_mapped_to_regions/{sample}_to_{bed}_mapped.bed.gz"
    shell:
        "sh scripts/intersect_fragments_to_loci.sh {input.overlapping_or_adjacent}"
        " {input.regions_bed} {output.mapped_to_regions}" 
rule select_fragments_covering_flanks:
    input:
        mapped_to_regions = "fragments_mapped_to_regions/{sample}_to_{bed}_mapped.bed.gz"
    params:
    output:
        fragments_covering_flank = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bed.gz"
    shell:
        "zcat {input.mapped_to_regions} | python scripts/fragments_spanning_flanks.py"
        " {wildcards.lf} {wildcards.rf} | gzip - > {output.fragments_covering_flank}"

rule occluded_edges_on_fragments:
    input:
        fragments_covering_flank = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bed.gz"
    params:
    output:
        occluded_edges = "occluded_edges/{sample}_to_{bed}_lf_{lf}_rf_{rf}_occluded.tsv",
        occluded_pkl = "occluded_edges/{sample}_to_{bed}_lf_{lf}_rf_{rf}_occluded.pkl",
    shell:
        "zcat {input.fragments_covering_flank} | python scripts/identify_occluded_edges.py"
        " {output.occluded_edges} {output.occluded_pkl}"  

######### Mapping fragments to bedpe - useful for independent visualization #############

rule map_fragments_to_flanked_bedpe:
    input:
        overlapping_or_adjacent = "footprints_on_fragments/{sample}_with_wobble_1_min_fp_10.bed.gz",
        regions_bed = "input_bed/{bed}.bedpe"
    params:
    output:
        fragments_covering_flank = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bedpe.gz"
    shell:
        "sh scripts/intersect_fragments_to_loci_bedpe.sh {input.overlapping_or_adjacent}"
        " {input.regions_bed} {output.fragments_covering_flank} {wildcards.lf} {wildcards.rf}" 

