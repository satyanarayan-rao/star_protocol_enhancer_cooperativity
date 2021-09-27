rule get_peak_pairs:
    input:
        mnase_peaks = "input_bed/{bed}.bed", 
        mapped_dsmf = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bed.gz", 
        genome_size_file = "metadata/dm3.chrom.sizes"
    params:
        distance_th = 30,
        flank_from_peak = 15
    output:
        mnase_peak_pairs = "region_for_cobinding_analysis/{sample}_to_{bed}_lf_{lf}_rf_{rf}_pairs_for_cobinding.bed", 
        mapped_dsmf_reads_to_pairs = "dsmf_reads_to_peak_pairs/{sample}_to_{bed}_lf_{lf}_rf_{rf}_dsmf_pairs.bed.gz"
    shell: 
        "sh scripts/generate_mnase_peak_pairs.sh {input.mnase_peaks}"
        " {input.mapped_dsmf} {input.genome_size_file} {params.distance_th} {params.flank_from_peak}"
        " {output.mnase_peak_pairs} {output.mapped_dsmf_reads_to_pairs}"

rule extend_dsmf_reads_for_cobinding_analysis:
    input:
        mapped_dsmf_reads_to_pairs = "dsmf_reads_to_peak_pairs/{sample}_to_{bed}_lf_{lf}_rf_{rf}_dsmf_pairs.bed.gz",
        footprint_dict = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.pkl",
        regions_metadata = "regions_metadata/{bed}.pkl"
    output:
        extended_fragments = "extended_dsmf_reads_to_peak_pairs/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}.bed.gz", 
        verbose = "extended_dsmf_reads_to_peak_pairs/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_verbose.bed.gz", 
    shell:
        "sh scripts/extend_fragments_for_cobinding.sh {input.mapped_dsmf_reads_to_pairs}"
        " {input.footprint_dict} {input.regions_metadata} {wildcards.lf}"
        " {wildcards.rf} {wildcards.lextend} {wildcards.rextend}"
        " {output.extended_fragments} {output.verbose}" 
rule extend_bedpe_for_cobinding:
    input:
        fragments_covering_flank = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bedpe.gz"
    params:
    output:
        extended_fragments = "extended_dsmf_reads_bedpe/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}.bedpe.gz", 
        verbose = "extended_dsmf_reads_bedpe/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_verbose.bedpe.gz"
    shell:
        "sh scripts/extend_fragments_for_cobinding_bedpe.sh {input.fragments_covering_flank}"
        " {wildcards.lf} {wildcards.rf} {wildcards.lextend} {wildcards.rextend}"
        " {output.extended_fragments} {output.verbose}" 
