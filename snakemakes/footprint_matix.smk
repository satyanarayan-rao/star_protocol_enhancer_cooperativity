rule build_fooptrint_dict: 
    input:
        fragments_covering_flank = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bed.gz"
    params:
    output:
        footprint_dict = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.pkl"
    shell:
        "zcat {input.fragments_covering_flank}"
        " | python scripts/build_footprint_dict.py {output.footprint_dict}"
 

rule build_footprint_matrix:
    input:
        fragments_covering_flank = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.bed.gz"
    params:
    output:
        methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix.tsv" ,
        methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix_flen.pkl",
        methylation_matrix_footprint_vec_pkl = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix_footprint.pkl",
        strand_agnostic_footprint_vec_pkl = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix_footprint.strand.agnostic.pkl",
        strand_agnostic_footprint_matrix = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix.strand.agnostic.tsv",
        footprint_length_at_bp_resolution_pkl = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix.footprint_len.bp.res.pkl",
        footprint_length_at_bp_resolution_tsv = \
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix.footprint_len.bp.res.tsv"        
    shell: 
        "zcat {input.fragments_covering_flank} | python scripts/build_flank_methylation_matrix.py"
        " {wildcards.lf} {wildcards.rf} {output.methylation_matrix}"
        " {output.methylation_matrix_real_footprint_len_pkl}"
        " {output.methylation_matrix_footprint_vec_pkl}"
        " {output.strand_agnostic_footprint_vec_pkl}"
        " {output.strand_agnostic_footprint_matrix}"
        " {output.footprint_length_at_bp_resolution_pkl}"
        " {output.footprint_length_at_bp_resolution_tsv}" 
rule get_footprint_and_percentage_methylation:
    input:
        methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix.tsv",
        methylation_matrix_real_footprint_len_pkl =\
              "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix_flen.pkl",
        footprint_dict = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.pkl",
        bed_file_for_strand = "input_bed/{bed}.bed"
    params:
        neighbor_left_from_peak_center = 30,
        neighbor_right_from_peak_center = 30, 
        
    output:
        footprint_length_and_percentage = \
          "methylation_levels_around_roi/{sample}_to_{bed}_lf_{lf}_rf_{rf}_footprint_fraction_and_length.tsv",
        methylation_level = \
          "methylation_levels_around_roi/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_level.tsv",
        methylation_level_dict = \
          "methylation_levels_around_roi/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_level.pkl",
    shell:
        "python scripts/get_methylation_level_and_footprint_length.py"
        " {input.methylation_matrix} {input.methylation_matrix_real_footprint_len_pkl}"
        " {input.footprint_dict} {input.bed_file_for_strand}"
        " {output.footprint_length_and_percentage}"
        " {output.methylation_level} {output.methylation_level_dict}" 
        " {params.neighbor_left_from_peak_center}"
        " {params.neighbor_right_from_peak_center}"
        " {wildcards.lf} {wildcards.rf}" 
