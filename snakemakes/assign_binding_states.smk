rule assign_binding_states:
    input:
        footprint_length_and_percentage = \
          "methylation_levels_around_roi/{sample}_to_{bed}_lf_{lf}_rf_{rf}_footprint_fraction_and_length.tsv",
        methylation_level_dict = \
          "methylation_levels_around_roi/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_level.pkl",
        methylation_matrix = "flank_footprint_matrix/{sample}_to_{bed}_lf_{lf}_rf_{rf}_methylation_matrix.tsv" ,
        occluded_pkl = "occluded_edges/{sample}_to_{bed}_lf_{lf}_rf_{rf}_occluded.pkl",
        footprint_dict = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.pkl",
        regions_metadata = "regions_metadata/{bed}.pkl"
    params:
    output:
        binding_state_for_fragment = "binding_states/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_binding_state.tsv" 
    shell:
        "python scripts/assign_binding_states.py {input.footprint_length_and_percentage}"
        " {input.methylation_level_dict} {input.methylation_matrix}"
        " {input.occluded_pkl} {input.footprint_dict} {input.regions_metadata}"
        " {output.binding_state_for_fragment} {wildcards.lf} {wildcards.rf}"
       
rule extend_binding_state_fragments:
    input:
        binding_state_for_fragment = "binding_states/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_binding_state.tsv",
        footprint_dict = "fragments_spanning_flanks/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}.pkl",
        regions_metadata = "regions_metadata/{bed}.pkl"
    params:
        
    output:
        extended_fragments = \
          "binding_states_fragments_extended/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}.tsv", 
        verbose = \
          "binding_states_fragments_extended/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}.verbose.tsv", 
    shell:
        "python scripts/extend_footprint.py {input.binding_state_for_fragment}"
        " {input.footprint_dict} {input.regions_metadata}"
        " {wildcards.lf} {wildcards.rf} {wildcards.lextend} {wildcards.rextend}"
        " {output.extended_fragments} {output.verbose}"
