rule assign_cobinding_states:
    input:
        extended_fragments = "extended_dsmf_reads_to_peak_pairs/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}.bed.gz", 
        verbose = "extended_dsmf_reads_to_peak_pairs/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_verbose.bed.gz", 
    params:
    output: 
        fragments_with_assigned_state = "cobinding_states/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_cobound_states.bed.gz", 
        fragments_with_assigned_state_verbose = "cobinding_states/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_cobound_states_verbose.bed.gz", 
        fragments_with_assigned_state_verbose_150bp = "cobinding_states/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_cobound_states_verbose_150bp.bed.gz", 
        
        
    shell:
        "sh scripts/assign_cobinding_states.sh"
        " {input.extended_fragments} {input.verbose}"
        " {wildcards.lf} {wildcards.rf} {wildcards.lextend} {wildcards.rextend}"  
        " {output.fragments_with_assigned_state} {output.fragments_with_assigned_state_verbose}" 
        " {output.fragments_with_assigned_state_verbose_150bp}"
