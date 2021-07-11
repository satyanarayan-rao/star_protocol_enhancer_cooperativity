rule make_percentage_occuppancy_table:
    input:
        binding_state_for_fragment = "binding_states/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_binding_state.tsv",
        peak_to_region_association = "input_bed/{bed}.bed"
      
    params:
    output:
        occupancy_long_table = "occupancy_table/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_occupancy.tsv"
    shell:
        "cat {input.binding_state_for_fragment}"
        " | python scripts/prepare_occupancy_table.py {input.peak_to_region_association}"
        " > {output.occupancy_long_table}" 
    
rule plot_occupancy_boxplot:
    input:
        occupancy_long_table = "occupancy_table/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_occupancy.tsv"
    params:
    output:
        occupancy_boxplot = "plots/occupancy_boxplot/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_occupancy.png"
    shell:
        "cat {input.occupancy_long_table} | Rscript scripts/occupancy_boxplot.R"
        " {output.occupancy_boxplot}" 
    
