rule order_fragments_single_binding:
    input:
        verbose = \
          "binding_states_fragments_extended/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}.verbose.tsv", 
    params:
    output:
        ordered_footprints = \
          "fragments_ordered_single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_fp.tsv",
        ordered_methylation = \
          "fragments_ordered_single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_methylation.tsv",
    shell:
        "sh scripts/order_footprints_single_binding.sh"
        " {input.verbose} {output.ordered_footprints}" 
        " {output.ordered_methylation} {wildcards.roi_id}"

rule plot_footprints:
    input:
        ordered_footprints = \
          "fragments_ordered_single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_fp.tsv",
        ordered_methylation = \
          "fragments_ordered_single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_methylation.tsv",
        mnase_data = "mnase_peaks/peak_229.tsv",
        gnuplt_mnase_params = "utils/gnuplot_base_files/mnase_params.gplt",
        gnuplt_footprint_params = "utils/gnuplot_base_files/footprint_params.gplt",
        gnuplt_methylation_params = "utils/gnuplot_base_files/methylation_params.gplt",        
    params:
        
    output:
        footprint_pdf = "plots/single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp.pdf",
        methylation_pdf = "plots/single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation.pdf",
        footprint_plt = "plots/single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp.gplt",
        methylation_plt = "plots/single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation.gplt",
        footprint_mat = "plots/single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp_mat.tsv", 
        methylation_mat = "plots/single_binding/{sample}_to_{bed}_spanning_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation_mat.tsv",
    shell:
        "sh scripts/plot_footprint_and_methylation.sh"
        " {input.ordered_footprints} {input.ordered_methylation}"
        " {input.mnase_data} {input.gnuplt_mnase_params}"
        " {input.gnuplt_footprint_params}"
        " {output.footprint_pdf} {output.methylation_pdf}"
        " {output.footprint_plt} {output.methylation_plt}"
        " {wildcards.lextend} {wildcards.rextend}"
        " {wildcards.lf} {wildcards.rf}"
        " {input.gnuplt_methylation_params}"
        " {output.footprint_mat} {output.methylation_mat}"
