rule order_fragments_cobinding:
    input:
        fragments_with_assigned_state_verbose_150bp = "cobinding_states/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_cobound_states_verbose_150bp.bed.gz", 
    params:
    output:
        ordered_footprints = "fragments_ordered_cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_fp.tsv",  
        ordered_methylation = "fragments_ordered_cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_methylation.tsv"
    shell:
        "sh scripts/order_footprints_cobinding.sh"
        " {input.fragments_with_assigned_state_verbose_150bp} {output.ordered_footprints}" 
        " {output.ordered_methylation} {wildcards.roi_id}"

rule plot_cobinding_footprints:
    input:
        ordered_footprints = "fragments_ordered_cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_fp.tsv",  
        ordered_methylation = "fragments_ordered_cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_methylation.tsv",
        mnase_data = "mnase_peaks/peak_229.tsv",
        gnuplt_mnase_params = "utils/gnuplot_base_files/mnase_params.gplt",
        gnuplt_footprint_params = "utils/gnuplot_base_files/footprint_params.gplt",
        gnuplt_methylation_params = "utils/gnuplot_base_files/methylation_params.gplt", 
    params:
        modified_lextend = lambda wildcards: str(int(wildcards.lextend) - 150),
        modified_rextend = lambda wildcards: str(int(wildcards.rextend) - 150),
    output:
        footprint_pdf = "plots/cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp.pdf",
        methylation_pdf = "plots/cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation.pdf",
        footprint_plt = "plots/cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp.gplt",
        methylation_plt = "plots/cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation.gplt",
        footprint_mat = "plots/cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp_mat.tsv", 
        methylation_mat = "plots/cobinding/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation_mat.tsv",
 
    shell:
        "sh scripts/plot_footprint_and_methylation_cobinding.sh"
        " {input.ordered_footprints} {input.ordered_methylation}"
        " {input.mnase_data} {input.gnuplt_mnase_params}"
        " {input.gnuplt_footprint_params}"
        " {output.footprint_pdf} {output.methylation_pdf}"
        " {output.footprint_plt} {output.methylation_plt}"
        " {params.modified_lextend} {params.modified_rextend}"
        " {wildcards.lf} {wildcards.rf}"
        " {input.gnuplt_methylation_params}"
        " {output.footprint_mat} {output.methylation_mat}"    

rule order_fragments_cobinding_bedpe:
    input:
        fragments_with_assigned_state_verbose_150bp = "cobinding_states_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_cobound_states_verbose_150bp.bedpe.gz", 
    params:
    output:
        ordered_footprints = "fragments_ordered_cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_fp.tsv",  
        ordered_methylation = "fragments_ordered_cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_methylation.tsv"
    shell:
        "sh scripts/order_footprints_cobinding.sh"
        " {input.fragments_with_assigned_state_verbose_150bp} {output.ordered_footprints}" 
        " {output.ordered_methylation} {wildcards.roi_id}"

rule plot_cobinding_footprints_bedpe:
    input:
        ordered_footprints = "fragments_ordered_cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_fp.tsv",  
        ordered_methylation = "fragments_ordered_cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.ordered_methylation.tsv",
        mnase_data = "mnase_peaks/peak_229.tsv",
        gnuplt_mnase_params = "utils/gnuplot_base_files/mnase_params.gplt",
        gnuplt_footprint_params = "utils/gnuplot_base_files/footprint_params.gplt",
        gnuplt_methylation_params = "utils/gnuplot_base_files/methylation_params.gplt", 
    params:
        modified_lextend = lambda wildcards: str(int(wildcards.lextend) - 150),
        modified_rextend = lambda wildcards: str(int(wildcards.rextend) - 150),
    output:
        footprint_pdf = "plots/cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp.pdf",
        methylation_pdf = "plots/cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation.pdf",
        footprint_plt = "plots/cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp.gplt",
        methylation_plt = "plots/cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation.gplt",
        footprint_mat = "plots/cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.fp_mat.tsv", 
        methylation_mat = "plots/cobinding_bedpe/{sample}_to_{bed}_lf_{lf}_rf_{rf}_extended_left_{lextend}_right_{rextend}_roi_{roi_id}.methylation_mat.tsv",
 
    shell:
        "sh scripts/plot_footprint_and_methylation_cobinding.sh"
        " {input.ordered_footprints} {input.ordered_methylation}"
        " {input.mnase_data} {input.gnuplt_mnase_params}"
        " {input.gnuplt_footprint_params}"
        " {output.footprint_pdf} {output.methylation_pdf}"
        " {output.footprint_plt} {output.methylation_plt}"
        " {params.modified_lextend} {params.modified_rextend}"
        " {wildcards.lf} {wildcards.rf}"
        " {input.gnuplt_methylation_params}"
        " {output.footprint_mat} {output.methylation_mat}"    


