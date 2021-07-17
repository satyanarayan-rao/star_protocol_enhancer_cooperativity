rule plot_methylation_percentages_on_reads: # related to figure 5a
    input:
        occluded_edges_open = "occluded_edges/{sample}_to_{open_enh}_lf_{lf}_rf_{rf}_occluded.tsv",
        occluded_edges_closed = "occluded_edges/{sample}_to_{closed_enh}_lf_{lf}_rf_{rf}_occluded.tsv",
        
        
    output:
        percentage_methylation_plot = "plots/percentage_methylation/{sample}_to_{open_enh}_vs_{closed_enh}_lf_{lf}_rf_{rf}.png"
    shell:
        "sh scripts/plot_methylation_percentage.sh {input.occluded_edges_open}"
        " {input.occluded_edges_closed} {output.percentage_methylation_plot}"
