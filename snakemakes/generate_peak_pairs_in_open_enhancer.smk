rule get_peak_pairs:
    input:
        mnase_peaks = "input_bed/mnase_peaks_in_open_and_closed_enhancers.bed", 
        mapped_dsmf = "fragments_spanning_flanks/suppressed_merged_S2_to_mnase_peaks_in_open_and_closed_enhancers_spanning_lf_15_rf_15.bed.gz", 
        genome_size_file = "metadata/dm3.chrom.sizes"
    params:
        distance_th = 30,
        flank_from_peak = 15
    output:
        mnase_peak_pairs = "input_bed/mnase_peak_pairs_in_open_enhancers.bed"
    shell: 
        "sh scripts/generated_mnase_peak_pairs.sh {input.mnase_peaks}"
        " {input.mapped_dsmf} {input.genome_size_file} {params.distance_th} {params.flank_from_peak}"
        " {output.mnase_peak_pairs}"

