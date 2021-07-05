rule create_roi_info_dict:
    input:
        regions_bed = "input_bed/{bed}.bed"
    output:
        regions_metadata = "regions_metadata/{bed}.pkl"
    shell:
        "cat {input.regions_bed} | python scripts/generate_roi_metadata.py"
        " {output.regions_metadata}"
