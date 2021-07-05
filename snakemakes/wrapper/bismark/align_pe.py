"""Snakemake wrapper for aligning using bismark """
__author__ = "Satyanarayan Rao"
__copyright__ = "Copyright 2019, Satyanarayan Rao"
__email__ = "satyanarayan.rao@cuanschutz.edu"
__license__ = "MIT"

from snakemake.shell import shell 
import os.path 

log = snakemake.log_fmt_shell()

# Check if input length is 2
n = len(snakemake.input)
assert n == 2, "Input must have two files. Got: %r. " % n

shell(
    "(bismark"
    " {snakemake.params.extra}"
    " {snakemake.params.genome}"
    " {snakemake.params.basename}"
    " -o {snakemake.params.out_dir}"
    " -1 {snakemake.input[0]}"
    " -2 {snakemake.input[1]})" 
    " {log}"
)
