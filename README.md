# SARS-CoV-2 intrahost diversity in hospitalized patients and infected healthcare workers

This repository holds code and data for our paper on SARS-CoV-2 intrahost diversity. Below is a map of the repository structure.

# Overview
--------

    project
    |- README          # The top level description of content. You are here!
    |
    |- data  
    |  |- metadata/   # Sample metadata for the specimens sequenced in this study.
    |  |- reference/  # Reference files used in data analysis pipelines.
    |  |- raw/        # Output data from pipelines, except files >100MB.
    |  |- processed/  # Intermediate files. Used for most analyses.
    |  |- tree/       # TreeTime output for consensus phylogenetic tree.
    |  |- mixing/     # Data from synthetic RNA mixture experiment
    |- scripts/       # Code for analyses presented in manuscript.
    |- pipelines/     # Snakemake files for alignment, consensus calling, and variant calling.
    
  --------

# Notes

Raw sequence reads are available at the Sequence Read Archive, BioProject accession PRJNA682212.

# Contact

If you have questions, please contact the [Lauring Lab](https://lauringlab.wordpress.com/contacts/).
