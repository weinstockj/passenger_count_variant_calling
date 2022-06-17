# Passenger count variant calling

Below is the passenger count variant calling procedure that was developed for [Clonal hematopoiesis is driven by aberrant activation of TCL1A](https://www.biorxiv.org/content/10.1101/2021.12.10.471810v1). 
For questions, contact [Josh Weinstock](jweinstk@umich.edu). 

The procedure assumes that [Mutect2](https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2.wdl) has already been run in 'tumor-only' mode across the entire genome
for samples of interest. For Mutect2 documentation, please refer to the above linked documentation
describing a WDL workflow for running Mutect2. The script assumes the suffix of these files
is `*-filtered.vcf.gz`, though changing this to bcf will also work. Storing these files
as raw (rather than bgzipped compressed or bcf) VCFs is strongly discouraged. 

Currently, this directory is harded coded in `create_sample_list.py`, and should be modified
to the appropriate user directory. In addition, these scripts assume that a sample manifest
with metadata on CHIP carriers has already been created. This file (hard-coded at the moment) should be a tab separated file with the following column headers and one row per sample:
 
 1. Sample (unique sample identifier)
 2. INFERRED\_SEX (a column coded as 1/2 for genotype inferred sex)
 3. Gene (a column indicating the mutated driver gene)
 4. haschip (a binary column where CHIP carriers are coded as 1)
 5. STUDY (indicating the cohort the sample is from)
 
In addition, the script assumes a second variant level metadata file (also hard-coded) in tab separated form. This file should contain one row per driver variant. This file should have the following headers:

 1. Sample
 2. Gene
 3. AD (allelic depths of REF and ALT reads coded as {REF},{ALT})
 4. VAF (alt / (ref + alt))

This file is primarly used to subset the CHIP carriers to those with a single driver. 

For the variant calling itself, some secondary files are needed:
 1. "bravo-dbsnp-all.bcf" A BCF of the TOPMed Bravo sitelist of pass variants to exclude. See [here](https://bravo.sph.umich.edu/freeze8/hg38/downloads) for how to download this file. 
 2. A file of low complexity regions (bed/mdust.bed/gz), available from UCSC
 3. A file of segmental dupliations (bed/genomicSuperDups.bed), availahle from UCSC

 ## Python (>= 3.6) dependencies
 1. pandas
 2. numpy
 3. [variantkey](https://github.com/Genomicsplc/variantkey)
 4. pyfaidx
 5. pyarrow
 6. cyvcf2
 7. logging

## Output
The output of `create_singleton_dump.py` is an Apache `parquet` file with one row
per variant. Several variant quality metrics are included in the output to 
facilitate downstream filtering. 

## Helper scripts
Helper scripts are provided for downloading the Bravo site list and installing variantkey. 

## Notes on portability
This analysis has only been tested on Ubuntu 18.04 with Python (>=3.6).  

## License
This code is dual-licensed. See the [license](LICENSE.md) for further details. 

