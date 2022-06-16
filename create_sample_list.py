import pandas as pd
import logging
import random
from glob import glob
import os
from IPython.core.debugger import set_trace

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

# directory to list of Mutect2 whole genome calls
vcf_dir = "../mutect2_whole_genome_output/"
def enumerate_vcfs(number_of_controls = 8000):
    suffix = "*-filtered.vcf.gz"
    paths = glob(vcf_dir + suffix)
    logging.info("now finding chip cases in file paths by searching for NWD_ID substring") 
    samples = create_samples_to_include(number_of_controls)
    paths_include = []
    for vcf_path in paths:
        nwd_id = os.path.basename(vcf_path).split(".")[0]
        if nwd_id in samples:
            paths_include.append(vcf_path)

    return paths_include

def create_chip_single_driver_case_lookup():
    """
    Returns a list of chip cases with a single driver
    Expects the following column names in the tsv (tab sep value)
    Sample
    Gene
    AD
    VAF
    """
    path = get_chip_variant_path()
    manifest = pd.read_table(path)
    manifest.rename(columns = {"Sample" : "NWD_ID", "Gene" : "chip_gene"}, inplace = True)
    # TODO: this logic does not work for multi-allelics - just taking first value here
    manifest["CHIP_REF_reads"] = manifest.AD.apply(lambda x: int(x.split(",")[0]))
    manifest["CHIP_ALT_reads"] = manifest.AD.apply(lambda x: int(x.split(",")[1]))
    # manifest["VAF"] = manifest.VAF.apply(lambda x: float(x.split(",")[0]))
    sample_counts = manifest.NWD_ID.value_counts() == 1
    samples_one_driver = sample_counts[sample_counts].index.to_list() 

    columns_to_select = ['NWD_ID', 'VAF', 'chip_gene', 'CHIP_REF_reads', 'CHIP_ALT_reads']
    subset_manifest = manifest[manifest.NWD_ID.isin(samples_one_driver)][columns_to_select]

    logging.info(f"case manifest has {subset_manifest.shape[0]} rows")

    return subset_manifest.set_index("NWD_ID")

def cohorts_to_include():
    # list TOPMed cohorts to include here (e.g.):
    # cohorts = [
    #     'whi',
    #     'mesa',
    #     'jackson',
    #     'atrial_fib',
    #     'ccdg-afib',
    #     'aric',
    #     'amish',
    #     'chs',
    #     'framingham',
    #     'genestar',
    #     'genoa',
    #     'copdgene',
    #     'biome-cad',
    #     'hemophilia',
    #     'hchs-sol',
    #     'cardia',
    #     'eclipse',
    #     'thrv',
    #     'reds-iii-scd',
    #     'sarp',
    #     'mayo-vte',
    #     'asthma_cr',
    #     'walk-phasst'
    # ]

    return cohorts

def create_samples_to_include(number_of_controls = 8000):
    path = get_chip_sample_path()
    manifest = pd.read_table(path)
    manifest.rename(columns = {"Sample" : "NWD_ID", "Gene" : "chip_gene"}, inplace = True)

    cohorts = cohorts_to_include()

    logging.info(f"now including {len(cohorts)} cohorts")

    chip_cases = manifest[manifest.haschip == 1].NWD_ID.values.tolist()
    chip_controls = manifest[manifest.haschip == 0].NWD_ID.values.tolist()
    cohort_samples = manifest[manifest.STUDY.isin(cohorts)].NWD_ID.values.tolist()

    # set_trace()
    chip_controls = random.sample(chip_controls, k = number_of_controls) # without replacement
    logging.info(f"returning {len(chip_cases)} cases and {len(chip_controls)} controls") 

    total_samples = chip_cases + chip_controls + cohort_samples
    total_samples = set(total_samples)

    exclude_samples = set(get_samples_to_exclude())

    total_samples = total_samples - exclude_samples

    return list(total_samples)

def get_samples_to_exclude_path():
    path = "input/TOPMed_CHIP_duplicates.txt"
    if not os.path.exists(path):
        IOError(f"chip path {path} does not exist")
    else:
        return path

def get_samples_to_exclude():
    sample_exlude_path = get_samples_to_exclude_path()
    exclude = []
    with open(sample_exlude_path, "r") as f:
        for line in f:
            exclude.append(line.strip())
    return exclude

def get_chip_variant_path():
    path = "input/TOPMed_variant_level_CHIPcalls_with_covariates_2020_08_31.tsv"
    if not os.path.exists(path):
        IOError(f"chip path {path} does not exist")
    else:
        return path

def get_chip_sample_path():
    path = "input/TOPMed_CHIPcalls_with_covariates_2020_08_31.tsv.gz"
    if not os.path.exists(path):
        IOError(f"chip path {path} does not exist")
    else:
        return path

def create_sex_lookup():
    path = get_chip_sample_path()
    manifest = pd.read_table(path)
    manifest.rename(columns = {"Sample" : "NWD_ID"}, inplace = True)

    manifest = manifest[["NWD_ID", "INFERRED_SEX"]]

    manifest = manifest.set_index("NWD_ID")

    return manifest

def lookup_sex(manifest, sample):
    return manifest.loc[sample, "INFERRED_SEX"]
