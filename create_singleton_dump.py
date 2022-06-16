from glob import glob
from collections import defaultdict 
import logging
import math
import random
import pandas as pd
import numpy as np
import variantkey as vk
import pyfaidx as pyfa
import pyarrow as pa
import pyarrow.parquet as pq
# import pybedtools as pybed
import cyvcf2 as cy
from multiprocessing import Pool
# from ray.util.multiprocessing import Pool
from datetime import date
import os
from IPython.core.debugger import set_trace
from create_sample_list import *

debug = False

num_batches = 270
num_cores = 55
# num_batches = 3
# num_cores = 3
acceptable_chroms = ["chr" + c for c in list(map(str, range(1, 23))) + ["X"]]

if debug:
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
else:
    logging.basicConfig(
        filename = f"logs/singleton_{date.today().strftime('%Y_%m_%d')}.log", 
        filemode ="w", format='%(asctime)s - %(message)s', 
        level=logging.INFO
    )


def chip_lookup_vaf_and_gene(chip_manifest, key):
    if key in chip_manifest.index:
        return chip_manifest.loc[key, ["VAF", "chip_gene", 'CHIP_REF_reads', 'CHIP_ALT_reads']]
    else:
        return np.nan, None, None, None

def bravo_parse_chrom(chrom):

    logging.info(f"now parsing {chrom} variants in bravo vcf")

    bravo_bcf = "bravo-dbsnp-all.bcf"
    if not os.path.exists(bravo_bcf):
        raise IOError(f"{bravo_bcf} does not exist")

    iterable = cy.VCF(bravo_bcf)
    keys = []

    count = 0
    for v in iterable(chrom):

        if count % 1e7 == 0:
            logging.info(f"now parsed {count} variants in {chrom} from bravo vcf")

        CHROM = v.CHROM
        POS = v.start + 1 # cyvcf2 start is 0 based
        REF = v.REF
        ALT = "".join(v.ALT) # is a list
        FILTER = v.FILTER

        # AF = v.INFO.get("AF") 
        # AC = v.INFO.get("AC")

        count += 1

        if FILTER is not None:
            # FILTER = "PASS" # none = pass
            continue

        is_snp = v.is_snp

        if not is_snp and parameters.only_snps:
            continue

        key = vk.variantkey(CHROM, POS, REF, ALT) # THIS WILL CONVERT chrom X to 23! When decoding will be ints

        keys.append(key)

    return keys

def get_bravo_cache_file():
    return "bravo.npy"

def create_bravo_lookup():

    cache_file = get_bravo_cache_file()

    if os.path.exists(cache_file):
        logging.info("now returning bravo lookup from cache")
        # return pd.read_parquet(cache_file, engine = "pyarrow")
        # with open(cache_file, "rb") as f:
        #     res = pickle.load(f)
        #     return res

        res = set(np.load(cache_file).tolist())
        return res
    else:
        logging.info("now creating bravo lookup")

        num_cores = 10
        pool = Pool(num_cores)

        key_collection = pool.map(bravo_parse_chrom, acceptable_chroms)
        keys = set([item for chrom_keys in key_collection for item in chrom_keys])

        pool.close()
        pool.join()

        logging.info("done creating bravo lookup")

        logging.info("now writing index to disk")
        # res.to_parquet(cache_file)
        # with open(cache_file, "wb") as f:
        #     pickle.dump(keys, f)

        np.save(cache_file, keys)

        logging.info("done writing to disk")

        return keys

def bravo_lookup(bravo, key):
    if key in bravo:
        return True
    else:
        # return np.nan
        # return 5e-8 # to be consistent with gnomad based coding in mutect2 germline resource
        return False

def split_into_batches(vcfs, num_batches):

    size_of_batch = math.ceil(len(vcfs) / num_batches)
    
    for i in range(0, len(vcfs), size_of_batch):
        yield vcfs[i:i + size_of_batch]


# def create_interval(chrom, start):
#     return pybed.Interval(chrom, start + 1, start + 1)

def enumerate_test_vcfs():
    test_dir = "test_data/"
    suffix = "*.vcf"
    return glob(test_dir + suffix)

def create_index(bed):
    index = dict()
    for chrom in acceptable_chroms:
        index[chrom] = np.stack((bed[bed.CHROM == chrom].START.values, bed[bed.CHROM == chrom].END.values), axis= 1)

    return index

def check_if_bed_contains_interval(chrom, start, index):
    # 316 us without jit
    # 386 us with jit
    # subtract 1 from VCF start because VCF are 1 based and bed are 0 based
    # this function only makes sense for SNPs
    return np.any((index[chrom][:, 0] <= start - 1) & (start < index[chrom][:, 1]))

def load_bed_file(fname):
    if not os.path.exists(fname):
        IOError(f"{fname} does not exist")

    # return pybed.BedTool(fname)
    df = pd.read_table(fname, header = None, names = ["CHROM", "START", "END"])
    df = df[df.CHROM.isin(acceptable_chroms)]
    return create_index(df)
    
def load_low_complexity():
    fname = "bed/mdust.bed.gz"
    df = load_bed_file(fname)
    return df

def load_seg_dups():
    fname = "bed/genomicSuperDups.bed"
    df = pd.read_table(fname, usecols = ["chrom", "chromStart", "chromEnd"])
    df.columns = ["CHROM", "START", "END"]
    df = df[df.CHROM.isin(acceptable_chroms)]
    return create_index(df)

# def check_if_bed_contains_interval(chrom, start, bed):
    # return bed_file.any_hits(interval) == 1 # returns 1 if there is a hit
    # subtract 1 from VCF start because VCF are 1 based and bed are 0 based
    # this function only makes sense for SNPs
    # return bed[(bed.CHROM == chrom) & (bed.START >= start - 1) & (bed.END <= start)].shape[0] > 0

class vcf_parsing_parameters(object):
   def __init__(
        self, 
        vaf_lower_bound, 
        vaf_upper_bound,
        dp_lower_bound,
        dp_upper_bound,
        only_snps,
        exclude_low_complexity,
        exclude_seg_dups,
        exclude_bravo
    ):
    
        self.vaf_lower_bound        = vaf_lower_bound
        self.vaf_upper_bound        = vaf_upper_bound
        self.dp_lower_bound         = dp_lower_bound
        self.dp_upper_bound         = dp_upper_bound
        self.only_snps              = only_snps
        self.exclude_low_complexity = exclude_low_complexity
        self.exclude_seg_dups       = exclude_seg_dups
        self.exclude_bravo          = exclude_bravo

def parse_vcf(vcf):
    iterable = cy.VCF(vcf)
    if len(iterable.samples) > 1:
        ValueError(f"vcf {vcf} has more than 1 sample")

    sample = iterable.samples[0]
    # logging.info(f"memory address of bravo is {hex(id(bravo))}")

    variants = defaultdict(tuple)
    count = 0

    chip_driver_VAF, chip_driver_gene, chip_ref_reads, chip_alt_reads = chip_lookup_vaf_and_gene(chip_manifest, sample)

    sex = lookup_sex(sex_manifest, sample)

    if chip_driver_gene:
        chip_case = True
    else:
        chip_case = False

    for v in iterable:
        CHROM = v.CHROM
        POS = v.start + 1
        REF = v.REF
        FILTER = v.FILTER

        # set_trace()
        if CHROM not in acceptable_chroms:
            continue

        if FILTER is not None:
            # FILTER = "PASS" # none = pass
            continue

        is_multiallelic = len(v.ALT) > 1
        if is_multiallelic:
            continue

        VAF = np.asscalar(v.format("AF")[0][0])
        if sex == 1 and CHROM == "chrX" : # half VAF if male for chrom X variants
            VAF = VAF / 2.0

        if VAF >= parameters.vaf_upper_bound: # only include variants with VAF < 30%
            continue

        DP = v.INFO.get("DP") 
        if DP < parameters.dp_lower_bound or DP > parameters.dp_upper_bound: # only include sites with depth >= 10
            continue

        # key = (CHROM, POS, REF, ALT, )
        ALT = "".join(v.ALT) # is a list

        sequence_context = get_reference_sequence(reference, CHROM, POS) # reference is a global variable

        if "N" in sequence_context:
            continue

        is_snp = v.is_snp and len(REF) == 1 #cyvcf2 seems to only check length of alt allele?

        if not is_snp and parameters.only_snps:
            continue

        # if sequence_context[4] != REF:
        #     set_trace()

        if is_snp:
            trinucleotide_substitution = f"{sequence_context[3]}{REF}>{ALT}{sequence_context[5]}"
        else:
            trinucleotide_substitution = None
        

        key = vk.variantkey(CHROM, POS, REF, ALT) # THIS WILL CONVERT chrom X to 23! When decoding will be ints

        variant_is_in_bravo = bravo_lookup(bravo, key)

        if parameters.exclude_bravo and variant_is_in_bravo:
            continue

        if parameters.exclude_low_complexity and check_if_bed_contains_interval(CHROM, POS, low_complexity):
            # logging.info(f"now excluding {CHROM}:{POS} because it is a low complexity")
            continue

        if parameters.exclude_seg_dups and check_if_bed_contains_interval(CHROM, POS, seg_dups):
            # logging.info(f"now excluding {CHROM}:{POS} because it is a seg dup")
            continue


        # Mutect2 bug with integer overflowhttps://gatkforums.broadinstitute.org/gatk/discussion/comment/60647
        # here I set to missing if less than math.pow(2, 15), i.e., limit of a signed int16
        MPOS = np.asscalar(v.format("MPOS")[0][0])      
        if MPOS < -32768:
            MPOS = np.nan

        count += 1

        # if count % 1e2 == 0:
        #     logging.info(f"now completed {count} variants")

        value = (sample,
                # FILTER,
                sequence_context,
                trinucleotide_substitution,
                v.is_transition,
                is_snp,
                DP, 
                v.INFO.get("ECNT"), 
                v.INFO.get("IN_PON"), 
                v.INFO.get("TLOD"), 
                v.INFO.get("STR"), 
                v.INFO.get("P_GERMLINE"), 
                v.INFO.get("P_CONTAM"),
                v.INFO.get("POP_AF"),
                v.format("AD")[0][1],
                VAF,
                np.asscalar(v.format("MMQ")[0][0]), # mean mapping quality 
                MPOS,
                np.asscalar(v.format("F1R2")[0][1]), 
                np.asscalar(v.format("F2R1")[0][1]), 
                chip_driver_VAF,
                chip_driver_gene,
                chip_ref_reads,
                chip_alt_reads,
                chip_case
                )
        if key not in variants:
            variants[key] = value
        else:
            raise ValueError("variant is already in single sample VCF")


    logging.info(f"found {count} variants in {vcf}")
    return  variants

def walk_batch(vcfs, num_cores = num_cores):
    
    pool = Pool(num_cores)

    variant_tables = pool.map(parse_vcf, vcfs)
    pool.close()
    pool.join()
    # variant_tables = [parse_vcf(v) for v in vcfs]

    return find_singletons(variant_tables)

def get_reference_file():
    path = "input/Homo_sapiens_assembly38.fasta"
    if not os.path.exists(path):
        raise IOError(f"{path} does not exist")
    return path

def get_reference_sequence(reference, chrom, start, kmer_size = 9):

    # inspiration from https://github.com/carjed/helmsman/blob/master/util.py
    if kmer_size % 2 == 0 or kmer_size <= 1:
        raise ValueError(f"kmer_size of {kmer_size} is not an odd number > 1")

    pad = kmer_size // 2
    sequence = reference[chrom][(start - pad - 1):(start + pad)]

    return str(sequence)

def find_singletons(variant_tables, grand = False):
    
    if grand:
        logging.info(f"now finding singletons out of {len(variant_tables)} batches")
    else:
        logging.info(f"now finding singletons out of {len(variant_tables)} vcfs")

    not_singletons = []
    singletons = defaultdict(tuple)
    count_of_singletons = 0
    count = 0

    # set_trace()
    for table in variant_tables:
        count += 1
        logging.info(f"now on table {count} out of {len(variant_tables)}")
        singletons = {k: singletons[k] if k in singletons else table[k] for k in singletons.keys() ^ table.keys()}


    logging.info(f"found {len(list(singletons.keys()))} variants that are singletons")
    return singletons

def save_to_parquet(singletons, test = False):

    logging.info("now constructing dataframe from singleton list")

    parsed_variants = []

    count = 0
    total = len(list(singletons.keys()))

    collector = dict()

    schema = pa.schema({
        "key"                        : pa.uint64(),
        "CHROM"                      : pa.int32(),
        "POS"                        : pa.int32(),
        "REF"                        : pa.string(),
        "ALT"                        : pa.string(),
        "NWD_ID"                     : pa.string(),
        # "FILTER"                   : FILTER,
        "sequence_context"           : pa.string(),
        "trinucleotide_substitution" : pa.string(),
        "is_transition"              : pa.bool_(),
        "is_snp"                     : pa.bool_(),
        "DP"                         : pa.int32(),
        "ECNT"                       : pa.int32(),
        "IN_PON"                     : pa.bool_(),
        "TLOD"                       : pa.float64(),
        "STR"                        : pa.bool_(),
        "P_GERMLINE"                 : pa.float64(),
        "P_CONTAM"                   : pa.float64(),
        "POP_AF"                     : pa.float32(),
        # "bravo_AF"                   : pa.float32(),
        "AD"                         : pa.int32(),
        "VAF"                         : pa.float32(),
        "MMQ"                        : pa.int32(),
        "MPOS"                       : pa.int32(), # yes, this can be negative or positive
        "F1R2"                       : pa.int32(),
        "F2R1"                       : pa.int32(),
        "chip_driver_VAF"            : pa.float32(),
        "chip_driver_gene"           : pa.string(),
        "chip_ref_reads"             : pa.int32(),
        "chip_alt_reads"             : pa.int32(),
        "chip_case"                  : pa.bool_()
    })

    for name in schema.names:
        collector[name] = []

    for key, item in singletons.items():
        chrom, pos, refalt = vk.decode_variantkey(key)
        ref, alt, _, __ = vk.decode_refalt(refalt)

        collector["key"].append(key)
        collector["CHROM"].append(chrom)
        collector["POS"].append(pos)
        collector["REF"].append(ref)
        collector["ALT"].append(alt)

        for name, val in zip(schema.names[5:], item):
            collector[name].append(val)

        count += 1
        if count % int(1e5) == 0:
            logging.info(f"now on {count} out of {total} singletons")


    logging.info("now converting to pyarrow")
    # parsed_variants = pa.Table.from_pydict(
    #    collector, schema  = schema
    # )
    parsed_variants = pa.Table.from_pandas(
       pd.DataFrame(collector), schema  = schema
    )

    logging.info("now writing to parquet")

    if test:
        pq.write_table(parsed_variants, "test_data/singletons_{}.parquet".format(date.today().strftime("%Y_%m_%d")))
    else:
        pq.write_table(parsed_variants, "singletons_{}.parquet".format(date.today().strftime("%Y_%m_%d")))
        pq.write_table(parsed_variants.drop(["key"]), "singletons_spark_{}.parquet".format(date.today().strftime("%Y_%m_%d")))

    logging.info("done writing.")
    return parsed_variants

        
if __name__ == "__main__":
    random.seed(1)
    variants = defaultdict(tuple)
    # chip_manifest = create_chip_case_lookup()
    chip_manifest = create_chip_single_driver_case_lookup()
    sex_manifest = create_sex_lookup()
    vcf_files = enumerate_vcfs()
    # vcf_files = enumerate_test_vcfs()
    batches = split_into_batches(vcf_files, num_batches)

    logging.info("now loading reference now")
    reference = pyfa.Fasta(get_reference_file())
    # load global variables to be referenced in parse_vcf function
    bravo = create_bravo_lookup()
    # logging.info(f"memory address of bravo is {hex(id(bravo))}")
    low_complexity = load_low_complexity()
    seg_dups = load_seg_dups()

    parameters = vcf_parsing_parameters(
        vaf_lower_bound        = 0.0, 
        vaf_upper_bound        = 0.35,
        dp_lower_bound         = 25,
        dp_upper_bound         = 100,
        only_snps              = True,
        exclude_low_complexity = True,
        exclude_seg_dups       = True,
        exclude_bravo          = True
    )

    logging.info(f"now starting on {num_batches} batches")

    batch_tables = []
    count = 1
    for batch in batches:
        logging.info(f"now on batch {count} out of {num_batches}")
        batch_tables.append(walk_batch(batch))
        count = count + 1

    singletons = find_singletons(batch_tables)
    logging.info("now finding singletons among batches")

    df = save_to_parquet(singletons)
