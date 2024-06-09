import os
import yaml

['AMD_Epyc' 'Ampere_Q8030' 'Intel_SPR_DDR' 'Intel_SPR_HBM' 'archer_rome' 'snellius-genoa' 'snellius-rome']
['AMD Bergamo', 'Ampere Altra', 'Intel SPR DDR', 'Intel SPR HBM', 'AMD Rome 7742', 'AMD Genoa', 'AMD Rome 7H12']

platform_cores = {
        "AMD Bergamo": 256,
        "Ampere Altra": 80,
        "Intel SPR DDR": 112,
        "Intel SPR HBM": 112,
        "AMD Rome 7742": 128,
        "AMD Genoa": 192,
        "AMD Rome 7H12": 128
        }


translate_exp = {
        "fixed": "S1",
        "Fixed": "S1",
        "relative": "W1",
        "Relative-scale": "W1-scale",
        "fixed-large": "S2",

        "cube-fixed": "S1",
        "cube-relative": "W1",
        "cube-fixed-large": "S2",

        "fixed-empty": "S1-no-iter",
        "relative-empty": "W1-no-iter",
        "fixed-large-empty": "S2-no-iter",

        "cube-fixed-empty": "S1-no-iter",
        "cube-relative-empty": "W1-no-iter",
        "cube-fixed-large-empty": "S2-no-iter",
        
        "fixed-norbc": "S1-no-RBC",
        "fixed-large-norbc": "S2-no-RBC",
        "relative-norbc": "W1-no-RBC",
        "atomic-block": "atomic",
        None: None
        }

translate_exp_to_internal = {
        "S1": "fixed",
        "W1": "relative",
        "S2": "fixed-large",
        
        None: None
        }

translate_platform = {
        "AMD_Epyc": "AMD Bergamo",
        "Ampere_Q8030": "Ampere Altra",
        "Intel_SPR_DDR": "Intel SPR DDR",
        "Intel_SPR_HBM": "Intel SPR HBM",
        "archer_rome": "AMD Rome 7742",
        "snellius-genoa": "AMD Genoa",
        "snellius-rome": "AMD Rome 7H12",
        "snellius-rome-likwid": "rome-likwid",
        None: None
        }

translate_platform_rev = {
        "AMD Bergamo": "AMD_Epyc" ,
        "Ampere Altra": "Ampere_Q8030",
        "Intel SPR DDR": "Intel_SPR_DDR",
        "Intel SPR HBM": "Intel_SPR_HBM",
        "AMD Rome 7742": "archer_rome",
        "AMD Genoa": "snellius-genoa",
        "AMD Rome 7H12": "snellius-rome",
        "rome-likwid": "snellius-rome-likwid",
        None: None
        }

names_for_hw_events = {
        "instructions": "Instructions",
        "cycles": "Cycles",
        "l1_accesses": "L1 request",
        "l1_misses": "L1 misses",
        "l1_load_accesses": "L1 loads",
        "l1_load_misses": "L1 load misses",
        "l1_store_accesses": "L1 stores",
        "l1_store_misses": "L1 store misses",
        "CPI": "CPI",
        "data volume": "Data Volume",
        "runtime": "Runtime [s]",
        "mem bandwidth": "Bandwidth",
        "loads" : "Memory Loads",
        "stores" : "Memory Stores",
        "LLC_store_accesses": "LLC stores",
        "LLC_load_accesses": "LLC loads",
        "LLC_load_misses": "LLC load misses",
        "LLC_store_misses": "LLC store misses",
        "LLC_misses": "LLC misses",
        "LLC_accesses": "LLC requests",
        "total": "total",
        "mpi": "mpi",
        "comp": "comp",
        }

translate_hw_events = {
        "ex_ret_instr": names_for_hw_events["instructions"],
        "instructions": names_for_hw_events["instructions"],
        "INSTR_RETIRED_ANY": names_for_hw_events["instructions"],
        "RETIRED_INSTRUCTIONS": names_for_hw_events["instructions"],
        "INST_RETIRED": names_for_hw_events["instructions"],
        "CPU_CLOCKS_UNHALTED": names_for_hw_events["cycles"],
        "cpu-cycles": names_for_hw_events["cycles"],
        "CPU_CLK_UNHALTED_CORE": names_for_hw_events["cycles"],
        "CPU_CYCLES": names_for_hw_events["cycles"],
        "DATA_CACHE_ACCESSES": names_for_hw_events["l1_accesses"],
        "DATA_CACHE_REFILLS_ALL": names_for_hw_events["l1_misses"],
        "l1_data_cache_fills_all": names_for_hw_events["l1_misses"],
        'all_data_cache_accesses': names_for_hw_events["l1_accesses"],
        "L1-dcache-load-misses": names_for_hw_events["l1_load_misses"],
        "L1-dcache-loads": names_for_hw_events["l1_load_accesses"],
        "l2_cache_accesses_from_dc_misses": names_for_hw_events["l1_misses"],
        "all_dc_accesses": names_for_hw_events["l1_accesses"],
        "l1d_cache": names_for_hw_events["l1_accesses"],
        "l1d_cache_refill": names_for_hw_events["l1_misses"],
        "HBM data volume": names_for_hw_events["data volume"],
        "DDR data volume": names_for_hw_events["data volume"],
        "Runtime unhalted": names_for_hw_events["runtime"],
        "runtime": names_for_hw_events["runtime"],
        "Memory data volume": names_for_hw_events["data volume"],
        "MEM_INST_RETIRED_ALL_LOADS": names_for_hw_events["loads"],
        "MEM_INST_RETIRED_ALL_STORES": names_for_hw_events["stores"],
        "LD_SPEC" : names_for_hw_events["loads"],
        "ST_SPEC" : names_for_hw_events["stores"],
        'l2_request_g1.all_no_prefetch' : names_for_hw_events["l1_misses"],
        'l2_request_g1' : names_for_hw_events["l1_misses"],
        "LLC-load-misses": names_for_hw_events["LLC_load_misses"],
        "LLC-store-misses": names_for_hw_events["LLC_store_misses"],
        "LLC-loads": names_for_hw_events["LLC_load_accesses"],
        "LLC-stores": names_for_hw_events["LLC_store_accesses"],
        "l3d_cache": names_for_hw_events["LLC_accesses"],
        "l3d_cache_refill": names_for_hw_events["LLC_misses"],
        "CYCLES_NOT_IN_HALT": names_for_hw_events["cycles"],
        "perf::CYCLES": names_for_hw_events["cycles"],
        "perf::PERF_COUNT_HW_CACHE_L1D:MISS": names_for_hw_events["l1_misses"],
        "perf::PERF_COUNT_HW_CACHE_L1D:ACCESS": names_for_hw_events["l1_accesses"],
        "execution": names_for_hw_events["total"],
        "mpi": names_for_hw_events["mpi"],
        "comp": names_for_hw_events["comp"],
        }


def load_experiment_meta_data(experiment_dictionary):
    d = {}
    metafile = ""
    if os.path.isfile(f"{experiment_dictionary}/meta.yml"):
        metafile = f"{experiment_dictionary}/meta.yml"

    elif os.path.isfile(f"{experiment_dictionary}/meta.yaml"):
        metafile = f"{experiment_dictionary}/meta.yaml"
    else:
        return None

    with open(metafile) as f:
        d = yaml.load(f, Loader=yaml.FullLoader)
    return d
