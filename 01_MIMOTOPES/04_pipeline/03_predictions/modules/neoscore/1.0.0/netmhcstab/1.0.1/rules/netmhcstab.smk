import scripts.utils as utils
import pandas as pd
#from scripts.utils.path_builder import *
#writer = Writer(config["benchmark"], directory = "results", pipeline = "neoscore", module = "netmhcpan")
import importlib.util
from pathlib import Path

# Define paths using os.path.join for absolute consistency
base_dir = os.getcwd()
pipeline_dir = os.path.join(base_dir, "pipelines", "neoscore", "1.0.0")
utils_path = os.path.join(pipeline_dir, "scripts", "utils.py")

# Import the module dynamically
spec = importlib.util.spec_from_file_location("local_utils", utils_path)
local_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(local_utils)

LENGTHS = local_utils.get_pep_lengths(config["input_peptide_table"])
NUM_SPLITS = local_utils.set_num_splits(config["input_peptide_table"])
SPLITS = local_utils.generate_splits(NUM_SPLITS)

rulename = "run_netmhcstab"
use rule run_netmhcstab from neoscoreT as run_netmhcstab with:
    input:
        peps = os.path.join("results", "netmhcpan", "mhc1peps_len{length}_split{mhc1_split}.txt")
    output:
        raw_netmhcstab_out = temp(os.path.join("results", "netmhcpan", "raw_netmhcstab_output_{mhc1_allele}_len{length}_split{mhc1_split}.txt"))

rulename = "process_netmhcstabpan_output"
use rule process_netmhcstabpan_output from neoscoreT as process_netmhcstabpan_output with:
    input:
        raw_netmhc_out = os.path.join("results", "netmhcpan", "raw_netmhcstab_output_{mhc1_allele}_len{length}_split{mhc1_split}.txt")
    output:
        netmhcstab_pred = temp(os.path.join("results", "netmhcpan", "netmhcstab_predictions_{mhc1_allele}_len{length}_split{mhc1_split}.csv"))

rulename = "merge_netmhcstabpan_predictions"
use rule merge_netmhcstabpan_predictions from neoscoreT as merge_netmhcstabpan_predictions with:
    input:
        main_table = os.path.join("results", "netmhcpan", "annotated_summary_table.tsv"),
        peps = expand(os.path.join("results", "netmhcpan", "netmhcstab_predictions_{mhc1_allele}_len{length}_split{mhc1_split}.csv"), mhc1_allele = ["A"], length = LENGTHS, mhc1_split = SPLITS)
    output:
        merged = temp(os.path.join("results","netmhcpan","netmhcstab_predictions_merged.csv")),
        full_predictions = os.path.join("results","netmhcpan","netmhcpanstab_predictions_full_table.csv"),
        combined_table = os.path.join("results", "netmhcpan", "predictions_table.csv")
