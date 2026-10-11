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

rulename = "write_mhc1_peptides"
use rule write_mhc1_peptides from neoscoreT as write_mhc1_peptides with:
    input:
        peps = config["input_peptide_table"]
    output:
        splits = temp(expand(os.path.join("results", "netmhcpan", "mhc1peps_len{length}_split{mhc1_split}.txt"), length=LENGTHS, mhc1_split=SPLITS))

rulename = "run_netmhcpan"
use rule run_netmhcpan from neoscoreT as run_netmhcpan with:
    input:
        peps = os.path.join("results", "netmhcpan", "mhc1peps_len{length}_split{mhc1_split}.txt")
    output:
        raw_netmhc_out = temp(os.path.join("results", "netmhcpan", "raw_netmhc_output_{mhc1_allele}_len{length}_split{mhc1_split}.txt"))

rulename = "process_netmhcpan_output"
use rule process_netmhcpan_output from neoscoreT as process_netmhcpan_output with:
    input:
        raw_netmhc_out = rules.run_netmhcpan.output.raw_netmhc_out
    output:
        netmhc_pred = temp(os.path.join("results", "netmhcpan", "netmhc_predictions_{mhc1_allele}_len{length}_split{mhc1_split}.csv"))

rulename = "merge_netmhcpan_predictions"
use rule merge_netmhcpan_predictions from neoscoreT as merge_netmhcpan_predictions with:
    input:
        main_table = config["input_peptide_table"],
        peps = expand(os.path.join("results", "netmhcpan", "netmhc_predictions_{mhc1_allele}_len{length}_split{mhc1_split}.csv"), mhc1_allele = ["A"], length = LENGTHS, mhc1_split = SPLITS)
    output:
        merged = temp(os.path.join("results","netmhcpan","netmhc_predictions_merged.csv")),
        full_predictions = os.path.join("results","netmhcpan","netmhcpan_predictions_full_table.csv"),
        combined_table = temp(os.path.join("results", "netmhcpan", "annotated_summary_table.tsv"))
