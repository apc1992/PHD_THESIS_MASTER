import scripts.utils as utils
import pandas as pd
#from scripts.utils.path_builder import *
#writer = Writer(config["benchmark"], directory = "results", pipeline = "neoscore", module = "antigen_scoring")
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

rulename = "pep_to_fasta"
use rule pep_to_fasta from neoscoreT as pep_to_fasta with:
    input:
        peps = os.path.join("results", "netmhcpan", "mhc1peps_len{length}_split{mhc1_split}.txt")
    output:
        peps_fa = temp(os.path.join("results", "neo_score", "mhc1peps_len{length}_split{mhc1_split}.fa"))

rulename = "iedb_alignment"
use rule iedb_alignment from neoscoreT as iedb_alignment with:
    input:
        peps_fa = rules.pep_to_fasta.output.peps_fa
    output:
        alignment = temp(os.path.join("results", "neo_score", "mhc1peps_len{length}_split{mhc1_split}.iedb.csv"))

rulename = "human_alignment"
use rule human_alignment from neoscoreT as human_alignment with:
    input:
        peps_fa = rules.pep_to_fasta.output.peps_fa
    output:
        alignment = temp(os.path.join("results", "neo_score", "mhc1peps_len{length}_split{mhc1_split}.human.csv"))

rulename = "calculate_scores"
use rule calculate_scores from neoscoreT as calculate_scores with:
    input:
        peps = os.path.join("results", "netmhcpan", "mhc1peps_len{length}_split{mhc1_split}.txt"),
        human_alignment = rules.human_alignment.output.alignment,
        iedb_alignment = rules.iedb_alignment.output.alignment
    output:
        pep_scores = temp(os.path.join("results", "neo_score", "mhc1peps_len{length}_split{mhc1_split}.scores.csv"))

rulename = "merge_score_files"
use rule merge_score_files from neoscoreT as merge_score_files with:
    input:
        pep_scores = expand(os.path.join("results", "neo_score", "mhc1peps_len{length}_split{mhc1_split}.scores.csv"), length=LENGTHS, mhc1_split=SPLITS)
    output:
        merged_pep_scores = os.path.join("results", "neo_score", "neo_scores.csv")
