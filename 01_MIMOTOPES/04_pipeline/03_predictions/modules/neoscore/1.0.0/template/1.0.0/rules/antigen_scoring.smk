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
rule pep_to_fasta:
    message: "Converting peptides to fasta for {wildcards.length} and {wildcards.mhc1_split}"
    version: "1.0.0"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/antigen_scoring/pep_to_fasta/{length}_{mhc1_split}.log",
        stderr = "logs/neoscore/antigen_scoring/pep_to_fasta/{length}_{mhc1_split}.err"
    benchmark:
        "benchmarks/neoscore/antigen_scoring/pep_to_fasta/{length}_{mhc1_split}.benchmark.tsv"
    run:
        with open(input.peps, "r") as infile, open(output.peps_fa, "w") as outfile:
            for i, line in enumerate(infile):
                pep = line.strip()
                if pep:
                    # Use the line index as the header
                    outfile.write(f">{i}\n{pep}\n")

rulename = "iedb_alignment"
rule iedb_alignment:
    message: "Aligning to IEDB -- {wildcards.length} and {wildcards.mhc1_split}"
    version: "1.0.0"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/antigen_scoring/iedb_alignment/{length}_{mhc1_split}.log",
        stderr = "logs/neoscore/antigen_scoring/iedb_alignment/{length}_{mhc1_split}.err"
    params:
        database = config["blast_databases"]["iedb"],
        args = "-evalue 100000000 -matrix BLOSUM62 -gapopen 11 -gapextend 1",
        task = "-task blastp-short",
        output_format = "'10 qseqid sseqid qseq sseq'",
        container = config["containers"]["blast"]
    benchmark:
        "benchmarks/neoscore/antigen_scoring/iedb_alignment/{length}_{mhc1_split}.benchmark.tsv"
    wrapper:
       "file:tools/blast/2.15.0"

rulename = "human_alignment"
rule human_alignment:
    message: "Aligning to HUMAN -- {wildcards.length} and {wildcards.mhc1_split}"
    version: "1.0.0"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/antigen_scoring/human_alignment/{length}_{mhc1_split}.log",
        stderr = "logs/neoscore/antigen_scoring/human_alignment/{length}_{mhc1_split}.err"
    params:
        database = config["blast_databases"]["human"],
        args = "-evalue 100000000 -matrix BLOSUM62 -gapopen 11 -gapextend 1",
        task = "-task blastp-short",
        output_format = "'10 qseqid sseqid qseq sseq'",
        container = config["containers"]["blast"]
    benchmark:
        "benchmarks/neoscore/antigen_scoring/human_alignment/{length}_{mhc1_split}.benchmark.tsv"
    wrapper:
       "file:tools/blast/2.15.0"

rulename = "calculate_scores"
rule calculate_scores:
    message: "Calculating Neoantigen scores -- {wildcards.length} and {wildcards.mhc1_split}"
    version: "1.0.0"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/antigen_scoring/calculate_scores/{length}_{mhc1_split}.log",
        stderr = "logs/neoscore/antigen_scoring/calculate_scores/{length}_{mhc1_split}.err"
    benchmark:
        "benchmarks/neoscore/antigen_scoring/calculate_scores/{length}_{mhc1_split}.benchmark.tsv"
    params:
        container = config["containers"]["neo_score"]
    wrapper:
       "file:tools/neo_scores/1.0.0"

rulename = "merge_score_files"
rule merge_score_files:
    message: "Merging score files"
    version: "1.0.0"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/antigen_scoring/merge_score_files/run.log",
        stderr = "logs/neoscore/antigen_scoring/merge_score_files/run.err"
    benchmark:
        "benchmarks/neoscore/antigen_scoring/merge_score_files/benchmark.tsv"
    run:
        merged_scores = local_utils.merge_files(input.pep_scores, output.merged_pep_scores)
