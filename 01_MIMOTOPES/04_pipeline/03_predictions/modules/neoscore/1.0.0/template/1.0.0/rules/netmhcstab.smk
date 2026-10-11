import re
import os
import pandas as pd
import numpy as np
from datetime import datetime
#from scripts.utils.path_builder import *
#writer = Writer(config["benchmark"], directory = "results", pipeline = "neoscore", module = "netmhcstab")
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
rule run_netmhcstab:
    message: "Running NetMHCstabPan for {wildcards.mhc1_allele}, {wildcards.length} and {wildcards.mhc1_split}"
    version: "4.1.0"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/netmhcstab/run_netmhcstab/{mhc1_allele}_{length}_{mhc1_split}.log",
        stderr = "logs/neoscore/netmhcstab/run_netmhcstab/{mhc1_allele}_{length}_{mhc1_split}.err"
    benchmark:
        "benchmarks/neoscore/netmhcstab/run_netmhcstab/{mhc1_allele}_{length}_{mhc1_split}.benchmark.tsv"
    params:
        alleles = lambda wildcards: local_utils.get_mhc1_alleles(wildcards.mhc1_allele, config),
        container = config["containers"]["netmhcstab"]
    wrapper:
       "file:tools/netmhcstab/1.0.1"

rulename = "process_netmhcstabpan_output"
rule process_netmhcstabpan_output:
    message: "Processing netmhcstab output for {wildcards.mhc1_allele}, {wildcards.length} and {wildcards.mhc1_split}"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/netmhcstab/process_netmhcstabpan_output/{mhc1_allele}_{length}_{mhc1_split}.log",
        stderr = "logs/neoscore/netmhcstab/process_netmhcstabpan_output/{mhc1_allele}_{length}_{mhc1_split}.err"
    benchmark:
        "benchmarks/neoscore/netmhcstab/process_netmhcstabpan_output/{mhc1_allele}_{length}_{mhc1_split}.benchmark.tsv"
    run:
        # Initialize output content with header line
        header = "Pos,MHC,Peptide,Identity,Pred,Thalf(h),%Rank_Stab,BindLevel\n"
        lines = [header]

        # Process each line in the raw output file
        with open(input.raw_netmhc_out) as f:
            for line in f:
                if "PEPLIST " in line:
                    # Remove "<=" and unnecessary whitespace, replace with commas
                    line = re.sub(r'<=', '', line)
                    line = re.sub(r'[ \t]+', ',', line)
                    line = line[1:]
                    lines.append(line)

        # Write the processed output to the specified CSV file
        with open(output.netmhcstab_pred, "w") as out_f:
            out_f.writelines(lines)

rulename = "merge_netmhcstabpan_predictions"
rule merge_netmhcstabpan_predictions:
    message: "Merging NetMHCPan predictions"
    threads: lambda wildcards, rulename=rulename: config["cluster"][rulename]["n"]
    resources:
        mem = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"],
        mem_mb = lambda wildcards, rulename=rulename: config["cluster"][rulename]["mem"] * 1024,
        time_min = lambda wildcards, rulename=rulename: config["cluster"][rulename]["time_min"]
    log:
        stdout = "logs/neoscore/netmhcpan/merge_netmhcstabpan_predictions/run.log",
        stderr = "logs/neoscore/netmhcpan/merge_netmhcstabpan_predictions/run.err"
    benchmark:
        "benchmarks/neoscore/netmhcpan/merge_netmhcstabpan_predictions/benchmark.tsv"
    run:
        # Merge MHC-I predictions
        merged_df = local_utils.merge_files(input.peps, output.merged)

        # Parse merged predictions
        merged_df = merged_df[["MHC", "Peptide", "Pred", "Thalf(h)", "%Rank_Stab", "BindLevel"]]
        merged_df = merged_df.set_index(['Peptide', 'MHC']).unstack('MHC')
        merged_df.columns = [f'{col[1]}_{col[0]}' for col in merged_df.columns]
        merged_df = merged_df.reset_index()

        # Write all the columns from the MHC predictions into a table
        merged_df.to_csv(output.full_predictions, index=False)

        # Identify Stab columns
        stab_columns = [col for col in merged_df.columns if '%Rank_Stab' in col]

        # Find binders for each row to find EL and BA binders
        merged_df['Stab_Binders'] = merged_df.apply(local_utils.find_binders, columns=stab_columns, axis=1)

        # Merge predictions into the main table
        main_df = pd.read_csv(input.main_table, sep="\t")
        combined_table = pd.merge(main_df, merged_df[['Peptide', 'Stab_Binders']], on='Peptide', how='left')
        combined_table.to_csv(output.combined_table, index=False)
