import os
import yaml
import sys
import pandas as pd
import numpy as np
import math

##############################################################################
#################### File System Operations and Management ###################
##############################################################################

def merge_files(input_files, output_file):
    """
    Merge multiple CSV or TSV files into a single file, maintaining the same format for the output.

    Args:
        input_files (list): A list of file paths for the input CSV or TSV files.
        output_file (str): The file path for the output merged file.

    Returns:
        pd.DataFrame: Merged DataFrame of all input files.
    """
    dfs = []
    for file in input_files:
        # Detect if the file is CSV or TSV based on the file extension
        if file.endswith('.csv'):
            df = pd.read_csv(file)
        elif file.endswith('.tsv'):
            df = pd.read_csv(file, delimiter='\t')
        else:
            raise ValueError(f"Unsupported file format: {file}. Only CSV and TSV are supported.")
        
        dfs.append(df)

    # Concatenate all DataFrames
    merged_df = pd.concat(dfs, ignore_index=True)

    # Determine output format based on the output file extension
    if output_file.endswith('.csv'):
        merged_df.to_csv(output_file, index=False)
    elif output_file.endswith('.tsv'):
        merged_df.to_csv(output_file, sep='\t', index=False)
    else:
        raise ValueError(f"Unsupported output file format: {output_file}. Only CSV and TSV are supported.")
    
    return merged_df

##############################################################################
############################# Pipeline Utilities #############################
##############################################################################

def set_num_splits(input_file, peptides_per_chunk=2500):
    """
    Sets the number of split chunks based on a target number of 
    peptides per file.
    """
    df = pd.read_csv(input_file, sep="\t")
    
    # Identify peptide column
    peptide_col = next((col for col in df.columns if 'peptide' in col.lower() or 'pep' in col.lower()), None)
    
    if peptide_col is None:
        raise ValueError("No peptide column found.")
    
    # Filter for valid MHC-I lengths (8-12) and get unique count
    peptides = df[peptide_col].astype(str).str.strip()
    unique_count = peptides[peptides.str.len().between(8, 12)].nunique()

    # Calculate splits: Total Unique / Target Peptides per Chunk
    # Using 2500-5000 is usually ideal for parallel bioinformatics pipelines
    splits = math.ceil(unique_count / peptides_per_chunk)

    return max(1, splits)

def get_pep_lengths(input_file):
    """
    Search for a column named 'pep' or 'peptide' case-insensitive, and extract the different peptide lengths.

    Args:
        input_file (str): File path to the input file.

    Returns:
        List of peptide lengths (list->str).
    """
    # Load data
    df = pd.read_csv(input_file, sep="\t")
    
    # Find the peptide column
    peptide_col = next((col for col in df.columns if 'peptide' in col.lower() or 'pep' in col.lower()), None)
    
    if peptide_col is None:
        raise ValueError("No peptide column found.")

    # 1. Calculate lengths
    lengths = df[peptide_col].astype(str).str.len()

    # 2. Filter for lengths between 8 and 12 (inclusive)
    filtered_lengths = lengths[lengths.between(8, 12)]

    # 3. Get unique, sorted values
    unique_lengths = sorted(filtered_lengths.unique())
    
    # Return as a list of strings
    return [str(length) for length in unique_lengths]

def extract_and_split_peptides(input_file, n_files, lengths, prefix):
    """
    Search for a column named 'pep' or 'peptide' case-insensitive, group peptides by length, and split the unique peptides across N files.

    Args:
        input_file (str): File path to the input file.
        n_files (int): Number of splits.
        prefix (str): File prefix for each split file.

    Returns:
        None.
    """
    # Load data and identify the peptide column
    n_files = int(n_files)
    lengths = [int(l) for l in lengths]
    lengths = [l for l in lengths if 8 <= l <= 12]  # keep only valid MHC-I lengths
    df = pd.read_csv(input_file, sep="\t")
    
    # Find the peptide column (case-insensitive search for 'peptide' or 'pep')
    peptide_col = next((col for col in df.columns if 'peptide' in col.lower() or 'pep' in col.lower()), None)
    if peptide_col is None:
        raise ValueError("No peptide column found. Ensure a column contains 'peptide' or 'pep'.")
    
    # Filter peptides by length and remove duplicates
    peptides = df[peptide_col].astype(str).str.strip()  # Convert to string and trim whitespace
    peptides = peptides[(peptides.str.len() >= 8) & (peptides.str.len() <= 12)].unique()  # Length filter and deduplication
    
    for length in lengths:
        # Filter for specific length
        peps_at_len = [p for p in peptides if len(p) == length]

        # Split peptides into nearly equal parts using numpy.array_split
        peptide_chunks = np.array_split(peps_at_len, n_files)

        # Write each chunk to a separate file
        for i, chunk in enumerate(peptide_chunks):
            output_path = f"{prefix}_len{length}_split{i}.txt"
            with open(output_path, 'w') as f:
                f.write('\n'.join(chunk))

def get_mhc1_alleles(hla_gene, config):
    """
    Get all the alleles in the configuration given an HLA-I gene name.

    Args:
        input_file (str): HLA-I gene name, i.e. "A".

    Returns:
        List of HLA-I alleles in the configuration for that gene.
    """
    try:
        return config["netmhcpan"]["alleles"][hla_gene]
    except KeyError:
        return f"HLA-I gene name key '{hla_gene}' not found in configuration."

def check_mhc1_alleles(config):
    """
    Check all the alleles in the configuration against the NetMHCPan whitelist.
    """

    # Extract only non-empty alleles from config
    alleles_in_config = [
        allele.strip()
        for key in ["A", "B", "C"]
        for allele in config["netmhcpan"]["alleles"].get(key, "").split(",")
        if allele.strip()
    ]

    # Load whitelist alleles
    with open(config["netmhcpan"]["available_alleles"], "r") as infile:
        alleles_in_whitelist = {line.strip() for line in infile}

    # Find missing alleles
    missing_alleles = [
        allele
        for allele in alleles_in_config
        if allele not in alleles_in_whitelist
    ]

    if missing_alleles:
        print(
            f"Error: The following alleles are not in the whitelist: "
            f"{', '.join(missing_alleles)}"
        )
        sys.exit(1)
    else:
        print("All alleles in config are available!")

def find_binders(row, columns, threshold=2):
    """
    Identify binding alleles for a given row based on prediction scores.

    Args:
        row (pd.Series): A single row from the DataFrame containing MHC prediction scores.
        columns (List[str]): A list of column names to check, typically containing scores for 
                             each MHC allele's binding prediction.
        threshold (float, optional): The threshold value for binding prediction. Alleles with 
                                     scores below this value are considered binders. Default is 2.

    Returns:
        str: A comma-separated string of allele names that are considered binders.
    """
    binders = [col.split('_')[0] for col in columns if row[col] < threshold]
    return ','.join(binders)

def generate_splits(num_splits):
    """
    Generate a list of zero-padded split identifiers based on the number of splits defined in the configuration.

    Returns:
        list: A list of zero-padded split identifiers as strings.
    """
    splits = ["{}".format(str(i)) for i in range(num_splits)]
    return splits
