#!/usr/bin/env python3
import os
import yaml
import pandas as pd

def main():
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    wt_seq = config["wt_protein_seq"]
    wt_pep = config["wt_peptide"]
    lengths = config["kmer_lengths"]
    candidates_file = config["input_candidates"]

    # Read candidate mimotopes and include WT control
    candidates_df = pd.read_csv(candidates_file)
    candidates = candidates_df["peptide"].dropna().unique().tolist()
    if wt_pep not in candidates:
        candidates.append(wt_pep)

    # Locate WT peptide in BIRC5 sequence (1-based coordinates)
    pep_len = len(wt_pep)
    p_start_wt = wt_seq.find(wt_pep)
    if p_start_wt == -1:
        raise ValueError(f"WT peptide {wt_pep} not found in BIRC5 sequence.")
    
    p_start_wt_1based = p_start_wt + 1
    p_end_wt_1based = p_start_wt + pep_len

    kmers_list = []

    for mimotope in candidates:
        # Construct mutant BIRC5 sequence
        mut_seq = wt_seq[:p_start_wt] + mimotope + wt_seq[p_end_wt_1based:]

        for k in lengths:
            # Sliding window start positions ensuring mimotope is fully enclosed
            min_start = max(1, p_end_wt_1based - k + 1)
            max_start = min(p_start_wt_1based, len(wt_seq) - k + 1)

            for s in range(min_start, max_start + 1):
                e = s + k - 1
                long_pep = mut_seq[s - 1 : e]
                kmers_list.append({
                    "mimotope": mimotope,
                    "long_peptide": long_pep,
                    "length_peptide": k,
                    "p.start": s,
                    "p.end": e
                })

    kmers_df = pd.DataFrame(kmers_list)
    os.makedirs("output", exist_ok=True)
    
    # Save coordinate map
    kmers_df.to_csv("output/mimotope_flanking_kmers.csv", index=False)

    # Save unique long peptides for NetMHCIIpan input
    unique_peps = kmers_df["long_peptide"].unique()
    with open("output/netmhc2pan_input_peptides.txt", "w") as f:
        for pep in unique_peps:
            f.write(f"{pep}\n")

    print(f"Generated {len(kmers_df)} total k-mers across {len(candidates)} candidates.")
    print(f"Saved {len(unique_peps)} unique peptides for NetMHCIIpan.")

if __name__ == "__main__":
    main()
