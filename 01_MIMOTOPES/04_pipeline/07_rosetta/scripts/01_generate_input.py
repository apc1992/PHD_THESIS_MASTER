# generate_inputs.py
import os
import glob
import pandas as pd

# 1. Load the selected candidates list
candidates_csv = "../06_alphafold3/output/rosetta_candidates_selected.csv"
if not os.path.exists(candidates_csv):
    raise FileNotFoundError(f"Missing candidates list: {candidates_csv}")

candidates_df = pd.read_csv(candidates_csv)
selected_peptides = set(candidates_df["peptide"].tolist())

# 2. Find all matching CIF files in your results directory
# Matches: results/HLA_A0201_B2M_Candidate_<ID>_<PEPTIDE>/.../*.cif
cif_pattern = "../06_alphafold3/results/HLA_A0201_B2M_Candidate_*/*/*/*/*/*.cif"
cif_files = glob.glob(cif_pattern)

records = []
for path in cif_files:
    filename = os.path.basename(path)
    # Check if this CIF corresponds to one of our selected peptides
    for pep in selected_peptides:
        if f"_{pep}/" in path or f"_{pep}_" in filename:
            records.append({
                "peptide": pep,
                "cif": os.path.abspath(path)
            })
            break

# 3. Write out to the designated csv format
input_df = pd.DataFrame(records).drop_duplicates(subset=["peptide"])
input_df.to_csv("input_libraries.csv", index=False)

print(f"[+] Successfully mapped {len(input_df)}/{len(selected_peptides)} candidates to CIF paths.")
print("[+] Saved file to: input_libraries.csv")
