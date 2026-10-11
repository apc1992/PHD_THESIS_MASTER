#!/usr/bin/env python3
import re
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# Paths
HLA_FILE = Path("leukopaks_hla_typing.csv")
PRED_FILE = Path("results/netmhc2pan_parsed_results.csv")
MIMO_FILE = Path("output/mimotope_flanking_kmers.csv")
FIG_DIR = Path("figures")
FIG_DIR.mkdir(exist_ok=True)

# 1. Load Data
if not HLA_FILE.exists() or not PRED_FILE.exists() or not MIMO_FILE.exists():
    sys.exit("Error: One or more required input files are missing.")

hla_df = pd.read_csv(HLA_FILE, index_col=0)
pred_df = pd.read_csv(PRED_FILE)
mimo_df = pd.read_csv(MIMO_FILE)

# Ensure Rank_EL is numeric
pred_df["Rank_EL"] = pd.to_numeric(pred_df["Rank_EL"], errors="coerce")
pred_df = pred_df.dropna(subset=["Rank_EL"])

# Map long peptides back to core mimotopes
seq_to_mimo = dict(zip(mimo_df["long_peptide"].astype(str).str.strip(), mimo_df["mimotope"].astype(str).str.strip()))
pred_df["Peptide"] = pred_df["Peptide"].astype(str).str.strip()
pred_df["Mimotope"] = pred_df["Peptide"].map(seq_to_mimo)
pred_df = pred_df.dropna(subset=["Mimotope"])

valid_netmhc_alleles = set(pred_df["MHC"].unique())

def get_4digit_code(val):
    if pd.isna(val):
        return None
    val = str(val).strip()
    if "*" in val:
        val = val.split("*")[-1]
    parts = val.split(":")
    if len(parts) >= 2:
        return f"{parts[0].zfill(2)}{parts[1].zfill(2)}"
    clean = re.sub(r"[^\d]", "", val)
    return clean[:4] if len(clean) >= 4 else None

def format_drb(val):
    code = get_4digit_code(val)
    return f"DRB1_{code}" if code else None

def format_dp(a, b):
    code_a, code_b = get_4digit_code(a), get_4digit_code(b)
    return f"HLA-DPA1{code_a}-DPB1{code_b}" if (code_a and code_b) else None

def format_dq(a, b):
    code_a, code_b = get_4digit_code(a), get_4digit_code(b)
    return f"HLA-DQA1{code_a}-DQB1{code_b}" if (code_a and code_b) else None

def get_matching_rows(df, prefix):
    return [r for r in df.index if str(r).strip().upper().startswith(prefix.upper())]

drb_rows = get_matching_rows(hla_df, "DRB")
dpa_rows = get_matching_rows(hla_df, "DPA")
dpb_rows = get_matching_rows(hla_df, "DPB")
dqa_rows = get_matching_rows(hla_df, "DQA")
dqb_rows = get_matching_rows(hla_df, "DQB")

# Target donors requested
target_donors = ["LEUKOPAK1-AC", "LEUKOPAK2-AC", "LEUKOPAK4-AC", "LEUKOPAK9-SC", "LEUKOPAK11-SC"]
available_donors = [d for d in target_donors if d in hla_df.columns]

donor_alleles = {}
for d in available_donors:
    alleles = set()
    for r in drb_rows:
        formatted = format_drb(hla_df.loc[r, d])
        if formatted and formatted in valid_netmhc_alleles:
            alleles.add(formatted)
    for r_a in dpa_rows:
        for r_b in dpb_rows:
            formatted = format_dp(hla_df.loc[r_a, d], hla_df.loc[r_b, d])
            if formatted and formatted in valid_netmhc_alleles:
                alleles.add(formatted)
    for r_a in dqa_rows:
        for r_b in dqb_rows:
            formatted = format_dq(hla_df.loc[r_a, d], hla_df.loc[r_b, d])
            if formatted and formatted in valid_netmhc_alleles:
                alleles.add(formatted)
    donor_alleles[d] = alleles

# Extract best hits preserving Long Peptide and Best MHC
donor_results = []
for d, alleles in donor_alleles.items():
    if not alleles:
        continue
    sub = pred_df[pred_df["MHC"].isin(alleles)]
    if sub.empty:
        continue
    
    idx_min = sub.groupby("Mimotope")["Rank_EL"].idxmin()
    best_df = sub.loc[idx_min, ["Mimotope", "MHC", "Peptide", "Rank_EL"]].copy()
    best_df.rename(columns={"MHC": "Best_MHC", "Peptide": "Best_Long_Peptide", "Rank_EL": "Min_Rank_EL"}, inplace=True)
    best_df["Donor"] = d
    donor_results.append(best_df)

if not donor_results:
    sys.exit("Error: No matching donor alleles found in NetMHCIIpan output.")

summary_df = pd.concat(donor_results, ignore_index=True)
summary_df.to_csv("results/donor_mimotope_detailed_best_binders.csv", index=False)

# Pivot for Filtered Heatmap
pivot_df = summary_df.pivot(index="Donor", columns="Mimotope", values="Min_Rank_EL")

plt.figure(figsize=(12, 4.5))
sns.heatmap(
    pivot_df,
    annot=True,
    fmt=".2f",
    cmap="YlOrRd_r",
    vmax=10.0,
    vmin=0.5,
    cbar_kws={'label': 'Best %Rank_EL (Lower = Stronger Binding)'}
)
plt.title("Class II Best Binding Rank (%Rank_EL) — Selected Donors (1, 2, 4, 9, 11)")
plt.xlabel("Mimotope")
plt.ylabel("Donor")
plt.tight_layout()
plt.savefig(FIG_DIR / "donor_mimotope_rank_heatmap_filtered.png", dpi=300)
plt.close()

print(f"Detailed binder table saved to results/donor_mimotope_detailed_best_binders.csv")
print(f"Filtered heatmap saved to figures/donor_mimotope_rank_heatmap_filtered.png")
