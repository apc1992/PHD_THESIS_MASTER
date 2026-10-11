#!/usr/bin/env python3
import sys
import pandas as pd
from pathlib import Path

raw_file = Path("results/netmhc2pan_raw_output.txt")
out_file = Path("results/netmhc2pan_parsed_results.csv")

if not raw_file.exists():
    sys.exit(f"Error: {raw_file} not found.")

records = []

with open(raw_file, "r") as f:
    for line in f:
        line = line.strip()
        
        # Skip comments, empty lines, and table separators
        if not line or line.startswith("#") or line.startswith("-"):
            continue
            
        parts = line.split()
        
        # Header or invalid line check
        if parts[0] == "Pos" or len(parts) < 13:
            continue
            
        # Parse standard NetMHCIIpan EL+BA output line
        # Handle potential line formatting where Pos is element 0 and MHC is element 1
        try:
            # Check if line starts with integer position
            if not parts[0].isdigit():
                continue
                
            pos = parts[0]
            mhc = parts[1]
            peptide = parts[2]
            core = parts[4] if len(parts) > 4 else ""
            score_el = parts[8] if len(parts) > 8 else ""
            rank_el = parts[9] if len(parts) > 9 else ""
            score_ba = parts[11] if len(parts) > 11 else ""
            rank_ba = parts[12] if len(parts) > 12 else ""
            affinity = parts[13] if len(parts) > 13 else ""
            bind_level = parts[14] if len(parts) > 14 else ""
            
            records.append({
                "Pos": pos,
                "MHC": mhc,
                "Peptide": peptide,
                "Core": core,
                "Score_EL": score_el,
                "Rank_EL": rank_el,
                "Score_BA": score_ba,
                "Rank_BA": rank_ba,
                "Affinity_nM": affinity,
                "BindLevel": bind_level
            })
        except ValueError:
            continue

df = pd.DataFrame(records)
df.to_csv(out_file, index=False)
print(f"Successfully processed {len(df)} predictions into {out_file}")
