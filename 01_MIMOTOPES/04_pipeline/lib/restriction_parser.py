import re
import pandas as pd
from pathlib import Path
from typing import Dict, List

STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def parse_restrictions(restrictions_path: Path, wt_seq: str) -> Dict[int, List[str]]:
    """
    Parses restrictions.csv for a target peptide.
    Returns a dictionary mapping ALL 0-indexed positions (0 to len(wt_seq)-1) 
    to their allowed amino acids. Unconstrained positions default strictly to [WT].
    """
    if not restrictions_path.exists():
        raise FileNotFoundError(f"Restrictions file not found at: {restrictions_path}")

    # Initialize all positions with only the Wild-Type amino acid
    parsed_constraints: Dict[int, List[str]] = {
        i: [wt_seq[i]] for i in range(len(wt_seq))
    }

    df = pd.read_csv(restrictions_path)
    if "peptide" in df.columns:
        df = df[df["peptide"] == wt_seq]

    for _, row in df.iterrows():
        raw_aa = str(row["aminoacids"]).strip()
        
        if "position" in row and pd.notna(row["position"]):
            pos_val = int(row["position"])
        else:
            match = re.search(r'\.(\d+)$', raw_aa)
            if match:
                pos_val = int(match.group(1))
                raw_aa = raw_aa[:match.start()]
            else:
                raise ValueError(f"Could not parse position from row: {row.to_dict()}")

        clean_aa_set = set(re.findall(r'[A-Z]', raw_aa.upper())) & STANDARD_AMINO_ACIDS
        zero_idx = pos_val - 1
        
        # Add WT amino acid to the allowed set
        clean_aa_set.add(wt_seq[zero_idx])
        parsed_constraints[zero_idx] = sorted(list(clean_aa_set))

    return parsed_constraints
