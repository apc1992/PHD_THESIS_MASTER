import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple

EPSILON = 1e-5

def validate_and_extract_cytokine_cols(df: pd.DataFrame) -> List[str]:
    """Validates that activity columns follow the CDXX naming convention."""
    cyto_cols = [c for c in df.columns if c.startswith(('CD4', 'CD8')) and not c.endswith('_score')]
    if not cyto_cols:
        raise ValueError("Critical Error: No cytokine columns starting with 'CD4' or 'CD8' were found in input CSV.")
    logging.info(f"Identified {len(cyto_cols)} valid cytokine activity columns: {cyto_cols}")
    return cyto_cols

def compute_wt_relative_ratios(
    cytokines_df: pd.DataFrame, 
    wt_seq: str, 
    cyto_cols: List[str]
) -> Dict[Tuple[int, str, str, str], float]:
    """Computes donor- and position-specific relative activity ratios against WT baselines."""
    ratios = {}
    donors = cytokines_df['donor'].unique()

    for donor in donors:
        donor_data = cytokines_df[cytokines_df['donor'] == donor]
        
        for pos in range(1, len(wt_seq) + 1):
            wt_aa = wt_seq[pos - 1]
            wt_row = donor_data[(donor_data['position'] == pos) & (donor_data['aminoacid'] == wt_aa)]
            
            for col in cyto_cols:
                wt_val = wt_row[col].values[0] if not wt_row.empty else EPSILON
                wt_val = wt_val if wt_val > 0 else EPSILON
                
                pos_data = donor_data[donor_data['position'] == pos]
                for _, row in pos_data.iterrows():
                    aa = row['aminoacid']
                    act_val = row[col]
                    if pd.notna(act_val):
                        ratios[(donor, col, pos, aa)] = float(act_val) / wt_val

    return ratios

def score_mimotopes(
    variants_df: pd.DataFrame, 
    wt_seq: str, 
    ratios: Dict[Tuple[int, str, str, str], float], 
    cyto_cols: List[str]
) -> pd.DataFrame:
    """Calculates multiplicative functional scores across mutated positions for all mimotopes."""
    scored_records = []

    for _, row in variants_df.iterrows():
        seq = row['mimotope']
        mutations = [(i + 1, aa) for i, (aa, wt) in enumerate(zip(seq, wt_seq)) if aa != wt]
        
        rec = {'mimotope': seq}
        donors = set(key[0] for key in ratios.keys())
        
        for col in cyto_cols:
            donor_scores = []
            for donor in donors:
                donor_ratio = 1.0
                for pos, aa in mutations:
                    donor_ratio *= ratios.get((donor, col, pos, aa), 1.0)
                donor_scores.append(donor_ratio)
            
            rec[f"{col}_score"] = np.mean(donor_scores) if donor_scores else np.nan
            
        scored_records.append(rec)

    return pd.DataFrame(scored_records)
