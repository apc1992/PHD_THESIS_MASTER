import logging
import numpy as np
import pandas as pd
from typing import Tuple, Dict

logger = logging.getLogger(__name__)

# Literature Cutoffs (%Rank EL)
SB_EL_CUTOFF = 0.5         # Strong Binder %Rank EL
WB_EL_CUTOFF = 2.0         # Weak Binder %Rank EL
BORDERLINE_EL_CUTOFF = 5.0  # Borderline Binder %Rank EL

def compute_composite_scores(df: pd.DataFrame) -> pd.DataFrame:
    """Computes geometric mean composite scores directly from summary cytokine metrics."""
    df = df.copy()
    
    cd8_bb = np.clip(df.get('CD8_41BB_score', 0.0).fillna(0.0), 0, None)
    cd8_ifn = np.clip(df.get('CD8_IFNg_score', 0.0).fillna(0.0), 0, None)
    cd4_bb = np.clip(df.get('CD4_41BB_score', 0.0).fillna(0.0), 0, None)
    cd4_ifn = np.clip(df.get('CD4_IFNg_score', 0.0).fillna(0.0), 0, None)
    
    df['CD8_Composite'] = np.sqrt(cd8_bb * cd8_ifn)
    df['CD4_Composite'] = np.sqrt(cd4_bb * cd4_ifn)
    
    return df

def assign_unified_tier(row: pd.Series) -> str:
    """Assigns structural tiers driven primarily by %Rank EL."""
    el = row.get('HLA-A*02:01_%Rank_EL', row.get('EL_rank_A0201', np.nan))

    if pd.isna(el):
        return 'Tier 4: Non-Binder'

    if el <= SB_EL_CUTOFF:
        return 'Tier 1: Strong Binder'
    elif el <= WB_EL_CUTOFF:
        return 'Tier 2: Weak Binder'
    elif el <= BORDERLINE_EL_CUTOFF:
        return 'Tier 3: Borderline Binder'
    else:
        return 'Tier 4: Non-Binder'

def select_candidates_struct_first(
    df: pd.DataFrame, 
    target_per_cat: int = 50,
    pos_pair_cap: int = 2
) -> Tuple[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Selects top 50 candidates per variant category, completely excluding WT."""
    
    # Exclude reference WT rows from selection processing
    df_var = df[df['category'].str.upper() != 'WT'].copy()
    
    el_col = 'HLA-A*02:01_%Rank_EL' if 'HLA-A*02:01_%Rank_EL' in df_var.columns else 'EL_rank_A0201'
    aff_col = 'HLA-A*02:01_Aff(nM)' if 'HLA-A*02:01_Aff(nM)' in df_var.columns else 'BA_nM_A0201'

    # 1. 1-AA Arm: Priority sort by %Rank EL (asc) -> Aff (asc) -> CD8 Composite (desc)
    df_1aa = df_var[df_var['category'].str.upper().isin(['1AA', '1-AA'])].copy()
    selected_1aa = df_1aa.sort_values(
        by=[el_col, aff_col, 'CD8_Composite'], 
        ascending=[True, True, False]
    ).head(target_per_cat)

    # 2. 2-AA Arm: Structural priority with positional diversity cap
    df_2aa = df_var[df_var['category'].str.upper().isin(['2AA', '2-AA'])].copy()
    df_2aa_sorted = df_2aa.sort_values(
        by=[el_col, aff_col, 'CD8_Composite'], 
        ascending=[True, True, False]
    )
    if 'position' in df_2aa_sorted.columns:
        selected_2aa = df_2aa_sorted.groupby('position').head(pos_pair_cap).head(target_per_cat)
    else:
        selected_2aa = df_2aa_sorted.head(target_per_cat)

    # 3. PS-CPL Arm: Structural priority (%Rank EL -> Aff -> CD8 Composite)
    df_pscpl = df_var[df_var['category'].str.upper().isin(['PS-CPL', 'PSCPL'])].copy()
    selected_pscpl = df_pscpl.sort_values(
        by=[el_col, aff_col, 'CD8_Composite'], 
        ascending=[True, True, False]
    ).head(target_per_cat)

    # Combine selected cohorts and deduplicate using variant sequence identifier ('mimotope')
    combined_selected = pd.concat([selected_1aa, selected_2aa, selected_pscpl]).drop_duplicates(subset=['mimotope'])
    
    category_dict = {
        '1AA': selected_1aa,
        '2AA': selected_2aa,
        'PS-CPL': selected_pscpl
    }
    
    return combined_selected, category_dict
