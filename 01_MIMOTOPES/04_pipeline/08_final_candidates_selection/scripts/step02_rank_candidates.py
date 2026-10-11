#!/usr/bin/env python3
"""step02_rank_candidates.py

Calculates Z-scores strictly on Step 1 passers using the original formula:
Composite Score = 0.30 * Z_EL_Rank + 0.30 * Z_AF3_PAE + 0.25 * Z_Rosetta_dG - 0.15 * Z_T1/2
"""

import logging
import numpy as np
import pandas as pd

logger = logging.getLogger("08_final_candidates")


def calculate_zscores_and_rank(
    df: pd.DataFrame, clip_min: float = -3.0, clip_max: float = 3.0
) -> tuple[pd.DataFrame, pd.DataFrame]:
  df_master = df.copy()

  # Isolate all candidates that passed Step 1 gates
  df_passed = df_master[df_master["passed_all_gates"] == True].copy()

  if df_passed.empty:
    logger.warning(
        "No candidates passed Step 1 gates. Falling back to entire master dataset."
    )
    df_passed = df_master.copy()

  def get_zscore(series: pd.Series) -> pd.Series:
    std = series.std()
    if std == 0 or np.isnan(std):
      return pd.Series(0.0, index=series.index)
    return ((series - series.mean()) / std).clip(clip_min, clip_max)

  # Z-score normalization computed STRICTLY across the Step 1 passed subset
  df_passed["z_el_rank"] = get_zscore(df_passed["elution_rank"])
  df_passed["z_af3_pae"] = get_zscore(df_passed["af3_mean_pae"])
  df_passed["z_rosetta_dg"] = get_zscore(df_passed["rosetta_dG"])
  df_passed["z_t12"] = get_zscore(df_passed["half_life"])

  # Original Composite Score formula:
  # Composite = 0.30*Z_EL_Rank + 0.30*Z_AF3_PAE + 0.25*Z_Rosetta_dG - 0.15*Z_T1/2
  df_passed["composite_score"] = (
      0.30 * df_passed["z_el_rank"]
      + 0.30 * df_passed["z_af3_pae"]
      + 0.25 * df_passed["z_rosetta_dg"]
      - 0.15 * df_passed["z_t12"]
  )

  df_ranked_validation = df_passed.sort_values("composite_score").reset_index(
      drop=True
  )

  logger.info(
      f"Ranked {len(df_ranked_validation)} Step 1 passed leads with zero-centered Z-scores."
  )
  return df_master, df_ranked_validation
