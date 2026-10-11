#!/usr/bin/env python3
"""step01_compile_master_table.py

Compiles outputs from Module 05, AF3, and Rosetta into a master candidate table.
Filters strictly for valid mimotope categories (1-AA, 2-AA, PS-CPL) and AF3 passers.
"""

import logging
from pathlib import Path
import numpy as np
import pandas as pd

logger = logging.getLogger("08_final_candidates")

VALID_CATEGORIES = {"1-AA", "2-AA", "PS-CPL", "1AA", "2AA", "PSCPL"}


def standardize_peptide_column(df: pd.DataFrame) -> pd.DataFrame:
  """Ensures candidate variant sequence column is named 'peptide'."""
  df = df.copy()

  if "mimotope" in df.columns:
    if "Peptide" in df.columns:
      df = df.rename(columns={"Peptide": "wt_peptide"})
    if "peptide" in df.columns:
      df = df.drop(columns=["peptide"])
    df = df.rename(columns={"mimotope": "peptide"})
    return df

  for col in [
      "variant_sequence",
      "variant_peptide",
      "mimotope_sequence",
      "sequence",
  ]:
    if col in df.columns:
      if "peptide" in df.columns:
        df = df.drop(columns=["peptide"])
      if "Peptide" in df.columns:
        df = df.rename(columns={"Peptide": "wt_peptide"})
      df = df.rename(columns={col: "peptide"})
      return df

  if "peptide" in df.columns:
    if "Peptide" in df.columns:
      df = df.rename(columns={"Peptide": "wt_peptide"})
    return df

  if "Peptide" in df.columns:
    df = df.rename(columns={"Peptide": "peptide"})

  return df


def extract_metric(
    df: pd.DataFrame, candidate_cols: list, default_val: float
) -> pd.Series:
  """Extracts the first valid numeric column from a candidate list."""
  for col in candidate_cols:
    if col in df.columns:
      series = pd.to_numeric(df[col], errors="coerce")
      if series.notna().any():
        return series.fillna(default_val)
  return pd.Series(default_val, index=df.index)


def compile_master_candidate_table(
    m05_file: Path,
    af3_compiled_file: Path,
    af3_pae_file: Path,
    af3_plddt_file: Path,
    af3_selected_file: Path,
    rosetta_summary_file: Path,
    rosetta_final_file: Path,
) -> pd.DataFrame:
  if not Path(m05_file).exists():
    raise FileNotFoundError(f"Module 05 input file missing: {m05_file}")

  df_master = pd.read_csv(m05_file)
  df_master = standardize_peptide_column(df_master)

  # 1. Filter out non-mimotope categories
  if "category" in df_master.columns:
    df_master["category_clean"] = (
        df_master["category"]
        .astype(str)
        .str.strip()
        .str.upper()
        .str.replace("_", "-")
    )
    df_master = df_master[
        df_master["category_clean"].isin(VALID_CATEGORIES)
    ].copy()
    logger.info(
        f"Filtered for valid categories (1-AA, 2-AA, PS-CPL): {len(df_master)}"
        " candidates remain."
    )

  # 2. Merge AF3 compiled results
  if Path(af3_compiled_file).exists():
    df_af3 = pd.read_csv(af3_compiled_file)
    df_af3 = standardize_peptide_column(df_af3)
    if "peptide" in df_af3.columns:
      overlap = [
          c
          for c in df_af3.columns
          if c != "peptide" and c in df_master.columns
      ]
      for c in overlap:
        if df_master[c].isna().all():
          df_master = df_master.drop(columns=[c])
      df_master = pd.merge(
          df_master, df_af3, on="peptide", how="left", suffixes=("", "_af3")
      )

  # 3. Merge AF3 PAE results
  if Path(af3_pae_file).exists():
    df_pae = pd.read_csv(af3_pae_file)
    df_pae = standardize_peptide_column(df_pae)
    if "peptide" in df_pae.columns:
      overlap = [
          c for c in df_pae.columns if c != "peptide" and c in df_master.columns
      ]
      for c in overlap:
        if df_master[c].isna().all():
          df_master = df_master.drop(columns=[c])
      df_master = pd.merge(
          df_master, df_pae, on="peptide", how="left", suffixes=("", "_pae")
      )

  # 4. Restrict to AF3 passed subset
  if Path(af3_selected_file).exists():
    df_af3_sel = pd.read_csv(af3_selected_file)
    df_af3_sel = standardize_peptide_column(df_af3_sel)
    if "peptide" in df_af3_sel.columns:
      af3_passed_peptides = set(df_af3_sel["peptide"].dropna().unique())
      df_master = df_master[
          df_master["peptide"].isin(af3_passed_peptides)
      ].copy()
      logger.info(
          f"Filtered using {af3_selected_file.name}: {len(df_master)}"
          " AF3-passed candidates remain."
      )
  elif "ranking_score" in df_master.columns:
    df_master = df_master[df_master["ranking_score"] >= 0.8].copy()
    logger.info(
        "Filtered using AF3 ranking_score >= 0.8:"
        f" {len(df_master)} candidates remain."
    )

  # 5. Merge Rosetta results
  if Path(rosetta_summary_file).exists():
    df_rosetta = pd.read_csv(rosetta_summary_file)
    df_rosetta = standardize_peptide_column(df_rosetta)
    if "peptide" in df_rosetta.columns:
      overlap = [
          c
          for c in df_rosetta.columns
          if c != "peptide" and c in df_master.columns
      ]
      for c in overlap:
        if df_master[c].isna().all():
          df_master = df_master.drop(columns=[c])
      df_master = pd.merge(
          df_master,
          df_rosetta,
          on="peptide",
          how="left",
          suffixes=("", "_rosetta"),
      )

  # Metric Extraction
  df_master["binding_affinity_nM"] = extract_metric(
      df_master,
      [
          "binding_affinity_nM",
          "HLA-A*02:01_Aff(nM)",
          "phla_affinity_nM",
          "Aff(nM)",
      ],
      default_val=100.0,
  )
  df_master["elution_rank"] = extract_metric(
      df_master,
      [
          "elution_rank",
          "HLA-A*02:01_%Rank_EL",
          "phla_elution_rank",
          "%Rank_EL",
      ],
      default_val=2.0,
  )
  df_master["half_life"] = extract_metric(
      df_master,
      [
          "half_life",
          "stability_score",
          "T1/2",
          "t12",
          "half_life_Thalf",
          "HLA-A*02:01_Thalf(h)",
      ],
      default_val=1.0,
  )
  df_master["rosetta_dG"] = extract_metric(
      df_master,
      ["rosetta_dG", "mean_dG", "best_dG", "mean_norm_dG"],
      default_val=-25.0,
  )

  df_master["af3_mean_pae"] = extract_metric(
      df_master,
      ["af3_mean_pae", "raw_inter_pae", "pae", "mean_pae", "af3_pae"],
      default_val=np.nan,
  )
  if df_master["af3_mean_pae"].isna().all():
    if (
        "ranking_score" in df_master.columns
        and df_master["ranking_score"].notna().any()
    ):
      rs = pd.to_numeric(df_master["ranking_score"], errors="coerce").fillna(
          0.5
      )
      df_master["af3_mean_pae"] = (1.0 - rs) * 20.0
    else:
      df_master["af3_mean_pae"] = 5.0
  else:
    df_master["af3_mean_pae"] = df_master["af3_mean_pae"].fillna(5.0)

  # Updated Tier 1 Selection Gate flag (Plot 1 quadrant: Rank <= 2.5%, Aff <= 500 nM)
  df_master["passed_all_gates"] = (
      (df_master["elution_rank"] <= 2.5)
      & (df_master["binding_affinity_nM"] <= 500.0)
  )

  logger.info(
      f"Compiled master table: {len(df_master)} AF3 structural candidates"
      f" ({df_master['passed_all_gates'].sum()} passed Plot 1 gates)."
  )
  return df_master
