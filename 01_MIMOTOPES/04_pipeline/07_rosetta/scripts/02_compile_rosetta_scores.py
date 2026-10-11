import glob
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

# Upstream metadata files to inspect for missing candidate categories/tiers
UPSTREAM_METADATA_SOURCES = [
    Path("../06_alphafold3/output/rosetta_candidates_selected.csv"),
    Path("../06_alphafold3/output/af3_compiled_results.csv"),
    Path("../05_in_silico_variants_selection/output/in_silico_variants_selected.csv"),
    Path("../04_data_summary/output/in_silico_variants_metrics.csv"),
    Path("output/input_libraries.csv"),
]


def load_upstream_metadata() -> pd.DataFrame:
  """Loads and unifies peptide metadata (category, tier, mimotope) across upstream pipeline modules."""
  meta_dfs = []

  for path in UPSTREAM_METADATA_SOURCES:
    if path.exists():
      try:
        df = pd.read_csv(path)
        if "peptide" in df.columns:
          # Standardize column naming across modules
          if (
              "rosetta_selection" in df.columns
              and "Selection_Tier" not in df.columns
          ):
            df["Selection_Tier"] = df["rosetta_selection"]

          cols = [
              c
              for c in ["peptide", "category", "Selection_Tier", "mimotope"]
              if c in df.columns
          ]
          sub_df = (
              df[cols]
              .dropna(subset=["peptide"])
              .drop_duplicates(subset=["peptide"])
          )
          meta_dfs.append(sub_df)
          logging.info(f"Retrieved candidate metadata from {path}")
      except Exception as e:
        logging.warning(f"Could not parse metadata from {path}: {e}")

  if not meta_dfs:
    return pd.DataFrame()

  # Merge metadata sources prioritizing earlier entries
  combined = meta_dfs[0]
  for df in meta_dfs[1:]:
    combined = pd.merge(
        combined, df, on="peptide", how="outer", suffixes=("", "_drop")
    )
    for col in ["category", "Selection_Tier", "mimotope"]:
      if col in combined.columns and f"{col}_drop" in combined.columns:
        combined[col] = combined[col].fillna(combined[f"{col}_drop"])
        combined.drop(columns=[f"{col}_drop"], inplace=True)

  return combined.drop_duplicates(subset=["peptide"])


def parse_score_file(file_path: Path) -> pd.DataFrame:
  """Parses a Rosetta score file (score.sc or interface_score.sc) into a DataFrame."""
  peptide = file_path.parent.name
  data_rows = []
  header = None

  with open(file_path, "r") as f:
    for line in f:
      parts = line.strip().split()
      if not parts:
        continue
      if parts[0] == "SCORE:" and "total_score" in parts:
        header = parts[1:]
      elif (
          parts[0] == "SCORE:" and header is not None and parts[1] != "total_score"
      ):
        data_rows.append(parts[1:])

  if not header or not data_rows:
    return pd.DataFrame()

  df = pd.DataFrame(data_rows, columns=header)
  df["peptide"] = peptide

  num_cols = [c for c in df.columns if c not in ["description", "peptide"]]
  for col in num_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

  return df


def main(
    rosetta_dir: Path = Path("rosetta_outputs"),
    output_dir: Path = Path("output"),
):
  if not rosetta_dir.exists():
    rosetta_dir = Path("results/rosetta_outputs")

  logging.info(f"Scanning Rosetta score files under: {rosetta_dir}")

  # 1. Parse energy score files
  score_files = list(rosetta_dir.glob("*/score.sc"))
  interface_files = list(rosetta_dir.glob("*/interface_score.sc"))

  if not score_files:
    raise FileNotFoundError(f"No score.sc files found under {rosetta_dir}")

  score_dfs = [
      parse_score_file(f) for f in score_files if not parse_score_file(f).empty
  ]
  if_dfs = [
      parse_score_file(f)
      for f in interface_files
      if not parse_score_file(f).empty
  ]

  relax_df = pd.concat(score_dfs, ignore_index=True)
  if_df = (
      pd.concat(if_dfs, ignore_index=True) if if_dfs else pd.DataFrame()
  )

  # 2. Merge decoy-level scores handling Rosetta suffix differences
  if not if_df.empty and "description" in if_df.columns:
    relax_df["merge_key"] = relax_df["description"].astype(str)
    if_df["merge_key"] = if_df["description"].astype(str)

    # If keys don't directly match, strip trailing worker suffix (_0001) added by InterfaceAnalyzer
    common_keys = set(relax_df["merge_key"]).intersection(
        set(if_df["merge_key"])
    )
    if not common_keys:
      if_df["merge_key"] = if_df["merge_key"].str.replace(
          r"_\d{4}$", "", regex=True
      )
      common_keys = set(relax_df["merge_key"]).intersection(
          set(if_df["merge_key"])
      )

    if common_keys:
      decoy_df = pd.merge(
          relax_df,
          if_df.drop(columns=["peptide", "description"], errors="ignore"),
          on="merge_key",
          how="left",
          suffixes=("", "_if"),
      ).drop(columns=["merge_key"])
      logging.info(
          "Successfully merged score.sc and interface_score.sc on decoy"
          " identifiers."
      )
    else:
      # Failsafe: match by peptide and decoy index position within group
      relax_df["decoy_idx"] = relax_df.groupby("peptide").cumcount()
      if_df["decoy_idx"] = if_df.groupby("peptide").cumcount()
      decoy_df = pd.merge(
          relax_df.drop(columns=["merge_key"], errors="ignore"),
          if_df.drop(
              columns=["peptide", "description", "merge_key"], errors="ignore"
          ),
          on=["peptide", "decoy_idx"],
          how="left",
          suffixes=("", "_if"),
      ).drop(columns=["decoy_idx"])
      logging.info(
          "Merged score.sc and interface_score.sc using peptide decoy position"
          " index."
      )
  else:
    decoy_df = relax_df

  # 3. Pull missing metadata directly from upstream module outputs
  meta_df = load_upstream_metadata()
  if not meta_df.empty:
    decoy_df = pd.merge(decoy_df, meta_df, on="peptide", how="left")
    logging.info(
        "Successfully annotated decoy scores with metadata from upstream"
        " modules."
    )

  # Apply safe fallbacks for required metadata fields
  if "category" not in decoy_df.columns:
    decoy_df["category"] = "Uncategorized"
  else:
    decoy_df["category"] = decoy_df["category"].fillna("Uncategorized")

  if "Selection_Tier" not in decoy_df.columns:
    decoy_df["Selection_Tier"] = "Tier 1"
  else:
    decoy_df["Selection_Tier"] = decoy_df["Selection_Tier"].fillna("Tier 1")

  # Save raw compiled decoy scores
  output_dir.mkdir(parents=True, exist_ok=True)
  compiled_csv = output_dir / "rosetta_compiled_scores.csv"
  decoy_df.to_csv(compiled_csv, index=False)
  logging.info(f"Saved compiled decoy scores -> {compiled_csv}")

  # 4. Build dynamic aggregation mapping
  agg_dict = {
      "category": ("category", "first"),
      "selection_tier": ("Selection_Tier", "first"),
      "n_decoys": ("total_score", "count"),
      "best_total_score": ("total_score", "min"),
      "mean_total_score": ("total_score", "mean"),
  }

  optional_metrics = {
      "best_dG": ("dG_separated", "min"),
      "mean_dG": ("dG_separated", "mean"),
      "mean_norm_dG": ("dG_separated/dSASAx100", "mean"),
      "mean_sc_value": ("sc_value", "mean"),
      "mean_dSASA": ("dSASA_int", "mean"),
      "mean_packstat": ("packstat", "mean"),
  }

  for out_col, (in_col, func) in optional_metrics.items():
    if in_col in decoy_df.columns:
      agg_dict[out_col] = (in_col, func)

  summary = decoy_df.groupby("peptide").agg(**agg_dict).reset_index()

  summary_csv = output_dir / "rosetta_candidate_summary.csv"
  summary.to_csv(summary_csv, index=False)
  logging.info(
      f"Saved candidate summary ({len(summary)} entries) -> {summary_csv}"
  )


if __name__ == "__main__":
  main()
