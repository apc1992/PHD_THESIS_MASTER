#!/usr/bin/env python3
"""main.py

CLI orchestrator for Module 08: Final Candidates Selection.
"""

import argparse
from pathlib import Path
import sys

# Dynamic path configuration
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(PIPELINE_ROOT))

from lib.logging_utils import setup_logger
from step01_compile_master_table import compile_master_candidate_table
from step02_rank_candidates import calculate_zscores_and_rank
from step03_publication_plots import plot_selection_funnel, plot_zscore_breakdown


def parse_args():
  parser = argparse.ArgumentParser(
      description="Module 08: Final Candidate Selection & Multi-Pillar Prioritization"
  )
  parser.add_argument(
      "--m05-file",
      type=Path,
      default=PIPELINE_ROOT
      / "05_in_silico_variants_selection/output/in_silico_variants_selected.csv",
      help="Path to Module 05 selection output CSV",
  )
  parser.add_argument(
      "--af3-dir",
      type=Path,
      default=PIPELINE_ROOT / "06_alphafold3/output",
      help="Path to Module 06 AlphaFold 3 output directory",
  )
  parser.add_argument(
      "--rosetta-dir",
      type=Path,
      default=PIPELINE_ROOT / "07_rosetta/output",
      help="Path to Module 07 Rosetta output directory",
  )
  parser.add_argument(
      "--output-dir",
      type=Path,
      default=PIPELINE_ROOT / "08_final_candidates_selection/output",
      help="Output directory for generated CSV tables",
  )
  parser.add_argument(
      "--figures-dir",
      type=Path,
      default=PIPELINE_ROOT / "08_final_candidates_selection/figures",
      help="Output directory for publication figures",
  )
  return parser.parse_args()


def main():
  args = parse_args()
  args.output_dir.mkdir(parents=True, exist_ok=True)
  args.figures_dir.mkdir(parents=True, exist_ok=True)

  log_file = PIPELINE_ROOT / "logs/08_final_candidates_selection.log"
  log_file.parent.mkdir(parents=True, exist_ok=True)

  logger = setup_logger("08_final_candidates", log_file)
  logger.info("=== Starting Module 08: Final Candidates Selection ===")

  try:
    # 1. Compile Master Dataset
    df_master = compile_master_candidate_table(
        m05_file=args.m05_file,
        af3_compiled_file=args.af3_dir / "af3_compiled_results.csv",
        af3_pae_file=args.af3_dir / "mimotope_screening_pae_compiled.csv",
        af3_plddt_file=args.af3_dir / "peptide_anchor_stability.csv",
        af3_selected_file=args.af3_dir / "rosetta_candidates_selected.csv",
        rosetta_summary_file=args.rosetta_dir / "rosetta_candidate_summary.csv",
        rosetta_final_file=args.rosetta_dir / "rosetta_final_candidates.csv",
    )

    # 2. Compute Z-Scores and Rank Candidates
    df_ranked_master, df_ranked_validation = calculate_zscores_and_rank(
        df_master
    )

    # Export annotated tables
    master_out = args.output_dir / "master_150_candidates_annotated.csv"
    validation_out = args.output_dir / "final_candidates_for_validation.csv"

    df_ranked_master.to_csv(master_out, index=False)
    df_ranked_validation.to_csv(validation_out, index=False)
    logger.info(f"Exported master table: {master_out}")
    logger.info(f"Exported validation leads table: {validation_out}")

    # 3. Generate Publication Plots
    plot_selection_funnel(df_ranked_master, args.figures_dir)
    plot_zscore_breakdown(df_ranked_validation, args.figures_dir, top_n=25)

    logger.info("=== Module 08 Execution Completed Successfully ===")

  except Exception as e:
    logger.error(f"Module 08 failed with error: {str(e)}", exc_info=True)
    sys.exit(1)


if __name__ == "__main__":
  main()
