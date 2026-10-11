import sys
import logging
from pathlib import Path
import pandas as pd

MODULE_DIR = Path(__file__).resolve().parent.parent
PIPELINE_DIR = MODULE_DIR.parent
sys.path.append(str(PIPELINE_DIR))

from lib.logging_utils import setup_logger
from tier_classifier import compute_composite_scores, assign_unified_tier, select_candidates_struct_first
from selection_plots import plot_selection_dashboard

def main():
    log_file = PIPELINE_DIR / "logs" / "05_in_silico_variants_selection.log"
    logger = setup_logger("05_selection", log_file)
    logger.info("=== INITIALIZING MODULE 05: IN SILICO VARIANTS SELECTION ===")

    input_file = PIPELINE_DIR / "04_data_summary" / "output" / "in_silico_variants_metrics.csv"
    if not input_file.exists():
        fallback_files = list((PIPELINE_DIR / "04_data_summary" / "output").glob("*_in_silico_variants_metrics.csv"))
        if fallback_files:
            input_file = fallback_files[0]
        else:
            logger.error(f"Input metric file missing at {input_file}")
            sys.exit(1)

    logger.info(f"Loading summary metrics from: {input_file}")
    df = pd.read_csv(input_file)

    # 1. Compute Composite Scores & Structural Tiers
    df = compute_composite_scores(df)
    df['Selection_Tier'] = df.apply(assign_unified_tier, axis=1)

    # 2. Select Top 50 Variants per Category (WT completely excluded)
    logger.info("Extracting top 50 variants per category...")
    selected_df, cat_dict = select_candidates_struct_first(df, target_per_cat=50, pos_pair_cap=2)

    # 3. Export Output Files
    out_dir = MODULE_DIR / "output"
    fig_dir = MODULE_DIR / "figures"
    out_dir.mkdir(exist_ok=True, parents=True)
    fig_dir.mkdir(exist_ok=True, parents=True)

    selected_df.to_csv(out_dir / "in_silico_variants_selected.csv", index=False)
    for cat_name, cat_df in cat_dict.items():
        cat_df.to_csv(out_dir / f"selected_{cat_name.replace('-', '')}.csv", index=False)

    # 4. Generate Publication Plots
    plot_selection_dashboard(df, selected_df, fig_dir / "selection_dashboard_unified.png")

    logger.info(f"Successfully exported {len(selected_df)} candidates across 1AA, 2AA, and PS-CPL arms.")
    logger.info("=== MODULE 05 COMPLETE ===")

if __name__ == "__main__":
    main()
