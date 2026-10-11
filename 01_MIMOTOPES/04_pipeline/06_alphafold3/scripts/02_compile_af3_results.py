import json
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def parse_summary_confidences(json_path: Path) -> dict:
    """Extracts confidence metrics from AF3 summary JSON."""
    with open(json_path, "r") as f:
        data = json.load(f)
    return {
        "iptm": data.get("iptm"),
        "ptm": data.get("ptm"),
        "ranking_score": data.get("ranking_score"),
    }

def main(
    results_dir: Path = Path("results"),
    metadata_csv: Path = Path("../05_in_silico_variants_selection/output/in_silico_variants_selected.csv"),
    output_csv: Path = Path("output/af3_compiled_results.csv"),
):
    logging.info("Starting AlphaFold 3 results compilation...")
    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")

    meta_df = pd.read_csv(metadata_csv)
    if "Peptide" in meta_df.columns:
        meta_df["Peptide"] = meta_df["Peptide"].astype(str).str.strip().str.upper()
    if "mimotope" in meta_df.columns:
        meta_df["mimotope"] = meta_df["mimotope"].astype(str).str.strip().str.upper()

    compiled_records = []
    summary_files = list(results_dir.glob("**/*_summary_confidences.json"))
    logging.info(f"Found {len(summary_files)} summary confidence files.")

    for summary_path in summary_files:
        folder_name = summary_path.parent.name
        parts = folder_name.split("_")
        extracted_seq = parts[-1].strip().upper()
        candidate_id = parts[-2] if len(parts) >= 2 else "Unknown"

        metrics = parse_summary_confidences(summary_path)
        metrics.update({
            "job_name": folder_name,
            "candidate_id": candidate_id,
            "folder_sequence": extracted_seq,
        })
        compiled_records.append(metrics)

    if not compiled_records:
        logging.warning("No completed AF3 results were found.")
        return

    af3_df = pd.DataFrame(compiled_records)

    # Match extracted folder sequence against 'mimotope' column if present, else 'Peptide'
    merge_key = "mimotope" if "mimotope" in meta_df.columns else "Peptide"
    merged_df = pd.merge(meta_df, af3_df, left_on=merge_key, right_on="folder_sequence", how="inner")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(output_csv, index=False)
    logging.info(f"Successfully compiled {len(merged_df)} runs into: {output_csv}")

if __name__ == "__main__":
    main()
