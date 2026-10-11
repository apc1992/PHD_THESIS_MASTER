import glob
import json
import logging
from pathlib import Path
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

def parse_residue_pae(json_path: Path) -> dict:
    """Parses AF3 confidence JSON to extract peptide-MHC PAE metrics."""
    with open(json_path, "r") as f:
        data = json.load(f)

    if "pae" not in data or "token_chain_ids" not in data:
        return None

    pae_matrix = np.array(data["pae"])
    chain_ids = np.array(data["token_chain_ids"])
    unique_chains = np.unique(chain_ids)

    if len(unique_chains) < 2:
        return None

    # Separate peptide (shortest chain) from MHC/B2M chains
    chain_indices = [np.where(chain_ids == c)[0] for c in unique_chains]
    chain_indices.sort(key=len)

    pep_idx = chain_indices[0]
    mhc_idx = np.concatenate(chain_indices[1:])

    # Peptide-to-MHC & MHC-to-Peptide PAE
    pae_pep_to_mhc = pae_matrix[np.ix_(pep_idx, mhc_idx)]
    pae_mhc_to_pep = pae_matrix[np.ix_(mhc_idx, pep_idx)]

    residue_means = np.mean(pae_pep_to_mhc, axis=1).tolist()
    mean_inter_pae = float((np.mean(pae_pep_to_mhc) + np.mean(pae_mhc_to_pep)) / 2.0)
    mean_mhc_internal = float(np.mean(pae_matrix[np.ix_(mhc_idx, mhc_idx)]))
    aggregator_score = float(mean_inter_pae - mean_mhc_internal)

    folder_name = json_path.parent.name
    peptide_seq = folder_name.split("_")[-1].strip().upper()

    return {
        "job_name": folder_name,
        "peptide": peptide_seq,
        "raw_inter_pae": mean_inter_pae,
        "phla_aggregator_score": aggregator_score,
        "residue_pae_scores": json.dumps(residue_means),
    }

def main(
    results_dir: Path = Path("results"),
    compiled_csv: Path = Path("output/af3_compiled_results.csv"),
    output_csv: Path = Path("output/mimotope_screening_pae_compiled.csv"),
):
    logging.info("Scanning AF3 confidence files for PAE extraction...")
    search_pattern = str(results_dir / "**/*_confidences.json")
    json_files = [Path(f) for f in glob.glob(search_pattern, recursive=True)]

    valid_files = [
        f for f in json_files
        if "summary_confidences" not in f.name and "seed-" not in f.name and "sample-" not in f.name
    ]

    logging.info(f"Processing {len(valid_files)} rank-1 confidence files...")
    records = []
    for json_path in valid_files:
        try:
            res = parse_residue_pae(json_path)
            if res:
                records.append(res)
        except Exception as e:
            logging.error(f"Error parsing {json_path.name}: {e}")

    if not records:
        logging.warning("No valid PAE records extracted.")
        return

    pae_df = pd.DataFrame(records)

    # Merge category metadata
    if compiled_csv.exists():
        comp_df = pd.read_csv(compiled_csv)
        meta_cols = [c for c in ["job_name", "category", "Selection_Tier", "mimotope"] if c in comp_df.columns]
        if "job_name" in comp_df.columns:
            meta_df = comp_df[meta_cols].drop_duplicates(subset=["job_name"])
            pae_df = pd.merge(pae_df, meta_df, on="job_name", how="left")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    pae_df.to_csv(output_csv, index=False)
    logging.info(f"Successfully compiled PAE metrics ({len(pae_df)} entries) to: {output_csv}")

if __name__ == "__main__":
    main()
