import json
import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

# Extracellular HLA-A*02:01 (275 AA) & Mature B2M (99 AA)
HLA_A0201_SEQ = (
    "GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQ"
    "THRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTK"
    "HKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMTHHAVSDHEATLRCWALSFYPAEITLT"
    "WQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWE"
)
B2M_SEQ = (
    "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTP"
    "TEKDEYACRVNHVTLSQPKIVKWDRDM"
)

def build_af3_payload(job_name: str, peptide_seq: str) -> dict:
    """Constructs the standard AF3 heterotrimer JSON schema."""
    return {
        "name": job_name,
        "modelSeeds": [1],
        "sequences": [
            {"protein": {"id": "A", "sequence": HLA_A0201_SEQ}},
            {"protein": {"id": "B", "sequence": B2M_SEQ}},
            {"protein": {"id": "C", "sequence": peptide_seq}},
        ],
        "dialect": "alphafold3",
        "version": 1,
    }

def main(
    input_csv: Path = Path("../05_in_silico_variants_selection/output/in_silico_variants_selected.csv"),
    output_dir: Path = Path("generated_jsons"),
    manifest_out: Path = Path("output/af3_job_manifest.csv"),
):
    logging.info("Initializing AF3 JSON payload generation...")
    if not input_csv.exists():
        raise FileNotFoundError(f"Input candidate sheet not found: {input_csv}")

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_out.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_csv)
    if "Peptide" not in df.columns:
        raise KeyError("Column 'Peptide' missing from input CSV.")

    unique_candidates = df.drop_duplicates(subset=["Peptide"]).reset_index(drop=True)
    manifest_records = []

    for idx, row in unique_candidates.iterrows():
        peptide = str(row["Peptide"]).strip().upper()
        candidate_id = f"Candidate_{idx+1:03d}"
        job_name = f"HLA_A0201_B2M_{candidate_id}_{peptide}"
        json_path = output_dir / f"{job_name}.json"

        payload = build_af3_payload(job_name, peptide)
        with open(json_path, "w") as f:
            json.dump(payload, f, indent=2)

        manifest_records.append({
            "candidate_id": candidate_id,
            "peptide": peptide,
            "category": row.get("category", "Unknown"),
            "selection_tier": row.get("Selection_Tier", "Unknown"),
            "json_path": str(json_path),
        })

    manifest_df = pd.DataFrame(manifest_records)
    manifest_df.to_csv(manifest_out, index=False)
    logging.info(f"Generated {len(manifest_df)} AF3 target payloads in '{output_dir}'")

if __name__ == "__main__":
    main()
