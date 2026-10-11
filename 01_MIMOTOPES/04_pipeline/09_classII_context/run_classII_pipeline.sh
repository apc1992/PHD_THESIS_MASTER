#!/usr/bin/env bash
set -euo pipefail

mkdir -p output results figures

CONTAINER_NETMHC2="/rsrch6/home/hema_bio-Malignan/adp_shukla/pipeline_resources/netmhc2pan/netmhc2pan-4.3b/netmhc2pan.4.3b.sif"
RAW_OUTPUT="results/netmhc2pan_raw_output.txt"

# Reset raw output file before appending results
> "${RAW_OUTPUT}"

echo "=== Running NetMHCIIpan 4.3b via Singularity across Allele Files ==="

# Array of generated allele list files (including DRB, DP, and DQ sub-chunks)
ALLELE_FILES=(
    "output/mhc2_alleles_drb.txt"
    "output/mhc2_alleles_dp.txt"
    output/mhc2_alleles_dq_part*.txt
)

for allele_file in "${ALLELE_FILES[@]}"; do
    if [[ -f "${allele_file}" ]]; then
        ALLELES=$(cat "${allele_file}")
        echo "Processing ${allele_file} (${#ALLELES} chars)..."

        singularity exec \
          -B "$(pwd)" \
          -B /rsrch4:/rsrch4 \
          -B /rsrch5:/rsrch5 \
          -B /rsrch6:/rsrch6 \
          -B /rsrch9:/rsrch9 \
          "${CONTAINER_NETMHC2}" \
          /opt/netMHCIIpan-4.3/netMHCIIpan \
          -inptype 1 \
          -f output/netmhc2pan_input_peptides.txt \
          -BA \
          -a "${ALLELES}" \
          >> "${RAW_OUTPUT}"
    fi
done

echo "=== Step 4: Parsing Predictions ==="
python3 scripts/03_parse_netmhc2pan_results.py

echo "=== Workflow Complete ==="
