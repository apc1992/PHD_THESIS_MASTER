#!/usr/bin/env python3
import os
import re
import pandas as pd

def parse_4digit_allele(val):
    if pd.isna(val):
        return None
    val_str = str(val).strip()
    match = re.search(r'([A-Z0-9]+)[*:]?(\d{2,3}):(\d{2})', val_str)
    if match:
        _, f1, f2 = match.groups()
        return f"{f1}{f2}"
    return None

def chunk_alleles(allele_list, max_len=800):
    """Splits allele list into comma-separated strings under max_len chars."""
    chunks = []
    current_chunk = []
    current_len = 0
    for a in allele_list:
        if current_len + len(a) + 1 > max_len:
            chunks.append(','.join(current_chunk))
            current_chunk = [a]
            current_len = len(a)
        else:
            current_chunk.append(a)
            current_len += len(a) + 1
    if current_chunk:
        chunks.append(','.join(current_chunk))
    return chunks

def main():
    os.makedirs('output', exist_ok=True)
    hla_df = pd.read_csv('leukopaks_hla_typing.csv', index_col=0)
    
    drb_alleles = set()
    dp_alleles = set()
    dq_alleles = set()

    for col in hla_df.columns:
        # 1. DRB1
        for drb_row in ['DRB1_1', 'DRB1_2']:
            if drb_row in hla_df.index:
                digits = parse_4digit_allele(hla_df.loc[drb_row, col])
                if digits:
                    drb_alleles.add(f"DRB1_{digits}")

        # 2. DP Heterodimers
        dpa_set = {parse_4digit_allele(hla_df.loc[r, col]) for r in ['DPA1_1', 'DPA1_2'] if r in hla_df.index} - {None}
        dpb_set = {parse_4digit_allele(hla_df.loc[r, col]) for r in ['DPB1_1', 'DPB1_2'] if r in hla_df.index} - {None}
        for a in dpa_set:
            for b in dpb_set:
                dp_alleles.add(f"HLA-DPA1{a}-DPB1{b}")

        # 3. DQ Heterodimers
        dqa_set = {parse_4digit_allele(hla_df.loc[r, col]) for r in ['DQA_1', 'DQA_2'] if r in hla_df.index} - {None}
        dqb_set = {parse_4digit_allele(hla_df.loc[r, col]) for r in ['DQB_1', 'DQB_2'] if r in hla_df.index} - {None}
        for a in dqa_set:
            for b in dqb_set:
                dq_alleles.add(f"HLA-DQA1{a}-DQB1{b}")

    loci_map = {
        'drb': sorted(list(drb_alleles)),
        'dp': sorted(list(dp_alleles)),
        'dq': sorted(list(dq_alleles))
    }

    for locus, alleles in loci_map.items():
        # Write primary locus file
        full_str = ','.join(alleles)
        out_file = f'output/mhc2_alleles_{locus}.txt'
        with open(out_file, 'w') as f:
            f.write(full_str)
        print(f"[{locus.upper()}] Exported {len(alleles)} alleles ({len(full_str)} chars) -> {out_file}")

        # If locus exceeds 800 chars (e.g. DQ at 1011 chars), export safe sub-chunks
        chunks = chunk_alleles(alleles, max_len=800)
        if len(chunks) > 1:
            for idx, chunk in enumerate(chunks, 1):
                chunk_file = f'output/mhc2_alleles_{locus}_part{idx}.txt'
                with open(chunk_file, 'w') as f:
                    f.write(chunk)
                print(f"  └─ Sub-chunk {idx}: {len(chunk)} chars -> {chunk_file}")

if __name__ == '__main__':
    main()
