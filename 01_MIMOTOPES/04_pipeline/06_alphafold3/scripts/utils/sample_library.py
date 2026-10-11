import sys
import os
sys.path.append(os.getcwd())
import pandas as pd
from typing import Dict, Optional
from scripts.utils import symbolic_handler as symbolic_handler

class SampleLibrary:
    def __init__(self, csv_path: str, pipeline, link_dir: Optional[str] = None):
        self.csv_path = csv_path
        self.link_dir = link_dir
        self.pipeline = pipeline
        self.samples_table = pd.read_csv(csv_path, sep=",")

        if link_dir:
            symbolic_handler.process_csv(csv_path, link_dir, pipeline)

        self.bam_inputs = self._extract_bam_inputs()
        self.fastq_inputs = self._extract_fastq_inputs()

    def _extract_bam_inputs(self) -> Dict[str, str]:
        return {
            row['sample_name']: row['bam_file']
            for _, row in self.samples_table.iterrows()
            if pd.notna(row.get('bam_file'))
        }

    def _extract_fastq_inputs(self) -> Dict[str, Dict[str, str]]:
        return {
            row['sample_name']: {
                "fastq1": row['fastqfile1'],
                "fastq2": row['fastqfile2']
            }
            for _, row in self.samples_table.iterrows()
            if pd.notna(row.get('fastqfile1')) and pd.notna(row.get('fastqfile2'))
        }

    def get_bam_input(self, sample: str, use_symlinked: bool = False) -> Optional[str]:
        path = self.bam_inputs.get(sample)
        if use_symlinked and self.link_dir and path:
            return f"{sample}/{self.link_dir}/preprocess/{sample}.bam"
        return path

    def get_fastq_input(self, sample: str, use_symlinked: bool = False) -> Optional[Dict[str, str]]:
        paths = self.fastq_inputs.get(sample)
        if use_symlinked and self.link_dir and paths:
            return {
                "fastq1": f"{sample}/{self.link_dir}/preprocess/fastq/{sample}_R1.fq.gz",
                "fastq2": f"{sample}/{self.link_dir}/preprocess/fastq/{sample}_R2.fq.gz"
            }
        return paths

    def all_sample_keys(self) -> list:
        return list(set(self.bam_inputs.keys()) | set(self.fastq_inputs.keys()))
