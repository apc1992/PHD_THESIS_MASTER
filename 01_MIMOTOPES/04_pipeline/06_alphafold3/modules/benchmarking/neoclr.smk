import os
import pandas as pd

include: os.path.join("rules", "benchmarks.smk")
rule all:
    """This is the entire NeoLR pipeline"""
    input:
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
            "refine"], 
            module = ["preprocessing"], pipeline=["neolr"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
            "isoseq_cluster",
            "cluster_bam_to_fasta",
            "cluster_bam_to_fastq",
            "summarize_cluster_counts",
            ], module =["clustering"], pipeline=["neolr"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
            "orfipy_translation",
            "clean_orfipy_output",
            "deduplicate_orfs",
            "clean_dedup_output",
            "orfs_counts_summary",
            ], module =["detectORF"], pipeline=["neolr"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
            "pbmm2",
            "index_aligned_bam",
            "alignment_metrics",
            "extract_mapped_reads",
            "extract_unmapped_reads",
            "extract_read_ids",
            ], module =["alignment"], pipeline=["neolr"]),
        # Add isoforms module
        expand(rules.collapse_plot_benchmarks.output.execution_time, 
            module =["sqanti3"],
            benchmark_rules = ["sqanti_qc", "sqanti_filter"],
            benchmark_submodule = ["isoseq_collapse"], pipeline=["neolr"]),
        # Missing a few from sqanti3
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
            "LRfusion"], 
            module = ["lrfusion"], pipeline=["neolr"])