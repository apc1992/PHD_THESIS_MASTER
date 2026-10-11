import os
import pandas as pd

include: os.path.join("rules", "benchmarks.smk")
rule all:
    """This is the entire NeoLR pipeline"""
    input:
        expand(rules.plot_benchmarks.output.execution_time, benchmark_rules=[
            "extracthifi",
            "lima",
            "segmentation"], module = ["preprocessing"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, benchmark_rules=[
            "hifi_quality_metrics",
            "quality_control_report"], module = ["qualitycontrol"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, benchmark_rules=[
            "combined_demultiplexing_report",
            "combined_hifi_report",
            "combined_lima_report",
            "combined_segmentation_report",
            "merge_pdf"], module = ["summary"], pipeline=["isosplit"])