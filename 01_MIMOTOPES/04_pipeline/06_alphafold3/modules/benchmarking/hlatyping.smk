import os
import pandas as pd


include: os.path.join("rules", "benchmarks.smk")
rule all:
    """This is the entire NeoLR pipeline"""
    input:
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "star"], 
            module = ["star"])