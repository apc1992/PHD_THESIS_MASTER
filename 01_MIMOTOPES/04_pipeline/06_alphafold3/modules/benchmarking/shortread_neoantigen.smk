import os
import pandas as pd

include: os.path.join("rules", "benchmarks.smk")
rule all:
    """This is the entire NeoLR pipeline"""
    input:
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "link_bam",
                # "link_fastq",
                "view",
                "sort",
                "fastq",
                "performance"], 
            module = ["preprocess"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "QC_byFastp"], 
            module = ["fastpQC"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "align_bySTAR_M100",
                "sam2bam_M100",
                "fetch_split_unmapRead_M100",
                "align_bySTAR_M1000",
                "sam2bam_M1000",
                "fetch_split_unmapRead_M1000",
                "index_and_merge_bam"], 
            module = ["alignment"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "stringtie_assembly_withGG",
                "stringtie_assembly_withoutGG",
                "spades_unmapRead_assembly",
                "trinity_unmapRead_assembly"], 
            module = ["assembly"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "gtf_compare"], 
            module = ["gtf_comparison"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "build_BWT_fr_SR"], 
            module = ["fmlrc2convert"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "fetch_fa_fr_novelTptGTF_combA",
                "run_fmlrc2"], 
            module = ["correction"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "translateBySixF_NovelTpt_combA",
                "rmdup_novelpep_combA"], 
            module = ["translation"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "fetch_NovelNmer_Ref"], 
            module = ["selection"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "repeatmasker_TEtpt"], 
            module = ["annotation"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "transcript2pep_TEtpt",
                "translateBySixF_TEtpt"], 
            module = ["fetch"], pipeline=["isosplit"]),
        expand(rules.plot_benchmarks.output.execution_time, 
            benchmark_rules = [
                "mergeAnduniq_pepfa",
                "blast2ORF1ORF2"], 
            module = ["blast"], pipeline=["isosplit"]),

