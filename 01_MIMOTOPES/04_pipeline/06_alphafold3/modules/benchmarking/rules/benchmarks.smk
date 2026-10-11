BASE_PREPROCESSING = "benchmarks/summary/"
COLLAPSE_PREPROCESSING = "benchmarks/summary_collapse/"

if config["benchmark_grouping"]:
    grouping = "benchmarks"
    rule gather_benchmarks:
        input:
            extracthifi_benchmarks=os.path.join("benchmarks", "{pipeline}", "{module}", "{benchmark_rules}")
        output:
            benchmark_summary = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}", "summary.csv")
        # shadow: "shallow"
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        group: grouping
        shell:
            """
            python scripts/utils/benchmark_handler.py consolidate {input.extracthifi_benchmarks} {output.benchmark_summary}
            """

    rule plot_benchmarks:
        input:
            benchmark_summary=rules.gather_benchmarks.output.benchmark_summary
        output:
            execution_time = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}",  "{benchmark_rules}_execution_time.png"),
            memory_usage = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_rules}_memory_usage.png")
        # shadow: "shallow"
        params: outdir = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}")
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        group: grouping
        shell:
            """
            python scripts/utils/benchmark_handler.py plot {input.benchmark_summary} {params.outdir} 
            """

    rule collapse_gather_benchmarks:
        input:
            extracthifi_benchmarks=os.path.join("benchmarks",  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}")
        output:
            benchmark_summary = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}", "summary.csv")
        # shadow: "shallow"
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        group: grouping
        shell:
            """
            python scripts/utils/benchmark_handler.py consolidate {input.extracthifi_benchmarks} {output.benchmark_summary}
            """

    rule collapse_plot_benchmarks:
        input:
            benchmark_summary=rules.collapse_gather_benchmarks.output.benchmark_summary
        output:
            execution_time = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}", "{benchmark_submodule}_execution_time.png"),
            memory_usage = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}", "{benchmark_submodule}_memory_usage.png")
        # shadow: "shallow"
        params: outdir = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}")
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        group: grouping
        shell:
            """
            python scripts/utils/benchmark_handler.py plot {input.benchmark_summary} {params.outdir} 
            """
else:
    rule gather_benchmarks:
        input:
            extracthifi_benchmarks=os.path.join("benchmarks", "{pipeline}", "{module}", "{benchmark_rules}")
        output:
            benchmark_summary = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}", "summary.csv")
        # shadow: "shallow"
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        shell:
            """
            python scripts/utils/benchmark_handler.py consolidate {input.extracthifi_benchmarks} {output.benchmark_summary}
            """

    rule plot_benchmarks:
        input:
            benchmark_summary=rules.gather_benchmarks.output.benchmark_summary
        output:
            execution_time = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}",  "{benchmark_rules}_execution_time.png"),
            memory_usage = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_rules}_memory_usage.png")
        # shadow: "shallow"
        params: outdir = os.path.join(BASE_PREPROCESSING, "{pipeline}", "{module}", "{benchmark_rules}")
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        shell:
            """
            python scripts/utils/benchmark_handler.py plot {input.benchmark_summary} {params.outdir} 
            """

    rule collapse_gather_benchmarks:
        input:
            extracthifi_benchmarks=os.path.join("benchmarks",  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}")
        output:
            benchmark_summary = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}", "summary.csv")
        # shadow: "shallow"
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        shell:
            """
            python scripts/utils/benchmark_handler.py consolidate {input.extracthifi_benchmarks} {output.benchmark_summary}
            """

    rule collapse_plot_benchmarks:
        input:
            benchmark_summary=rules.collapse_gather_benchmarks.output.benchmark_summary
        output:
            execution_time = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}", "{benchmark_submodule}_execution_time.png"),
            memory_usage = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}", "{benchmark_submodule}_memory_usage.png")
        # shadow: "shallow"
        params: outdir = os.path.join(COLLAPSE_PREPROCESSING,  "{pipeline}", "{module}", "{benchmark_rules}", "{benchmark_submodule}")
        threads: 1
        resources: mem_mb=2048, time_min=10, pipeline=lambda wildcards: wildcards.pipeline
        shell:
            """
            python scripts/utils/benchmark_handler.py plot {input.benchmark_summary} {params.outdir} 
            """