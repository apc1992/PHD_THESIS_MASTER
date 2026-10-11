import os
from types import SimpleNamespace

class PathBuilder:
    def __init__(self, benchmark, directory, pipeline, module):
        """
        Initialize PathBuilder with a configuration dictionary.
        """
        self.benchmark = benchmark
        self.directory = directory
        self.pipeline = pipeline
        self.module = module

    def get_input_path(self, bam_id, threads, module, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join(self.directory, bam_id, self.pipeline, module, rule, str(threads), file)
        else:
            return os.path.join(self.directory, bam_id, self.pipeline, module, rule, file)

    def get_input_path_with_pipeline(self, bam_id, threads, pipeline, module, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join(self.directory, bam_id, pipeline, module, rule, str(threads), file)
        else:
            return os.path.join(self.directory, bam_id, pipeline, module, rule, file)
    
    def get_output_path(self, bam_id, threads, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join(self.directory, bam_id, self.pipeline, self.module, rule, str(threads), file)
        else:
            return os.path.join(self.directory, bam_id, self.pipeline, self.module, rule, file)

    def get_output_path_labeled(self, bam_id, threads, rule, label, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join(self.directory, bam_id, self.pipeline, self.module, rule, label, str(threads), file)
        else:
            return os.path.join(self.directory, bam_id, self.pipeline, self.module, rule, label, file)

    def wo_get_output_path(self, threads, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join(self.directory, self.pipeline, self.module, rule, str(threads), file)
        else:
            return os.path.join(self.directory, self.pipeline, self.module, rule, file)

    def get_benchmark_path(self, bam_id, threads, rule):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, bam_id + "_" + str(threads) + "_threads.txt")
        else:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, bam_id + "." + rule + ".benchmark.tsv")

    def fget_benchmark_path(self, bam_id, threads, rule, modulename):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join("benchmarks", self.pipeline, modulename, rule, bam_id + "_" + str(threads) + "_threads.txt")
        else:
            return os.path.join("benchmarks", self.pipeline, modulename, rule, bam_id + "." + rule + ".benchmark.tsv")

    def wo_get_benchmark_path(self, threads, rule):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, str(threads) + "_threads.txt")
        else:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, rule + ".benchmark.tsv")

    def get_log_paths(self, bam_id, threads, rule):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, bam_id + "_" + str(threads) + "_threads.log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, bam_id + "_" + str(threads) + "_threads.err")
        else:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, bam_id + ".log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, bam_id + ".err")
        return stdout, stderr

    def fget_log_paths(self, bam_id, threads, rule, modulename):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            stdout = os.path.join("logs", self.pipeline, modulename, rule, bam_id + "_" + str(threads) + "_threads.log")
            stderr = os.path.join("logs", self.pipeline, modulename, rule, bam_id + "_" + str(threads) + "_threads.err")
        else:
            stdout = os.path.join("logs", self.pipeline, modulename, rule, bam_id + ".log")
            stderr = os.path.join("logs", self.pipeline, modulename, rule, bam_id + ".err")
        return stdout, stderr

    def wo_get_log_paths(self, threads, rule):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        if self.benchmark:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, str(threads) + "_threads.log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, str(threads) + "_threads.err")
        else:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, rule + ".log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, rule + ".err")
        return stdout, stderr

    def initialize_writers(self):
        """
        Initialize all writers (output, benchmark, log) with shared configuration.
        """
        owriter = OutputWriter(self.benchmark, self.directory, self.pipeline, self.module)
        bwriter = BenchmarkWriter(self.benchmark, self.directory, self.pipeline, self.module)
        lwriter = LogWriter(self.benchmark, self.directory, self.pipeline, self.module)
        return {"output": owriter, "benchmark": bwriter, "log": lwriter}

class Writer(PathBuilder):
    def inputs(self, bam_id, threads, module, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_input_path(bam_id, threads, module, rule, file)

    def inputs_pmr(self, bam_id, threads, pipeline, module, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_input_path_with_pipeline(bam_id, threads, pipeline, module, rule, file)

    def outputs(self, bam_id, threads, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_output_path(bam_id, threads, rule, file)

    def outputs_label(self, bam_id, threads, rule, label, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_output_path(bam_id, threads, rule, file)

    def wos_outputs(self, threads, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        return self.wo_get_output_path(threads, rule, file)

    def collapse_outputs(self, bam_id, threads, rule, file, secondary_directory=""):
        if self.benchmark:
            return os.path.join(self.directory, bam_id, self.pipeline, self.module, rule, secondary_directory, str(threads), file)
        else:
            return os.path.join(self.directory, bam_id, self.pipeline, self.module, rule, secondary_directory, file)

    def benchmarks(self, bam_id, threads, rule):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_benchmark_path(bam_id, threads, rule)

    def fbenchmarks(self, bam_id, threads, rule, modulename):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        return self.fget_benchmark_path(bam_id, threads, rule, modulename)

    def wo_benchmarks(self, threads, rule):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        return self.wo_get_benchmark_path(threads, rule)

    def collapse_benchmarks(self, bam_id, threads, rule, secondary_directory=""):
        if self.benchmark:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, secondary_directory, bam_id + "_" + str(threads) + "_threads.txt")
        else:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, secondary_directory, bam_id + "." + rule + ".benchmark.tsv")

    def logs(self, bam_id, threads, rule):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_log_paths(bam_id, threads, rule)

    def flogs(self, bam_id, threads, rule, modulename):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        return self.fget_log_paths(bam_id, threads, rule, modulename)

    def wo_logs(self, threads, rule):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        return self.wo_get_log_paths(threads, rule)

    def collapse_logs(self, bam_id, threads, rule, secondary_directory=""):
        if self.benchmark:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + "_" + str(threads) + "_threads.log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + "_" + str(threads) + "_threads.err")
        else:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + ".log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + ".err")
        return stdout, stderr
        
    def io(self, rulename):
        """
        Generate Snakemake-compatible benchmark and log paths for a rule using {sample} and {threads} wildcards.
        """
        return SimpleNamespace(
            benchmark=self.benchmarks("{sample}", "{threads}", rulename),
            log_stdout=self.logs("{sample}", "{threads}", rulename)[0],
            log_stderr=self.logs("{sample}", "{threads}", rulename)[1],
        )

    def io_ns(self, wildcard_parameter, rulename):
        """
        Generate Snakemake-compatible benchmark and log paths for a rule using {sample} and {threads} wildcards.
        """
        return SimpleNamespace(
            benchmark=self.benchmarks(wildcard_parameter, "{threads}", rulename),
            log_stdout=self.logs(wildcard_parameter, "{threads}", rulename)[0],
            log_stderr=self.logs(wildcard_parameter, "{threads}", rulename)[1],
        )

    def fio(self, rulename, modulename):
        """
        Generate Snakemake-compatible benchmark and log paths for a rule using {sample} and {threads} wildcards.
        """
        return SimpleNamespace(
            benchmark=self.fbenchmarks("{sample}", "{threads}", rulename, modulename),
            log_stdout=self.flogs("{sample}", "{threads}", rulename, modulename)[0],
            log_stderr=self.flogs("{sample}", "{threads}", rulename, modulename)[1],
        )

    def mio(self, rulename):
        """
        Generate Snakemake-compatible benchmark and log paths for a rule using {sample} and {threads} wildcards.
        """
        return SimpleNamespace(
            benchmark=self.benchmarks("{msample}", "{threads}", rulename),
            log_stdout=self.logs("{msample}", "{threads}", rulename)[0],
            log_stderr=self.logs("{msample}", "{threads}", rulename)[1],
        )

    def ws_io(self, rulename):
        """
        Generate Snakemake-compatible benchmark and log paths for a rule using {sample} and {threads} wildcards.
        """
        return SimpleNamespace(
            benchmark=self.wo_benchmarks("{threads}", rulename),
            log_stdout=self.wo_logs("{threads}", rulename)[0],
            log_stderr=self.wo_logs("{threads}", rulename)[1],
        )

    def collapse_io(self, rulename, secondary_directory=""):
        """
        Generate collapsed benchmark and log paths for summary or aggregate rules.
        """
        return SimpleNamespace(
            benchmark=self.collapse_benchmarks("{sample}", "{threads}", rulename, secondary_directory),
            log_stdout=self.collapse_logs("{sample}", "{threads}", rulename, secondary_directory)[0],
            log_stderr=self.collapse_logs("{sample}", "{threads}", rulename, secondary_directory)[1],
        )
    

class OutputWriter(PathBuilder):
    def paths(self, bam_id, threads, rule, file):
        """
        Generate the output path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_output_path(bam_id, threads, rule, file)

class BenchmarkWriter(PathBuilder):
    def paths(self, bam_id, threads, rule):
        """
        Generate the benchmark path for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_benchmark_path(bam_id, threads, rule)

    def collapse_benchmarks(self, bam_id, threads, rule, secondary_directory=""):
        if self.benchmark:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, secondary_directory, bam_id + "_" + str(threads) + "_threads.txt")
        else:
            return os.path.join("benchmarks", self.pipeline, self.module, rule, secondary_directory, bam_id + "." + rule + ".benchmark.tsv")


class LogWriter(PathBuilder):
    def paths(self, bam_id, threads, rule):
        """
        Generate the log paths (stdout and stderr) for a given BAM ID, respecting the benchmark setting.
        """
        return self.get_log_paths(bam_id, threads, rule)

    def collapse_logs(self, bam_id, threads, rule, secondary_directory = ""):
        if self.benchmark:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + "_" + str(threads) + "_threads.log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + "_" + str(threads) + "_threads.err")
        else:
            stdout = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + ".log")
            stderr = os.path.join("logs", self.pipeline, self.module, rule, secondary_directory, bam_id + ".err")
        return stdout, stderr