import os
import pandas as pd
import re
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

# I could simplify this code to produce for a module, but I wrote it in this format to run as a wildcard in a rule.

class BenchmarkConsolidator:
    """Class to consolidate benchmark files from a folder into a single CSV file."""

    def __init__(self, input_folder, output_file):
        self.input_folder = input_folder
        self.output_file = output_file

    def consolidate(self):
        """Reads benchmark files, consolidates them, and saves as a CSV."""
        folder_name = os.path.basename(os.path.normpath(self.input_folder))
        dataframes = []

        for filename in os.listdir(self.input_folder):
            if filename.endswith(".txt") and "_threads" in filename:
                # Use regex to extract sample and thread count
                match = re.match(r"^(.*)_(\d+)_threads\.txt$", filename)
                if not match:
                    continue  # Skip files that don't match the pattern

                sample = match.group(1)
                threads = match.group(2)

                file_path = os.path.join(self.input_folder, filename)
                df = pd.read_csv(file_path, sep="\t")

                df["rule"] = folder_name  # assuming folder_name is defined
                df["sample"] = sample
                df["threads"] = int(threads)

                dataframes.append(df)

        if not dataframes:
            print("No valid benchmark files found.")
            return

        combined_df = pd.concat(dataframes, ignore_index=True)
        columns = ["rule", "sample", "threads"] + [col for col in combined_df.columns if col not in ["rule", "sample", "threads"]]
        combined_df = combined_df[columns].sort_values(by=["sample", "threads"])
        os.makedirs(os.path.dirname(self.output_file), exist_ok=True)
        combined_df.to_csv(self.output_file, index=False)
        print(f"Consolidated CSV saved to {self.output_file}")
        print(folder_name)


class BenchmarkPlotter:
    """Class to generate benchmark plots from a consolidated CSV file."""

    def __init__(self, input_file, output_dir):
        self.input_file = input_file
        self.output_dir = output_dir

    def plot(self):
        """Generates and saves execution time and memory usage plots."""
        os.makedirs(self.output_dir, exist_ok=True)
        df = pd.read_csv(self.input_file)

        if "rule" not in df.columns or df["rule"].nunique() > 1:
            print("Error: CSV should contain a single unique 'rule'.")
            return

        rule_name = df["rule"].unique()[0]
        sns.set(style="whitegrid")

        plots = [
            {
                "y_column": "s",
                "y_label": "Time (seconds)",
                "title": f"{rule_name}: Execution Time vs. Number of Threads",
                "output_file": os.path.join(self.output_dir, f"{rule_name}_execution_time.png"),
            },
            {
                "y_column": "max_rss",
                "y_label": "Memory (MB)",
                "title": f"{rule_name}: Memory Usage (MB) vs. Number of Threads",
                "output_file": os.path.join(self.output_dir, f"{rule_name}_memory_usage.png"),
            }
        ]

        for plot in plots:
            plt.figure(figsize=(10, 6))
            sns.lineplot(
                data=df,
                x="threads",
                y=plot["y_column"],
                hue="sample",
                marker="o"
            )
            plt.title(plot["title"], fontsize=16)
            plt.xlabel("Number of Threads", fontsize=12)
            plt.ylabel(plot["y_label"], fontsize=12)
            plt.legend(title="Sample", fontsize=10)
            plt.xticks(fontsize=10)
            plt.yticks(fontsize=10)
            plt.tight_layout()
            os.makedirs(os.path.dirname(self.input_file), exist_ok=True)
            plt.savefig(plot["output_file"])
            print(f"Plot saved as {plot['output_file']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark Consolidation and Plotting Tool")
    subparsers = parser.add_subparsers(dest="command")

    # Subparser for consolidation
    consolidate_parser = subparsers.add_parser("consolidate", help="Consolidate benchmark files into a single CSV")
    consolidate_parser.add_argument("input_folder", type=str, help="Folder containing benchmark files")
    consolidate_parser.add_argument("output_file", type=str, help="Path to save the consolidated CSV file")

    # Subparser for plotting
    plot_parser = subparsers.add_parser("plot", help="Plot benchmark data from a consolidated CSV")
    plot_parser.add_argument("input_file", type=str, help="Path to the consolidated benchmark CSV file")
    plot_parser.add_argument("output_dir", type=str, help="Directory to save the output plots")

    args = parser.parse_args()

    if args.command == "consolidate":
        BenchmarkConsolidator(args.input_folder, args.output_file).consolidate()
    elif args.command == "plot":
        BenchmarkPlotter(args.input_file, args.output_dir).plot()
    else:
        parser.print_help()