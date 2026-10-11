import os
import pandas as pd

def link_bam(bam_file, sample_name, output_dir, pipeline):
    """
    Creates a symbolic link for a BAM file.
    """
    output_dir = os.path.join(output_dir, sample_name, pipeline, "preprocess", "fastq")
    os.makedirs(output_dir, exist_ok=True)
    link_path = os.path.join(output_dir, f"{sample_name}.bam")
    try:
        os.symlink(bam_file, link_path)
        print(f"Linked BAM: {bam_file} -> {link_path}")
    except FileExistsError:
        print(f"Link already exists: {link_path}")
    except Exception as e:
        print(f"Error linking BAM file: {e}")


def link_fastq(fastq_left, fastq_right, sample_name, output_dir, pipeline):
    """
    Creates symbolic links for FASTQ files.
    """
    output_dir = os.path.join(output_dir, sample_name,  pipeline, "preprocess", "fastq")
    os.makedirs(output_dir, exist_ok=True)
    left_link = os.path.join(output_dir, f"{sample_name}_R1.fq.gz")
    right_link = os.path.join(output_dir, f"{sample_name}_R2.fq.gz")

    try:
        os.symlink(fastq_left, left_link)
        os.symlink(fastq_right, right_link)
        print(f"Linked FASTQs: {fastq_left} -> {left_link}, {fastq_right} -> {right_link}")
    except FileExistsError:
        print(f"Links already exist: {left_link}, {right_link}")
    except Exception as e:
        print(f"Error linking FASTQ files: {e}")


def process_csv(csv_file, output_dir, pipeline):
    """
    Processes the CSV file to link BAM or FASTQ files. Raises an error if neither is found.
    """
    df = pd.read_csv(csv_file,sep=",",index_col=False)

    for _, row in df.iterrows():
        sample_name = row["sample_name"]
        bam_file = row["bam_file"]
        fastq_left = row["fastqfile1"]
        fastq_right = row["fastqfile2"]

        if pd.notna(bam_file):
            print(f"Processing BAM for sample: {sample_name}")
            link_bam(bam_file, sample_name, output_dir, pipeline)
        elif pd.notna(fastq_left) and pd.notna(fastq_right):
            print(f"Processing FASTQ for sample: {sample_name}")
            link_fastq(fastq_left, fastq_right, sample_name, output_dir, pipeline)
        else:
            raise ValueError(f"Input error for sample {sample_name}: No valid BAM or FASTQ files found.")


if __name__ == "__main__":
    import argparse

    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Link BAM or FASTQ files based on a CSV input.")
    parser.add_argument("csv_file", help="Path to the input CSV file.")
    parser.add_argument("pipeline", help="Name of pipeline")
    parser.add_argument("output_dir", help="Directory to store symbolic links.")
    args = parser.parse_args()

    # Process the CSV file
    process_csv(args.csv_file, args.output_dir, args.pipeline)

