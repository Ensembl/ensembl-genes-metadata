import logging
import logging.config
import os
import subprocess
import json
import re

def run_fastqc(fastqc_dir):
    # Get a list of all fastq files in the specified folder
    fastq_files = [file for file in os.listdir(fastqc_dir) if file.endswith((".fastq", ".gz"))]

    # Create output folder for FastQC results
    output_folder = os.path.join(fastqc_dir, "fastqc_results")
    os.makedirs(output_folder, exist_ok=True)

    for fastq_file in fastq_files:
        # Construct the full path to the input fastq file
        input_file_path = os.path.join(fastqc_dir, fastq_file)

        # Construct the full path to the output folder and file
        output_file_path = os.path.join(output_folder, fastq_file.replace(".fastq", "_fastqc"))

        # Run FastQC command using subprocess
        #cmd = ["LC_ALL=C singularity", "exec", "./fastqc_latest.sif", "fastqc", "--outdir", output_folder, input_file_path, "-q", "--extract"]
        cmd = "LC_ALL=C singularity exec ./fastqc_latest.sif fastqc --outdir {} {} -q --extract".format(output_folder, input_file_path)
        subprocess.run(cmd, shell=True)

        logging.info(f"FastQC completed for {fastq_file}. Results saved to {output_file_path}")

def parse_fastqc_summary(summary_path):
    # Read the summary.txt file and parse relevant information
    summary_data = {}
    with open(summary_path, 'r') as summary_file:
        for line in summary_file:
            test_results = line.strip().split('\t')
            if len(test_results) == 3:
                status, module, filename = test_results
                module = module.replace(" ","_").lower()
                summary_data[module] = status

    return summary_data


def parse_fastqc_data(fastqc_data_path):
    # Read the fastqc_data.txt file to get the total number of sequences and GC content
    fastqc_data = {}
    with open(fastqc_data_path, 'r') as f:
        data = f.read()
    # Regex patterns to extract total sequences and %GC
    total_sequences_pattern = r"Total Sequences\s+(\d+)"
    gc_content_pattern = r"%GC\s+(\d+)"
    # Search for total sequences
    total_sequences_match = re.search(total_sequences_pattern, data)
    if total_sequences_match:
        fastqc_data['total_sequences'] = int(total_sequences_match.group(1))
    # Search for %GC
    gc_content_match = re.search(gc_content_pattern, data)
    if gc_content_match:
        fastqc_data['gc_content'] = int(gc_content_match.group(1))

    return fastqc_data

def convert_to_json(folder_path):
    # Create a dictionary to store parsed summary information for each fastq file
    summary_data_dict = {}

    # Get a list of all fastq files in the specified folder
    fastq_files = [file for file in os.listdir(folder_path) if file.endswith(('.fastq', '.gz'))]

    for fastq_file in fastq_files:
        # Construct the full path to the summary.txt file
        summary_path = os.path.join(folder_path, 'fastqc_results', f'{fastq_file.replace(".fastq.gz", "_fastqc")}', 'summary.txt')
        fastqc_data_path = os.path.join(folder_path, 'fastqc_results', f'{fastq_file.replace(".fastq.gz", "_fastqc")}', 'fastqc_data.txt')

        # Parse the summary file and store the information in the dictionary
        summary_data = parse_fastqc_summary(summary_path)
        fastqc_data = parse_fastqc_data(fastqc_data_path)
        summary_data_dict[fastq_file] = summary_data, fastqc_data

    # Convert the dictionary to a JSON string
    json_data = json.dumps(summary_data_dict, indent=2)

    # Save the JSON data to a file
    output_json_path = os.path.join(folder_path, 'fastqc_results_summary.json')
    with open(output_json_path, 'w') as json_file:
        json_file.write(json_data)

    logging.info(f"Summary data converted to JSON. Results saved to {output_json_path}")

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    folder_path = "/hps/nobackup/flicek/ensembl/genebuild/swati/fastqc/version_12_test/"
    run_fastqc(folder_path)
    convert_to_json(folder_path)

