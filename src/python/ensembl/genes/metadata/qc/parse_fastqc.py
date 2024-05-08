import argparse
import logging
import logging.config
import os
import subprocess
import json
import re


def parse_fastqc_summary(summary_path:str):
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
    sequence_length_pattern = r"Sequence length\s+(\d+)"
    gc_content_pattern = r"%GC\s+(\d+)"
    # Search for total sequences
    total_sequences_match = re.search(total_sequences_pattern, data)
    if total_sequences_match:
        fastqc_data['total_sequences'] = int(total_sequences_match.group(1))
    # Search for %GC
    gc_content_match = re.search(gc_content_pattern, data)
    if gc_content_match:
        fastqc_data['gc_content'] = int(gc_content_match.group(1))
    sequence_length_match = re.search(sequence_length_pattern, data)
    if sequence_length_match:
        fastqc_data['sequence_length'] = int(sequence_length_match.group(1))

    return fastqc_data

def convert_to_json(folder_path):
    # Create a list to store parsed summary information for each fastq file
    data_files = []

    # Get a list of all fastq files in the specified folder
    fastq_files = [file for file in os.listdir(folder_path) if file.endswith(('.fastq', '.gz'))]

    for fastq_file in fastq_files:
        # Construct the full path to the summary.txt file
        summary_path = os.path.join(folder_path, 'fastqc_results', f'{fastq_file.replace(".fastq.gz", "_fastqc")}', 'summary.txt')
        fastqc_data_path = os.path.join(folder_path, 'fastqc_results', f'{fastq_file.replace(".fastq.gz", "_fastqc")}', 'fastqc_data.txt')

        # Parse the summary file and store the information in a dictionary
        summary_data = parse_fastqc_summary(summary_path)
        fastqc_data = parse_fastqc_data(fastqc_data_path)
        
        # Create a dictionary for each file's data
        file_data = {
            "file_id": "",  # You can add an ID here if needed
            "name": fastq_file,
            **summary_data,  # Merge summary data into the dictionary
            **fastqc_data     # Merge fastqc data into the dictionary
        }
        
        # Append file data to the list
        data_files.append(file_data)

    # Wrap the list of file data into a dictionary with key "data_files"
    output_data = {"data_files": data_files}
    
    # Convert the dictionary to a JSON string
    json_data = json.dumps(output_data, indent=2)

    # Save the JSON data to a file
    output_json_path = os.path.join(folder_path, 'fastqc_results_summary.json')
    with open(output_json_path, 'w') as json_file:
        json_file.write(json_data)

    logging.info(f"Summary data converted to JSON. Results saved to {output_json_path}")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse FASTQC output and create JSON with specified keys")
    parser.add_argument("--fastqc_results_path", required=True, type=str, help="Path to the FASTQC data file")
    parser.add_argument("--output_dir", required=True, type=str, help="Output dir")
    for key in keys_to_extract:
        parser.add_argument(f"--{key}", action="store_true", help=f"Include {key} in output JSON")
    # print(vars(parser))
    return parser.parse_args()

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.INFO)
    folder_path = "/hps/nobackup/flicek/ensembl/genebuild/swati/fastqc/version_12_test/"
    convert_to_json(folder_path)

