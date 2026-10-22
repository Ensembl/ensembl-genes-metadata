#!/usr/bin/env python3
# pylint: disable=missing-module-docstring
#  See the NOTICE file distributed with this work for additional information
#  regarding copyright ownership.
#
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

"""Download fastq files from NCBI SRA using wget and check md5 checksums.
Args:
    --taxon_id: Taxon ID of the organism.
    --gca: Genome assembly accession.
    --run_accession: Run accession number.
    --url1: URL for the first read.
    --url2: URL for the second read (optional).
    --md5_1: Expected md5 checksum for the first read.
    --md5_2: Expected md5 checksum for the second read (optional).
    --outDir: Output directory to save downloaded files.
    --genomeFile: Genome file.
    --paired: Flag indicating if the reads are paired-end.
"""

import os
import hashlib
import subprocess
import time
import sys
import argparse


def md5(file_path: str) -> str:
    """Calculate the MD5 checksum of a file."""
    with open(file_path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def download_file(url: str, dest: str)-> None:
    """Download a file using wget."""
    subprocess.run(["wget", "-q", "-c", "-O", dest, f"ftp://{url}"], check=True)


def main():
    """Main function to handle command-line arguments and download files."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxon_id")
    #parser.add_argument("--gca")
    parser.add_argument("--run_accession")
    parser.add_argument("--url1")
    parser.add_argument("--url2", default="")
    parser.add_argument("--md5_1")
    parser.add_argument("--md5_2", default="")
    parser.add_argument("--outDir")
    #parser.add_argument("--genomeFile")
    parser.add_argument("--paired", default=False,  action="store_true")
    args = parser.parse_args()

    run_dir = os.path.join(args.outDir, args.taxon_id, args.run_accession)
    os.makedirs(run_dir, exist_ok=True)

    pair1_path = os.path.join(run_dir, f"{args.run_accession}_1.fastq.gz")
    pair2_path = os.path.join(run_dir, f"{args.run_accession}_2.fastq.gz") if args.paired else None

    existing_gz = [f for f in os.listdir(run_dir) if f.endswith(".gz")]

    # Skip download if files are present
    if (args.paired and len(existing_gz) >= 2) or (not args.paired and len(existing_gz) >= 1):
        #print(
        #    f"{args.taxon_id}\t{args.genomeFile}\t{args.gca}\t{args.run_accession}\t{pair1_path}\t{pair2_path}"
        #)
        return

    max_retries = 3
    for _ in range(max_retries):
        try:
            download_file(args.url1, pair1_path)
            if args.paired:
                download_file(args.url2, pair2_path)
        except subprocess.CalledProcessError:
            continue

        try:
            md5_1_actual = md5(pair1_path)
            if args.paired:
                md5_2_actual = md5(pair2_path)
        except FileNotFoundError:
            continue

        if args.paired:
            if md5_1_actual == args.md5_1 and md5_2_actual == args.md5_2:
                break
        else:
            if md5_1_actual == args.md5_1:
                break

        # Retry
        if os.path.exists(pair1_path):
            os.remove(pair1_path)
        if args.paired and os.path.exists(pair2_path):
            os.remove(pair2_path)
        time.sleep(1)
    else:
        print("ERROR: MD5 checksums do not match after retries", file=sys.stderr)
        sys.exit(1)

    #print(f"{args.taxon_id}\t{args.genomeFile}\t{args.gca}\t{args.run_accession}\t{pair1_path}\t{pair2_path}")


if __name__ == "__main__":
    main()
