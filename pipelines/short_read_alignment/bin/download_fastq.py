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


def md5_old(file_path: str) -> str:
    """Calculate the MD5 checksum of a file."""
    with open(file_path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def md5(file_path: str) -> str:
    """modified to read in chunks to avoid memory issues with large files."""
    h = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_ok(path: str) -> bool:
    """Check if a gzip file is valid."""
    try:
        subprocess.run(
            ["gzip", "-t", path],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except subprocess.CalledProcessError:
        return False


def download_file(url: str, dest: str) -> None:
    """Download a file using wget."""
    subprocess.run(["wget", "-q", "-O", dest, f"ftp://{url}"], check=True)


def file_is_valid(path: str, expected_md5: str) -> bool:
    """Check if a file exists, has the expected md5 checksum, and is a valid gzip file."""
    return os.path.exists(path) and md5(path) == expected_md5 and gzip_ok(path)


def ensure_file(
    url: str, final_path: str, expected_md5: str, max_retries: int = 3
) -> None:
    """
    Use final_path if it already exists and is valid.
    Otherwise download to a temp file, validate, and atomically replace final_path.
    """
    if file_is_valid(final_path, expected_md5):
        return

    tmp_path = final_path + ".tmp"

    for attempt in range(max_retries):
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

            download_file(url, tmp_path)

            if not gzip_ok(tmp_path):
                raise ValueError(f"gzip check failed for {tmp_path}")

            if md5(tmp_path) != expected_md5:
                raise ValueError(f"MD5 mismatch for {tmp_path}")

            os.replace(tmp_path, final_path)
            return

        except (subprocess.CalledProcessError, ValueError, FileNotFoundError) as e:
            print(
                f"Validation/download failed for {final_path} "
                f"on attempt {attempt + 1}/{max_retries}: {e}",
                file=sys.stderr,
            )
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            time.sleep(1)

    raise RuntimeError(
        f"ERROR: failed to obtain a valid file for "
        f"{final_path} after {max_retries} attempts"
    )


__version__ = "1.0.0"


def main():
    """Main function to handle command-line arguments and download files."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    # parser.add_argument("--taxon_id")
    # parser.add_argument("--gca")
    parser.add_argument("--run_accession")
    parser.add_argument("--url1")
    parser.add_argument("--url2", default="")
    parser.add_argument("--md5_1")
    parser.add_argument("--md5_2", default="")
    parser.add_argument("--outDir")
    # parser.add_argument("--genomeFile")
    parser.add_argument("--paired", default=False, action="store_true")
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    args = parser.parse_args()

    pair1_path = f"{args.run_accession}_1.fastq.gz"
    pair2_path = f"{args.run_accession}_2.fastq.gz" if args.paired else None
    ensure_file(args.url1, pair1_path, args.md5_1)
    if args.paired:
        ensure_file(args.url2, pair2_path, args.md5_2)


if __name__ == "__main__":
    main()
