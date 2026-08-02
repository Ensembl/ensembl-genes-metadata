#!/usr/bin/env python3

import argparse
import json
import math
import sys
from pathlib import Path



def genome_stats(fasta_file: Path, read_length: int) -> dict:
    """Compute basic genome statistics from a FASTA file."""

    references = 0
    genome_length = 0

    with fasta_file.open() as fh:
        for line in fh:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                references += 1
            else:
                genome_length += len(line)

    if references == 0:
        raise ValueError(f"No FASTA records found in {fasta_file}")

    if genome_length == 0:
        raise ValueError(f"Genome length is zero in {fasta_file}")

    genomeSAindexNbases = min(
        14,
        int(math.log2(genome_length) / 2 - 1)
    )

    genomeChrBinNbits = min(
        18,
        int(math.log2(max(genome_length / references, read_length)))
    )

    return {
        "genomeSAindexNbases": genomeSAindexNbases,
        "genomeChrBinNbits": genomeChrBinNbits
    }

__version__ = "1.0.0"
def main():

    parser = argparse.ArgumentParser(
        description="Calculate STAR genome index parameters from a FASTA."
    )

    parser.add_argument(
        "fasta",
        type=Path,
        help="Genome FASTA (.fna)"
    )
    parser.add_argument(
    "--read_length",
    type=int,
    default=100,
    help="Expected read length (default: 100)"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args()

    if not args.fasta.exists():
        sys.exit(f"ERROR: File not found: {args.fasta}")

    stats = genome_stats(args.fasta, args.read_length)

    json.dump(stats, sys.stdout, indent=2)
    print()


if __name__ == "__main__":
    main()