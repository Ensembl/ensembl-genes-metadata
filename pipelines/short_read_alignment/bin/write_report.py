#!/usr/bin/env python3
"""Download the CSV file and enrich it with paths to BAM, CRAM, \
    and BigWig files based on the provided base directory."""
#pylint: disable=pointless-string-statement
import argparse
from pathlib import Path
import pandas as pd
import pymysql
from sqlalchemy import create_engine


def enrich_csv_with_paths(  # pylint: disable=too-many-statements, too-many-locals,too-many-arguments, too-many-positional-arguments,too-many-branches
    csv_path: Path,
    base_dir: Path,
    output_csv: Path,
    merge_tissue: bool,
    bam2bigWig: bool,#pylint: disable= invalid-name
    engine,  # pylint: disable=redefined-outer-name
):
    """
    Enrich the input CSV with paths to BAM, CRAM, and BigWig files
    based on the provided base directory and database connection."""
    base_dir = Path(base_dir).resolve()
    df = pd.read_csv(csv_path)
    accessions = df["run_accession"].dropna().unique().tolist()
    placeholders = ",".join(["%s"] * len(accessions))
    # Fetch run_accession -> biosample mappings
    query = f"""
    SELECT run_accession,sample_accession
    FROM run
    WHERE run_accession IN ({placeholders})
    """
    # Connect to the database
    # db_connection = connect_to_db(**db_config)
    # Clean the input text
    try:
        # if db_connection:
        # df = pd.read_sql(query, db_connection)
        biosample_df = pd.read_sql(query, engine, params=tuple(accessions))
        # print(f"Loaded {len(df)} rows.")
    except pymysql.MySQLError as e:
        print(f"Error connecting to MySQL: {e}")
    # Merge biosample into original DataFrame
    df = df.merge(biosample_df, on="run_accession", how="left")
    enriched_paths = []
    for idx, row in df.iterrows():
        taxon_id = str(row["taxon_id"])
        run_accession = str(row["run_accession"])
        platform = str(row["platform"])
        paired = bool(row["paired"])
        # print(paired)
        assembly_accession = str(row["assembly_accession"])
        biosample = str(row["sample_accession"])  # pylint: disable=unused-variable
        # sample_accession = str(row['sample_accession'])
        tissue = str(row["tissue"])
        genome_index_dir = Path(base_dir) / taxon_id / assembly_accession
        # print("Looking for .fna in:", list(genome_index_dir.glob("*genome.fna")))

        # genome_file = list(genome_index_dir.glob("*genome.fna"))[0]
        # relative_genome_file = genome_file.relative_to(base_dir)
        ##df["genome_file"] =  relative_genome_file
        # df.loc[idx, "genome_file"] = str(relative_genome_file)
        df.loc[idx, "genome_file"] = "genome.fna"
        # print(relative_genome_file)
        bam_path = Path(base_dir) / taxon_id / run_accession / "alignment"
        relative_splice_junction_file=[]
        if paired:
            star_index_genome = genome_index_dir / "Genome"
            # print(star_index_genome)
            relative_star_index_genome = star_index_genome.relative_to(base_dir)
            # print(relative_star_index_genome)
            # df["indexed_genome"] = relative_star_index_genome
            df.loc[idx, "indexed_genome"] = str(relative_star_index_genome)
            splice_junction_file = list(bam_path.glob("*SJ.out.tab"))[0]
            relative_splice_junction_file = splice_junction_file.relative_to(base_dir)
        else:
            minimap_files = list(genome_index_dir.glob("*.mmi"))

            if minimap_files:
                relative_minimap_index_genome = minimap_files[0].relative_to(base_dir)
                df.loc[idx, "indexed_genome"] = str(relative_minimap_index_genome)
            else:
                df.loc[idx, "indexed_genome"] = "" 
        # print(run_accession)
        #bam_path = Path(base_dir) / taxon_id / run_accession / "alignment"
        # print(f"{base_dir}/{taxon_id}/{run_accession}")
        # print(list(bam_path.glob("*.bam")))
        bam_file = list(bam_path.glob("*.bam"))[0]
        relative_bam_file = bam_file.relative_to(base_dir)
        # print(relative_bam_file)
        #splice_junction_file = list(bam_path.glob("*SJ.out.tab"))[0]
        #relative_splice_junction_file = splice_junction_file.relative_to(base_dir)
        # row["bam_file"] = relative_bam_file
        # row["splice_junction_file"] = relative_splice_junction_file
        df.loc[idx, "bam_file"] = str(relative_bam_file)
        df.loc[idx, "splice_junction_file"] = str(relative_splice_junction_file)
        # now consider the tissue

        if merge_tissue:
            tissue_dir = Path(base_dir) / taxon_id / platform / tissue / "alignment"
            print(str(tissue_dir))

            tissue_bams = list(tissue_dir.glob(f"*{tissue}.bam"))

            if tissue_bams:
                tissue_bam = tissue_bams[0]
                relative_tissue_bam = tissue_bam.relative_to(base_dir)
                df.loc[idx, "tissue_bam"] = str(relative_tissue_bam)
                enriched_paths.append(str(tissue_bam.resolve()))
            else:
                df.loc[idx, "tissue_bam"] = ""

            if bam2bigWig:
                # Stranded BigWigs
                file_fw = list(tissue_dir.glob("*forward_strand.bw"))
                file_rw = list(tissue_dir.glob("*reverse_strand.bw"))

                if file_fw or file_rw:
                    # Stranded data
                    df.loc[idx, "forward_bw_file"] = (
                        str(file_fw[0].relative_to(base_dir))
                        if file_fw else ""
                    )

                    df.loc[idx, "reverse_bw_file"] = (
                        str(file_rw[0].relative_to(base_dir))
                        if file_rw else ""
                    )

                    df.loc[idx, "bw_file"] = ""

                else:
                    # Unstranded data: one BigWig
                    bw_files = list(tissue_dir.glob("*.bw"))

                    df.loc[idx, "bw_file"] = (
                        str(bw_files[0].relative_to(base_dir))
                        if bw_files else ""
                    )

                    df.loc[idx, "forward_bw_file"] = ""
                    df.loc[idx, "reverse_bw_file"] = ""

        else:
            if bam2bigWig:
                file_fw = list(bam_path.glob("*forward_strand.bw"))
                file_rw = list(bam_path.glob("*reverse_strand.bw"))

                if file_fw or file_rw:
                    df.loc[idx, "forward_bw_file"] = (
                        str(file_fw[0].relative_to(base_dir))
                        if file_fw else ""
                    )
                    df.loc[idx, "reverse_bw_file"] = (
                        str(file_rw[0].relative_to(base_dir))
                        if file_rw else ""
                    )
                    df.loc[idx, "bw_file"] = ""

                else:
                    # Long-read / unstranded
                    bw_files = list(bam_path.glob("*.bw"))

                    df.loc[idx, "bw_file"] = (
                        str(bw_files[0].relative_to(base_dir))
                        if bw_files else ""
                    )

                    df.loc[idx, "forward_bw_file"] = ""
                    df.loc[idx, "reverse_bw_file"] = ""
            # if expected_path.exists():
        #    enriched_paths.append(str(expected_path.resolve()))
        # else:
        #    enriched_paths.append("")

    df.to_csv(output_csv, index=False)
    print(f"[✓] Output saved to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Enrich CSV with output paths.")
    parser.add_argument(
        "--host",
        type=str,
        default="mysql-ens-genebuild-prod-1",
        required=False,
        help="Host",
    )
    parser.add_argument(
        "--user", type=str, default="ensadmin", required=False, help="User"
    )
    parser.add_argument(
        "--password", type=str, default="ensembl", required=False, help="Password"
    )
    parser.add_argument(
        "--database",
        default="gb_transcriptomic_registry",
        type=str,
        required=False,
        help="Database",
    )
    parser.add_argument("--port", default=4527, type=int, required=False, help="Port")
    parser.add_argument("--csv", required=True, help="Input CSV path")
    parser.add_argument(
        "--base_dir", required=True, help="Base directory where outputs are stored"
    )
    parser.add_argument(
        "--out", default="enriched_output.csv", help="Output CSV filename"
    )
    parser.add_argument(
        "--merge_tissue",
        type=str,
        choices=["true", "false"],
        default="false",
        required=False,
        help="option to diplay report fot merged tissue",
    )
    parser.add_argument(
        "--bam2bigWig",
        type=str,
        choices=["true", "false"],
        default="true",
        required=False,
        help="option to diplay report fot merged tissue",
    )
    parser.add_argument(
        "--version",
        type=int,
        default=1,
        required=False,
        help="option to diplay report fot merged tissue",
    )

    args = parser.parse_args()
    engine = create_engine(
        f"mysql+pymysql://{args.user}:{args.password}@{args.host}:{args.port}/{args.database}"
    )
    enrich_csv_with_paths(
        args.csv, args.base_dir, args.out, args.merge_tissue, args.bam2bigWig, engine
    )
