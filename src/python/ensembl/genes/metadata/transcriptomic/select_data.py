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
"""Select assessed data from the database and create a report.
This script connects to a MySQL database, retrieves data related to sequencing runs,
applies quality checks using FastQC and STAR, and generates a report summarizing the results.
"""

import argparse
import re
import pandas as pd
import pymysql


def connect_to_db(host, user, password, database, port=3306) -> pymysql.connections.Connection:
    """
    Establishes a connection to the MySQL database.

    :param host: The hostname or IP address of the MySQL server.
    :param user: The username to connect to the database.
    :param password: The password for the user.
    :param db: The name of the database to connect to.
    :param port: The port of the MySQL server (default is 3306).
    :return: A pymysql connection object.
    """
    try:
        connection = pymysql.connect(host=host, user=user, password=password, database=database, port=port)
        return connection
    except pymysql.MySQLError as e:
        print(f"Error connecting to MySQL: {e}")
        return None


def fastqc_quality(row) -> bool:
    """Calculate FastQC quality based on the criteria.
    The function checks the FastQC quality criteria for each row and returns True if the criteria are met.


    Args:
        row (_type_): dataframe row

    Returns:
        bool: True if the criteria are met, False otherwise.
    """
    fastqc_criteria = [
        row["per_base_sequence_quality"] in ["PASS", "WARN"],
        row["overrepresented_sequences"] in ["PASS", "WARN"],
        row["per_base_n_content"] in ["PASS", "WARN"],
        row["per_sequence_quality_scores"] in ["PASS", "WARN"],
    ]
    return sum(fastqc_criteria) == 4


# Define STAR Quality classification
def star_quality(row) -> bool:
    """Calculate STAR quality based on the criteria.
    The function checks the STAR quality criteria for each row and returns True if the criteria are met.
    Args:
        row (_type_): dataframe row
    Returns:
        bool: True if the criteria are met, False otherwise.
    """
    return (
        row["uniquely_mapped_reads_percentage"] >= 50
        and row["percentage_reads_mapped_to_multiple_loci"] <= 30
        and row["percentage_reads_unmapped_too_short"] <= 20
    )


def check_fastqc_star_quality(df) -> pd.DataFrame:
    """Apply FastQC and STAR quality to each row
    Args:
        df (pd.DataFrame): DataFrame containing the data to be processed.
    Returns:
        pd.DataFrame: DataFrame with additional columns for FastQC and STAR quality.
    """
    # Apply FastQC and STAR quality to each row
    df["fastqc_pass"] = df.apply(fastqc_quality, axis=1)
    df["star_quality"] = df.apply(star_quality, axis=1)

    # Assign final quality labels

    df["final_fastqc_status"] = df.groupby("run_accession")["fastqc_pass"].transform(
        lambda x: "Failed" if not all(x) else "Passed"
    )
    df["final_star_status"] = df.groupby("run_accession")["star_quality"].transform(
        lambda x: "Failed" if not all(x) else "Passed"
    )
    df["passed_both"] = (df["final_fastqc_status"] == "Passed") & (df["final_star_status"] == "Passed")
    return df


def create_report(df) -> pd.DataFrame:
    """Create a report based on the DataFrame.
    This function summarizes the FastQC and STAR quality results for each taxon_id and run_accession.
    It calculates the number of runs that passed both quality checks and the percentage of runs that passed.
    The report is printed to the console.
    1. Group the DataFrame by taxon_id and run_accession.
    2. Aggregate the final_fastqc_status and final_star_status columns.
    3. Create a new column passed_both to indicate if both quality checks passed.
    4. Group the DataFrame by taxon_id and aggregate the results.
    5. Calculate the percentage of runs that passed both quality checks.
    6. Print the report to the console.
    7. Return the summarized DataFrame.
    8. The report includes the following columns:
        - taxon_id: The taxon ID.
        - passed_fastqc: The number of runs that passed FastQC.
        - passed_star: The number of runs that passed STAR.
        - passed_both: The number of runs that passed both quality checks.
        - total_runs: The total number of runs for the taxon_id.
        - percentage_passed_both: The percentage of runs that passed both quality checks.
    Args:
        df (pd.DataFrame): DataFrame containing the data to be processed.
    Returns:
        pd.DataFrame: DataFrame with the summarized report.
    """
    # Collapse to one row per run_accession
    run_summary = (
        df.groupby(["taxon_id", "run_accession"])
        .agg(
            final_fastqc_status=("final_fastqc_status", "first"),
            final_star_status=("final_star_status", "first"),
        )
        .reset_index()
    )
    run_summary["passed_both"] = (run_summary["final_fastqc_status"] == "Passed") & (
        run_summary["final_star_status"] == "Passed"
    )
    # REPORT
    taxon_pass_summary = (
        run_summary.groupby("taxon_id")
        .agg(
            passed_fastqc=("final_fastqc_status", lambda x: (x == "Passed").sum()),
            passed_star=("final_star_status", lambda x: (x == "Passed").sum()),
            passed_both=("passed_both", "sum"),
            total_runs=("run_accession", "nunique"),
        )
        .reset_index()
    )

    # Compute percentage outside of .agg()
    taxon_pass_summary["percentage_passed_both"] = (
        taxon_pass_summary["passed_both"] / taxon_pass_summary["total_runs"] * 100
    ).round(2)
    print(taxon_pass_summary.head())


def filter_data(df) -> list:
    df = df.drop_duplicates(subset="run_accession", keep="first")
    df_tissue = df[df["tissue_prediction"].notnull()].reset_index(drop=True)
    df_no_tissue = df[df["tissue_prediction"].isna()].reset_index(
        drop=True
    )  # .isnull()].reset_index(drop=True)
    print(df_tissue.head())
    print(df_no_tissue.head())
    # Apply FastQC and STAR quality to each row

    run_accessions = []
    # for i in range(len(df_tissue)):
    for current_tissue, group in df_tissue.groupby("tissue_prediction"):
        # Check if the run has passed both FastQC and STAR
        total_length = 0
        df_both = group[group["passed_both"]]
        if not df_both.empty:
            # for j in range(len(df_both)):
            for _, row in group[group["passed_both"]].iterrows():
                if total_length + row["total_sequences"] <= 250 * 10**6:
                    total_length += row["total_sequences"]
                    run_accessions.append(row["run_accession"])
                    print("both ", row["run_accession"], total_length)

            if total_length <= 250 * 10**6:
                df_fastq = group[(~group["passed_both"]) & (group["final_fastqc_status"] == "Passed")]
                for _, row in df_fastq.iterrows():
                    if total_length + row["total_sequences"] <= 250 * 10**6:
                        print("both fastqc ", row["run_accession"])
                        total_length += row["total_sequences"]
                        run_accessions.append(row["run_accession"])

            if total_length <= 250 * 10**6:
                df_star = group[(~group["passed_both"]) & (group["final_star_status"] == "Passed")]
                for _, row in df_star.iterrows():
                    if total_length + row["total_sequences"] <= 250 * 10**6:
                        print("both star ", row["run_accession"])
                        total_length += row["total_sequences"]
                        run_accessions.append(row["run_accession"])

        else:
            df_fastq = group[(~group["passed_both"]) & (group["final_fastqc_status"] == "Passed")]
            for j in range(len(df_fastq)):
                if total_length + row["total_sequences"] <= 250 * 10**6:
                    total_length += df_fastq["total_sequences"].iloc[j]
                    run_accessions.append(df_fastq["run_accession"].iloc[j])

            if total_length <= 250 * 10**6:
                df_star = group[(~group["passed_both"]) & (group["final_star_status"] == "Passed")]
                for j in range(len(df_star)):
                    if total_length + row["total_sequences"] <= 250 * 10**6:
                        total_length += df_star["total_sequences"].iloc[j]
                        run_accessions.append(df_star["run_accession"].iloc[j])

    for run, group in df_no_tissue.groupby("run_accession"):
        total_length = 0
        if group["passed_both"].any():
            for _, row in group[group["passed_both"]].iterrows():
                if total_length + row["total_sequences"] <= 250 * 10**6:
                    total_length += row["total_sequences"]
                    print(row["run_accession"])
                    run_accessions.append(row["run_accession"])

    return run_accessions


def main() -> None:
    """Module's entry-point."""
    parser = argparse.ArgumentParser(prog="llm_prediction.py", description="Predict tissue using LLMs")

    parser.add_argument("--taxon_id", default="10116", type=str, required=True, help="Taxonomy ID")
    parser.add_argument("--host", type=str, required=True, help="Host")
    parser.add_argument("--user", type=str, required=True, help="User")
    parser.add_argument("--password", type=str, required=True, help="Password")
    parser.add_argument(
        "--database",
        default="gb_transcriptomic_registry",
        type=str,
        required=False,
        help="Database",
    )
    parser.add_argument("--port", default=4527, type=int, required=False, help="Port")
    parser.add_argument("--file_name", type=str, required=False, help="Output file name")
    parser.add_argument(
        "--read_type",
        type=str,
        choices=["short", "long"],
        default="short",
        required=False,
        help="Read type short/long",
    )
    args = parser.parse_args()
    db_config = {
        "host": args.host,
        "user": args.user,
        "password": args.password,
        "database": args.database,
        "port": args.port,
    }

    # Load the rows to fix (e.g., 150k rows)
    query = ""
    if args.read_type == "short":
        query = (
            "SELECT r.taxon_id, r.run_accession, r.sample_accession, r.platform, r.paired, r.tissue_prediction, \
            d.file_name, d.file_url,d.md5,d.per_base_sequence_quality, d.per_sequence_quality_scores, \
            d.per_base_n_content, d.overrepresented_sequences,d.total_sequences,\
            a.uniquely_mapped_reads_percentage, \
            a.percentage_reads_mapped_to_multiple_loci,a.percentage_reads_unmapped_too_short \
            FROM run r INNER JOIN data_files d ON r.run_id=d.run_id INNER JOIN align a ON r.run_id=a.run_id \
            WHERE  r.qc_status='ALIGNED' AND r.paired=1 AND r.platform='ILLUMINA' AND r.taxon_id = "
            + args.taxon_id
        )
    elif args.read_type == "long":
        query = (
            "SELECT r.taxon_id, r.run_accession, r.sample_accession, r.platform, r.paired, r.tissue_prediction, \
            d.file_name, d.file_url,d.md5,d.per_base_sequence_quality, d.per_sequence_quality_scores, \
            d.per_base_n_content, d.overrepresented_sequences,d.total_sequences, \
                a.uniquely_mapped_reads_percentage, \
            a.percentage_reads_mapped_to_multiple_loci,a.percentage_reads_unmapped_too_short \
            FROM run r INNER JOIN data_files d ON r.run_id=d.run_id INNER JOIN align a ON r.run_id=a.run_id \
            WHERE  r.qc_status='ALIGNED' AND r.paired=0 AND r.platform='PACBIO_SMRT' AND r.taxon_id = "
            + args.taxon_id
        )
    # Connect to the database
    db_connection = connect_to_db(**db_config)
    # Clean the input text

    if db_connection:
        df = pd.read_sql(query, db_connection)
        print(f"Loaded {len(df)} rows.")
    df = check_fastqc_star_quality(df)

    #df["tissue_prediction"] = df["tissue_prediction"].apply(
    #    lambda x: re.search(r"^(.*)", str(x)).group(1) if re.search(r"^(.*)", str(x)) else x
    #)
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: (m.group(1) if (m := re.search(r"^(.*)", str(x))) else x)
    )
    df["tissue_prediction"] = (
        df["tissue_prediction"]
        .astype(str)
        .str.replace('"', "", regex=True)
        .replace("Answer: ", "", regex=True)
        .replace("\n", "", regex=True)
        .replace("Output: ", "", regex=True)
    )
    df["tissue_prediction"] = df["tissue_prediction"].apply(
        lambda x: re.sub(
            r"[, ]+",
            "_",  # replace space/comma with underscore
            re.sub(r"\b\d+\b|\bday\b|\b[a-zA-Z]\b", "", str(x)),  # remove numbers, "day", and single letters
        ).strip(
            "_"
        )  # remove leading/trailing underscores
    )
    df["tissue_prediction"] = df["tissue_prediction"].astype(str).str.lower()

    df.loc[
        df["tissue_prediction"].astype(str).str.contains("NONE", case=False, na=False), "tissue_prediction"
    ] = None
    df.loc[
        df["tissue_prediction"].astype(str).str.contains("nan", case=False, na=False), "tissue_prediction"
    ] = None

    create_report(df)
    print(df.head())
    # Remove duplicates based on 'run_accession' while keeping the first row for each
    df_original = df.copy()
    run_accessions = filter_data(df)

    # Build a regex pattern from run_accession list
    pattern = "|".join(re.escape(acc) for acc in run_accessions)

    # Filter df_original where file_name contains any of the run_accessions
    # df_final = df_original[df_original["file_name"].str.contains(pattern, na=False)].copy()

    df_final = df_original[df_original["run_accession"].isin(run_accessions[0:250])].copy()

    df_final.drop_duplicates(inplace=True)

    df_final.loc[:, "file_name"] = df_final["file_name"].astype(str) + ".fastq.gz"
    df_final.loc[:, "col1"] = 1
    df_final.loc[:, "col_1"] = -1
    df_final.loc[:, "col0"] = 0
    df_final.loc[:, "ENA"] = "ENA"

    df_final["predicted_tissue"] = df_final.apply(
        lambda row: (
            row["sample_accession"] if pd.isnull(row["tissue_prediction"]) else row["tissue_prediction"]
        ),
        axis=1,
    )

    with open(args.file_name, "w") as f:
        df_final.to_csv(
            f,
            sep="\t",
            index=False,
            columns=[
                "predicted_tissue",
                "run_accession",
                "col1",
                "file_name",
                "col_1",
                "col1",
                "col0",
                "ENA",
                "platform",
                "sample_accession",
                "file_url",
                "md5",
            ],
            header=False,
        )

    print("✅ CSV READY!!!!.")


if __name__ == "__main__":
    main()
