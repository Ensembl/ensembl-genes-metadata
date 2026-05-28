import json
from pathlib import Path
import smtplib
from email.message import EmailMessage
import logging
from datetime import datetime
import argparse
import sys
import re
import requests
import pandas as pd
from io import StringIO


BIOPROJECT_ID_RE = re.compile(r"^PRJ[A-Z]{2}\d+$", re.IGNORECASE)


def setup_logging(log_folder: Path) -> None:
    """Set up logging to file and console."""
    log_folder.mkdir(parents=True, exist_ok=True)
    log_file = (
        log_folder / f"email_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    )

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )
    logging.info(f"Logging initialized. Log file: {log_file}")


def fetch_tables_from_api(
    api_base_url: str,
    bioproject_id=None,
    group_name=None,
):
    url = f"{api_base_url}/api/annotations/annotations/filter"

    payload = {
        "bioproject_id": bioproject_id,
        "group_name": group_name,
    }

    response = requests.post(url, json=payload, timeout=300)
    response.raise_for_status()
    return response.json()


def create_report_csv(project_key: str, project_info: dict, csv_folder: Path) -> Path:
    csv_folder.mkdir(parents=True, exist_ok=True)

    project_id = str(project_info["id"]).strip()

    if BIOPROJECT_ID_RE.match(project_id):
        bioproject_id = [project_id]
        group_name = None
    else:
        bioproject_id = None
        group_name = [project_id]

    api_result = fetch_tables_from_api(
        api_base_url="http://127.0.0.1:8000",
        bioproject_id=bioproject_id,
        group_name=group_name,
    )
    
    anno_wide_csv = api_result["downloadables_anno"]["anno_wide"]["csv"]
    anno_main = pd.read_csv(StringIO(anno_wide_csv))

    # Define statuses to exclude
    exclude_statuses = {"archive", "abandoned"}
    
    # Define renaming mapping
    status_mapping = {
        "check_busco": "low_protein_busco",
        "poor_genome_busco": "low_genome_busco",
        "completed": "pending_release",
        "handed_over": "pending_release",
        "coming_soon": "pending_release",
        "pre_released": "pending_release"
    }
    
    # Drop rows with unwanted statuses
    anno_main = anno_main[~anno_main["gb_status"].isin(exclude_statuses)]
    
    # Rename remaining statuses
    anno_main["gb_status"] = anno_main["gb_status"].replace(status_mapping)
    
    # Drop columns
    cols_to_keep = [
            "associated_project",
            "gca",
            "scientific_name",
            "common_name",
            "lowest_taxon_id",
            "release_date",
            "gb_status",
            "latest_annotated",
            "assembly_busco",
            "assembly_busco_lineage",
            "protein_busco",
            "protein_busco_lineage",
            "coding_genes"
            ]
    anno_main = anno_main[cols_to_keep]

    # Rename columns
    anno_main = anno_main.rename(
        columns={"lowest_taxon_id": "taxon_id"}
    )
    anno_main = anno_main.rename(
        columns={"latest_annotated": "latest_assembly_version_annotated"}
    )
    anno_main = anno_main.rename(
        columns={"coding_genes": "number_of_protein_coding_genes"}
    )
    anno_main = anno_main.rename(
        columns={"gb_status": "status"}
    )

    # Order table by status column
    anno_main = anno_main.sort_values("status")


    date_tag = datetime.now().strftime("%Y_%m")
    csv_path = csv_folder / f"{project_key}_{date_tag}.csv"
    anno_main.to_csv(csv_path, index=False)

    logging.info(f"CSV generated for project {project_key}: {csv_path}")
    return csv_path


def send_project_emails(
    base_folder: Path,
    from_email: str,
    smtp_host: str,
    smtp_port: int,
    body_text: str,
) -> None:
    """Send emails for each project in the base folder."""

    SCRIPT_FOLDER = Path(__file__).parent
    config_file = SCRIPT_FOLDER / "project_email_list.json"
    csv_folder = base_folder / "csv_reports"
    log_folder = base_folder / "logs"

    setup_logging(log_folder)

    logging.info(f"Starting multi-project email report (base folder: {base_folder})")

    if not config_file.exists():
        logging.error(f"Config file not found: {config_file}")
        sys.exit(1)

    with open(config_file, "r") as f:
        projects = json.load(f)  # nested dict: project_key -> {label, email, id}

    for project_key, project_info in projects.items():
        recipient = project_info["email"]
        project_label = project_info["label"]

        # Generate CSV for this project
        csv_file_path = create_report_csv(project_key, project_info, csv_folder)

        # Compose email
        msg = EmailMessage()
        msg["From"] = from_email
        msg["To"] = recipient
        msg["Subject"] = f"Monthly annotation report: {project_label}"
        msg.set_content(body_text)

        if csv_file_path.exists():
            with open(csv_file_path, "rb") as f:
                msg.add_attachment(
                    f.read(),
                    maintype="text",
                    subtype="csv",
                    filename=csv_file_path.name,
                )
            logging.info(f"Attached CSV for project {project_key}: {csv_file_path}")
        else:
            logging.warning(
                f"CSV for project {project_key} not found at {csv_file_path}. Sending email without attachment."
            )

        # Send email
        try:
            with smtplib.SMTP(smtp_host, smtp_port) as smtp:
                smtp.sendmail(
                msg["From"],
                [msg["To"]],
                msg.as_string()
                )
                logging.info(
                    f"Email successfully sent to {recipient} for project {project_label}"
                )
        except Exception as e:
            logging.error(
                f"Failed to send email to {recipient} for project {project_label}: {e}"
            )

    logging.info("Finished sending all project emails")


def main():
    parser = argparse.ArgumentParser(
        description="Send multi-project email reports with CSV attachments"
    )
    parser.add_argument(
        "-bf", "--base_folder", type=str, help="Base folder for CSV reports and logs"
    )
    args = parser.parse_args()

    BASE_FOLDER = Path(args.base_folder)

    # ---------- CONFIG ----------
    FROM_EMAIL = "genebuild-metadata@ebi.ac.uk"
    SMTP_HOST = "localhost"
    SMTP_PORT = 25
    BODY_TEXT = (
        BODY_TEXT
    ) = """Hello,

Please find attached the latest annotation report.

The `latest_assembly_version_annotated` indicates whether a newer version of the GCA chain is available but has not yet been annotated.

The `status` column indicates the current state of each genome annotation:

- live: Annotation is publicly available on beta.ensembl.org
- pending_release: Annotation is complete and  will be featured in an upcoming release on beta.ensembl.org. A draft version of the annotation can be found on https://ftp.ebi.ac.uk/pub/databases/ensembl/pre-release/ Please note that this may not be identical to the final annotation that will be made public on beta.ensembl.org
- in_progress: Annotation is currently underway
- low_genome_busco: Genome quality is low based on BUSCO assessment
- insufficient_data: Not enough data available to proceed with annotation
- low_protein_busco: Protein BUSCO score is low and/or differs from the genome BUSCO score by more than 10%

If you have any questions regarding these emails please send them to genebuild@ebi.ac.uk. This inbox is not monitored.

Regards,
Ensembl Genebuild
"""

    send_project_emails(
        base_folder=BASE_FOLDER,
        from_email=FROM_EMAIL,
        smtp_host=SMTP_HOST,
        smtp_port=SMTP_PORT,
        body_text=BODY_TEXT,
    )


if __name__ == "__main__":
    main()
