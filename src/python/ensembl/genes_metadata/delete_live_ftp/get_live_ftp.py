"""Discover live pre-release GCA directories and write a deletion manifest."""
import argparse
import ftplib
import logging
import sys
from datetime import datetime
from pathlib import Path

import pymysql
import pymysql.cursors

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

DEFAULT_DB = dict(
    host="mysql-ens-genebuild-prod-1",
    port=4527,
    user="ensro",
    password="",
    database="gb_assembly_metadata",
    cursorclass=pymysql.cursors.DictCursor,
)
FTP_HOST = "ftp.ebi.ac.uk"
FTP_BASE_DIR = "/pub/databases/ensembl/pre-release"


def get_latest_status_by_gca(db_config: dict) -> dict[str, str]:
    """Return the latest metadata status for each GCA accession."""
    log.info("Connecting to database %s@%s", db_config["database"], db_config["host"])
    conn = pymysql.connect(**db_config)
    try:
        with conn.cursor() as cur:
            cur.execute(
                """
                SELECT status.gca_accession, status.gb_status
                FROM genebuild_status AS status
                JOIN (
                    SELECT gca_accession,
                           MAX(genebuild_status_id) AS latest_id
                    FROM genebuild_status
                    WHERE gca_accession IS NOT NULL
                    GROUP BY gca_accession
                ) AS latest
                  ON latest.latest_id = status.genebuild_status_id
                """
            )
            rows = cur.fetchall()
    finally:
        conn.close()

    statuses = {
        row["gca_accession"].strip(): row["gb_status"] for row in rows
    }
    log.info(
        "Found latest statuses for %d GCA accessions (%d live).",
        len(statuses),
        sum(status == "live" for status in statuses.values()),
    )
    return statuses


def connect_ftp(host: str) -> ftplib.FTP:
    log.info("Connecting to FTP %s", host)
    ftp = ftplib.FTP(host, timeout=30)
    ftp.login()
    ftp.set_pasv(True)
    return ftp


def list_ftp_directory(ftp: ftplib.FTP, path: str) -> list[str]:
    try:
        return ftp.nlst(path)
    except ftplib.error_perm as exc:
        if "550" in str(exc):
            return []
        raise


def find_gca_directories_on_ftp(
    ftp: ftplib.FTP,
    base_dir: str,
    latest_status_by_gca: dict[str, str],
) -> list[dict[str, str]]:
    """Find pre-release GCA directories whose latest status is live."""
    records = []
    log.info("Scanning FTP base directory: %s", base_dir)

    for species_path in list_ftp_directory(ftp, base_dir):
        species = species_path.rstrip("/").split("/")[-1]
        for entry_path in list_ftp_directory(ftp, species_path):
            entry_name = entry_path.rstrip("/").split("/")[-1]
            if not entry_name.startswith("GCA_"):
                continue

            status = latest_status_by_gca.get(entry_name)
            if status != "live":
                log.debug("KEEP: %s (latest status: %s)", entry_name, status)
                continue

            ftp_relative_path = entry_path.removeprefix("/pub/")
            local_path = f"/nfs/ftp/public/{ftp_relative_path}"
            records.append(
                {
                    "path": local_path,
                    "species": species,
                    "gca_accession": entry_name,
                    "latest_status": status,
                }
            )
            log.info("DELETABLE: %s", local_path)

    log.info("FTP scan complete: %d deletable directories", len(records))
    return records


def write_manifest(records: list[dict[str, str]], output_file: Path) -> None:
    timestamp = datetime.utcnow().isoformat(timespec="seconds") + "Z"
    with output_file.open("w") as handle:
        handle.write(f"# Deletion manifest generated {timestamp}\n")
        handle.write(f"# Total paths: {len(records)}\n")
        handle.write("# Review this file before deletion.\n")
        handle.write("path\tspecies\tgca_accession\tlatest_status\taction\n")
        for record in sorted(records, key=lambda item: item["path"]):
            handle.write(
                "{path}\t{species}\t{gca_accession}\t{latest_status}\tdelete\n".format(
                    **record
                )
            )
    log.info("Manifest written to %s", output_file)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Discover pre-release directories whose latest status is live."
    )
    parser.add_argument("--db-host", default=DEFAULT_DB["host"])
    parser.add_argument("--db-port", type=int, default=DEFAULT_DB["port"])
    parser.add_argument("--db-user", default=DEFAULT_DB["user"])
    parser.add_argument("--db-password", default=DEFAULT_DB["password"])
    parser.add_argument("--db-name", default=DEFAULT_DB["database"])
    parser.add_argument("--ftp-host", default=FTP_HOST)
    parser.add_argument("--ftp-base", default=FTP_BASE_DIR)
    parser.add_argument("--output", type=Path, default=Path("paths_to_delete.tsv"))
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print candidate paths without writing a manifest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    db_config = dict(
        host=args.db_host,
        port=args.db_port,
        user=args.db_user,
        password=args.db_password,
        database=args.db_name,
        cursorclass=pymysql.cursors.DictCursor,
    )

    latest_status_by_gca = get_latest_status_by_gca(db_config)
    ftp = connect_ftp(args.ftp_host)
    try:
        records = find_gca_directories_on_ftp(
            ftp, args.ftp_base, latest_status_by_gca
        )
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    if not records:
        log.info("Nothing to delete")
        return

    if args.dry_run:
        log.info("DRY RUN - candidate paths:")
        for record in sorted(records, key=lambda item: item["path"]):
            print(record["path"])
        return

    write_manifest(records, args.output)


if __name__ == "__main__":
    main()
