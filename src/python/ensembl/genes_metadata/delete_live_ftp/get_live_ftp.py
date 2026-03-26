"""
Connects to registry DB to get live GCA accessions, checks FTP for matching directories,
and writes a manifest of paths that are safe to delete.

Run as: yourself (no special permissions needed)
"""

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


# Configuration — override via CLI args or environment
DEFAULT_DB = dict(
    host="mysql-ens-genebuild-prod-1",
    port=4527,
    user="ensro",
    password="",  # prefer env var MYSQL_PWD or ~/.my.cnf
    database="gb_assembly_metadata",
    cursorclass=pymysql.cursors.DictCursor,
)

FTP_HOST = "ftp.ebi.ac.uk"
FTP_BASE_DIR = "/pub/databases/ensembl/pre-release"


# DB helpers
def get_live_gca_accessions(db_config: dict) -> set[str]:
    """Return the set of GCA accessions currently live in the database."""
    log.info("Connecting to database %s@%s …", db_config["database"], db_config["host"])
    conn = pymysql.connect(**db_config)
    try:
        with conn.cursor() as cur:
            # ----------------------------------------------------------------
            # !! Adjust table/column names to match your schema !!
            # ----------------------------------------------------------------
            sql = """
                SELECT DISTINCT gca_accession
                FROM genebuild_status       
                WHERE  gb_status = 'live'
            """
            cur.execute(sql)
            rows = cur.fetchall()
    finally:
        conn.close()

    accessions = {row["gca_accession"].strip() for row in rows}
    log.info("Found %d live GCA accessions in DB.", len(accessions))
    return accessions


# FTP helpers
def connect_ftp(host: str) -> ftplib.FTP:
    log.info("Connecting to FTP %s …", host)
    ftp = ftplib.FTP(host, timeout=30)
    # Public FTP: anonymous login (no credentials needed)
    ftp.login()
    ftp.set_pasv(True)
    log.info("FTP anonymous login successful.")
    return ftp


def list_ftp_directory(ftp: ftplib.FTP, path: str) -> list[str]:
    """Return names of entries directly under *path*, or [] if it doesn't exist."""
    try:
        return ftp.nlst(path)
    except ftplib.error_perm as exc:
        if "550" in str(exc):  # No such file or directory
            return []
        raise


def find_gca_directories_on_ftp(
    ftp: ftplib.FTP,
    base_dir: str,
    live_accessions: set[str],
) -> list[str]:
    """
    Walk the FTP tree under *base_dir* looking one level deep for GCA directories.

    Structure assumed:
        base_dir/<species>/GCA_xxx...
    """
    deletable: list[str] = []

    log.info("Scanning FTP base dir (1 level deep): %s", base_dir)

    level1_dirs = list_ftp_directory(ftp, base_dir)

    for lvl1 in level1_dirs:
        # lvl1 is species directory, e.g. Ajuga_reptans
        gca_dirs = list_ftp_directory(ftp, lvl1)

        for entry_path in gca_dirs:
            entry_name = entry_path.split("/")[-1]

            if not entry_name.upper().startswith("GCA_"):
                continue

            parts = entry_name.split("_")
            if len(parts) >= 2:
                accession = f"{parts[0]}_{parts[1]}"
            else:
                accession = entry_name

            if accession in live_accessions:
                full_path = f"/nfs/ftp/public{entry_path}"
                log.info(
                    "  DELETABLE: %s  (accession %s live)",
                    full_path,
                    accession,
                )
                deletable.append(full_path)
            else:
                log.debug("  KEEP: %s", entry_name)

    log.info(
        "FTP scan complete. %d deletable director%s found.",
        len(deletable),
        "y" if len(deletable) == 1 else "ies",
    )
    return deletable


# Manifest writer
def write_manifest(paths: list[str], output_file: Path) -> None:
    timestamp = datetime.utcnow().isoformat(timespec="seconds") + "Z"
    with output_file.open("w") as fh:
        fh.write(f"# Deletion manifest generated {timestamp}\n")
        fh.write(f"# Total paths: {len(paths)}\n")
        fh.write("# Hand this file to the deletion script (phase 2).\n\n")
        for p in sorted(paths):
            fh.write(p + "\n")
    log.info("Manifest written to: %s", output_file)


# CLI
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Phase 1 — discover FTP dirs with live GCA accessions."
    )
    p.add_argument("--db-host", default=DEFAULT_DB["host"])
    p.add_argument("--db-port", type=int, default=DEFAULT_DB["port"])
    p.add_argument("--db-user", default=DEFAULT_DB["user"])
    p.add_argument("--db-password", default=DEFAULT_DB["password"])
    p.add_argument("--db-name", default=DEFAULT_DB["database"])
    p.add_argument("--ftp-host", default=FTP_HOST)
    p.add_argument("--ftp-base", default=FTP_BASE_DIR)
    p.add_argument(
        "--output",
        default="paths_to_delete.txt",
        help="Path to write the deletion manifest (default: paths_to_delete.txt)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print paths but do not write the manifest.",
    )
    return p.parse_args()


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

    live_accessions = get_live_gca_accessions(db_config)

    ftp = connect_ftp(args.ftp_host)
    try:
        deletable_paths = find_gca_directories_on_ftp(
            ftp, args.ftp_base, live_accessions
        )
    finally:
        try:
            ftp.quit()
        except Exception:
            pass

    if not deletable_paths:
        log.info("Nothing to delete — exiting.")
        return

    if args.dry_run:
        log.info("DRY RUN — paths that would be written to manifest:")
        for p in sorted(deletable_paths):
            print(p)
        return

    write_manifest(deletable_paths, Path(args.output))
    print()
    print("=" * 60)
    print("Next step — hand the manifest to the deletion script:")
    print()
    print("  1. Copy the manifest to a shared location accessible by genebuild.")
    print("  2. Switch to the genebuild user")
    print("  3. Request a datamover node")
    print("  4. Run:  python delete_live_from_ftp.py --manifest paths_to_delete.txt")
    print("=" * 60)


if __name__ == "__main__":
    main()
