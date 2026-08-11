from metadata_app.backend.app.services.annotations_service import generate_tables
import pandas as pd

STATUS_INFO = {
    "insufficient data": {
        "title": "Insufficient data",
        "description": (
            "There is currently not enough transcriptomic evidence to produce "
            "a high-quality annotation."
        ),
        "action": (
            "Waiting for additional transcriptomic data from the research " "community."
        ),
    },
    "abandoned": {
        "title": "Annotation abandoned",
        "description": (
            "Work on this assembly has stopped in favour of a better-quality "
            "assembly."
        ),
        "action": "A newer assembly is being annotated instead.",
    },
    "in progress": {
        "title": "Annotation in progress",
        "description": ("A genebuilder is actively working on this annotation."),
        "action": "Please wait while the annotation pipeline completes.",
    },
    "pre-released": {
        "title": "Pre-release",
        "description": (
            "The annotation has been completed and is awaiting final quality "
            "checks before release."
        ),
        "action": ("The annotation files are available on the pre-release FTP site."),
    },
    "coming soon": {
        "title": "Coming soon",
        "description": (
            "The annotation has passed quality control and is queued for the next Ensembl release."
        ),
        "action": "Wait for the next public release.",
    },
    "live": {
        "title": "Live",
        "description": ("The annotation has been officially released."),
        "action": "The annotation is available for download and browsing.",
    },
}


def build_reasons(row):
    reasons = []

    protein_busco = pd.to_numeric(
        row.get("protein_busco"),
        errors="coerce",
    )

    assembly_busco = pd.to_numeric(
        row.get("assembly_busco"),
        errors="coerce",
    )

    if pd.notna(protein_busco) and protein_busco < 90:
        reasons.append(
            {
                "title": "Low protein BUSCO",
                "description": (
                    "The predicted proteins have a low BUSCO completeness score, "
                    "indicating the annotation may be incomplete."
                ),
                "action": (
                    "Additional evidence or future pipeline improvements may "
                    "improve the annotation."
                ),
                "severity": "warning",
            }
        )

    if pd.notna(assembly_busco) and assembly_busco < 90:
        reasons.append(
            {
                "title": "Low genome BUSCO",
                "description": "The genome assembly itself appears incomplete.",
                "action": (
                    "A higher-quality assembly would likely produce a better annotation."
                ),
                "severity": "warning",
            }
        )

    return reasons


def build_timeline(row):
    timeline = []

    if row["date_started"]:
        timeline.append(
            {
                "date": row["date_started"],
                "event": "Annotation started",
            }
        )

    if row["last_genebuild_update"]:
        timeline.append(
            {
                "date": row["last_genebuild_update"],
                "event": "Last annotation update",
            }
        )

    if row["date_status_update"]:
        timeline.append(
            {
                "date": row["date_status_update"],
                "event": f"Status changed to {row['gb_status']}",
            }
        )

    if row["release_date"]:
        timeline.append(
            {
                "date": row["release_date"],
                "event": "Released",
            }
        )

    return sorted(timeline, key=lambda x: x["date"])


def get_gca_info(gcas: list[str]):
    results = []

    for gca in gcas:
        anno_wide, _, _ = generate_tables(
            annotation_date=None,
            taxon_id=None,
            bioproject_id=None,
            group_name=None,
            gca=gca,
        )

        result = anno_wide.loc[anno_wide["gca"] == gca]

        if result.empty:
            continue

        row = result.iloc[0]

        status = str(row["gb_status"]).lower()

        results.append(
            {
                "summary": row.to_dict(),
                "status": STATUS_INFO.get(
                    status,
                    {
                        "title": row["gb_status"],
                        "description": "",
                        "action": "",
                    },
                ),
                "reasons": build_reasons(row),
                "timeline": build_timeline(row),
            }
        )

    return results
