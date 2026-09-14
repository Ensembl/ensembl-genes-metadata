from metadata_app.backend.app.services.annotations_service import generate_tables
from metadata_app.backend.app.services.assembly_service import get_filtered_assemblies
from fastapi import HTTPException


def _to_number(value):
    if value is None:
        return None

    try:
        return float(value)
    except (TypeError, ValueError):
        return None


STATUS_INFO = {
    "insufficient_data": {
        "title": "Insufficient data",
        "description": (
            "There is currently not enough transcriptomic evidence to produce "
            "a high-quality annotation."
        ),
        "action": ("Waiting for additional transcriptomic data from the community."),
    },
    "abandoned": {
        "title": "Annotation abandoned",
        "description": (
            "Work on this assembly has stopped in favour of a better-quality "
            "assembly or transcriptomic data has been released since annotation was attempted"
        ),
        "action": "A newer assembly may be annotated instead or we can use evidence based pipeline.",
    },
    "in_progress": {
        "title": "Annotation in progress",
        "description": ("A genebuilder is actively working on this annotation."),
        "action": "Please wait while the annotation pipeline completes.",
    },
    "pre_released": {
        "title": "Pre-release",
        "description": (
            "The annotation has been completed and is awaiting final quality "
            "checks before release."
        ),
        "action": ("The annotation files are available on the pre-release FTP site."),
    },
    "coming_soon": {
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
    "check_busco": {
        "title": "Low protein BUSCO",
        "description": (
            "The predicted proteins have a low BUSCO completeness score, "
            "indicating the annotation may be incomplete."
        ),
        "action": (
            "Additional evidence or future pipeline improvements may improve "
            "the annotation."
        ),
    },
    "poor_genome_busco": {
        "title": "Low genome BUSCO",
        "description": "The genome assembly itself appears incomplete.",
        "action": (
            "A higher-quality assembly would likely produce a better annotation."
        ),
    },
}


def build_reasons(row):
    status = str(row.get("gb_status", "")).lower()

    if status == "check_busco":
        return [
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
        ]

    if status == "poor_genome_busco":
        return [
            {
                "title": "Low genome BUSCO",
                "description": "The genome assembly itself appears incomplete.",
                "action": (
                    "A higher-quality assembly would likely produce a better annotation."
                ),
                "severity": "warning",
            }
        ]

    return []


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


def build_assembly_fallback(gca: str):
    try:
        assembly_result, _, _ = get_filtered_assemblies(
            bioproject_id=None,
            metric_thresholds={},
            asm_level=None,
            asm_type=None,
            release_date=None,
            taxon_id=None,
            current=False,
            transc=False,
            transc_ena=True,
            non_annotated=False,
            group_name=None,
            gca=[gca],
        )
    except HTTPException:
        return None

    if assembly_result is None or getattr(assembly_result, "empty", True):
        return None

    row = assembly_result.iloc[0]
    asm_level = str(row.get("asm_level") or "").lower()
    contig_n50 = _to_number(row.get("contig_n50"))
    rnaseq_lowest = _to_number(row.get("short_read_paired_end_illumina_lowest"))
    rnaseq_genus = _to_number(row.get("short_read_paired_end_illumina"))

    reasons = []

    if asm_level in {"contig", "scaffold"}:
        reasons.append(
            {
                "title": "Assembly is too fragmented",
                "description": "The assembly level is contig or scaffold.",
                "action": "Assembly is too fragmented for good quality annotation.",
                "severity": "warning",
            }
        )
    elif contig_n50 is not None and contig_n50 < 100000:
        reasons.append(
            {
                "title": "Assembly is too fragmented",
                "description": "The contig N50 is below 100000.",
                "action": "Assembly is too fragmented for good quality annotation.",
                "severity": "warning",
            }
        )

    rnaseq_signal = max(rnaseq_lowest or 0, rnaseq_genus or 0)

    if rnaseq_signal > 3:
        reasons.append(
            {
                "title": "Sufficient RNA-seq evidence",
                "description": (
                    "Transcriptomic evidence from ENA is available at the lowest "
                    "or genus taxon level."
                ),
                "action": "We have enough RNA-seq evidence for annotation. We can start annotation.",
                "severity": "info",
            }
        )
    elif 1 <= rnaseq_signal <= 3:
        reasons.append(
            {
                "title": "Limited RNA-seq evidence",
                "description": (
                    "Some transcriptomic evidence is available from ENA, but it is limited."
                ),
                "action": "Annotation may not be good quality with this evidence, but we can try.",
                "severity": "warning",
            }
        )

    return {
        "summary": {
            "gca": row.get("gca"),
            "gb_status": row.get("gb_status"),
            "asm_level": row.get("asm_level"),
            "contig_n50": contig_n50,
            "short_read_paired_end_illumina_lowest": rnaseq_lowest,
            "short_read_paired_end_illumina": rnaseq_genus,
        },
        "status": {
            "title": "Assembly found",
            "description": "No annotation record was found. Assembly metadata was retrieved instead.",
            "action": "Checked ENA transcriptomic data to estimate annotation readiness.",
        },
        "reasons": reasons,
        "timeline": [],
    }


def get_gca_info(gcas: list[str]):
    results = []

    for gca in gcas:
        try:
            anno_wide, _, _ = generate_tables(
                annotation_date=None,
                taxon_id=None,
                bioproject_id=None,
                group_name=None,
                gca=gca,
            )
        except HTTPException:
            fallback = build_assembly_fallback(gca)
            if fallback:
                results.append(fallback)
            continue

        result = anno_wide.loc[anno_wide["gca"] == gca]

        if result.empty:
            fallback = build_assembly_fallback(gca)
            if fallback:
                results.append(fallback)
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
