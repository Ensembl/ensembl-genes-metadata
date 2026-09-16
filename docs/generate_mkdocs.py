#!/usr/bin/env python3
"""
Generate MkDocs Material documentation for the Ensembl Genes Nextflow pipelines.

Mirrors the layout used by ensembl-genes-nf's live documentation site
(mkdocs.yml built with mkdocs/zensical), reusing the same Nextflow
module/workflow/schema parsing as the Sphinx generator (see
docs/generate_docs.py). The content lives under docs_mkdocs/, a sibling of
docs/, so it never overlaps with the Sphinx source tree or this package.

Only generated pages are overwritten. Handwritten pages (docs_mkdocs/README.md
and each pipeline's index.md, input.md, output.md, troubleshooting.md) are
created once and then left alone.

Run before the MkDocs build:

    python -m docs.generate_mkdocs
"""

from __future__ import annotations

import logging
from pathlib import Path

from docs.generators.module_page import render_module
from docs.generators.parameters import load_schema
from docs.generators.parameters import render as render_parameters
from docs.generators.parser import parse_module
from docs.generators.parser import parse_workflow
from docs.generators.utils import write_file
from docs.generators.workflow_page import render_workflow

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
)

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[1]

PIPELINES_DIR = REPO_ROOT / "pipelines"

DOCS_DIR = REPO_ROOT / "docs_mkdocs"


def _title(name: str) -> str:
    return name.replace("_", " ").title()


def _ensure_page(path: Path, content: str) -> None:
    """
    Create a handwritten page only if it does not already exist.
    """

    if path.exists():
        return

    logger.info("Creating %s", path.relative_to(REPO_ROOT))

    path.parent.mkdir(parents=True, exist_ok=True)

    path.write_text(content, encoding="utf8")


def _generate_pipeline(pipeline_dir: Path) -> str | None:
    """
    Generate the MkDocs pages for one pipeline. Returns the pipeline name,
    or None if the directory does not look like a Nextflow pipeline.
    """

    name = pipeline_dir.name

    out_dir = DOCS_DIR / "pipelines" / name

    modules = []

    modules_dir = pipeline_dir / "modules"

    if modules_dir.exists():

        for nf in sorted(modules_dir.glob("*.nf")):

            module = parse_module(nf)

            module.source = nf.relative_to(REPO_ROOT)

            modules.append(module)

        modules.sort(key=lambda m: m.process)

    workflows = []

    main_nf = pipeline_dir / "main.nf"

    if main_nf.exists():

        workflow = parse_workflow(main_nf)

        workflow.source = main_nf.relative_to(REPO_ROOT)

        workflows.append(workflow)

    if not modules and not workflows:
        return None

    #
    # Handwritten skeleton pages.
    #

    _ensure_page(
        out_dir / "index.md",
        f"# {_title(name)}\n\nOverview of the **{_title(name)}** pipeline.\n",
    )

    _ensure_page(out_dir / "input.md", "# Input\n")

    _ensure_page(out_dir / "output.md", "# Output\n")

    _ensure_page(out_dir / "troubleshooting.md", "# Troubleshooting\n")

    #
    # Parameters.
    #

    schema_path = pipeline_dir / "nextflow_schema.json"

    if schema_path.exists():

        write_file(
            out_dir / "parameters.md",
            render_parameters(load_schema(schema_path)),
        )

    #
    # Modules.
    #

    if modules:

        modules_out = out_dir / "modules"

        for module in modules:

            write_file(
                modules_out / f"{module.slug}.md",
                render_module(module),
            )

        lines = [
            "# Modules",
            "",
            "Documentation for the modules used by this pipeline.",
            "",
        ]

        for module in modules:
            lines.append(f"- [{module.process}]({module.slug}.md)")

        write_file(modules_out / "index.md", "\n".join(lines))

    #
    # Workflows.
    #

    if workflows:

        workflows_out = out_dir / "workflows"

        for workflow in workflows:

            write_file(
                workflows_out / f"{workflow.name}.md",
                render_workflow(workflow),
            )

        lines = [
            "# Workflows",
            "",
            "Documentation for the workflows used by this pipeline.",
            "",
        ]

        for workflow in workflows:
            lines.append(f"- [{workflow.name.replace('_', ' ').title()}]({workflow.name}.md)")

        write_file(workflows_out / "index.md", "\n".join(lines))

    return name


def _homepage_template(names: list[str]) -> str:

    lines = [
        "# Ensembl Genes Metadata",
        "",
        "This documentation contains all Ensembl Genes Metadata pipelines.",
        "",
        "## Available Pipelines",
        "",
    ]

    for name in names:

        lines.append(f"### {_title(name)}")
        lines.append("")
        lines.append(f"- [{_title(name)}](pipelines/{name}/index.md)")
        lines.append("")

    return "\n".join(lines)


def main() -> None:
    """
    Generate the MkDocs documentation for every pipeline.
    """

    logger.info("Discovering pipelines...")

    names = []

    for pipeline_dir in sorted(PIPELINES_DIR.iterdir()):

        if not pipeline_dir.is_dir():
            continue

        name = _generate_pipeline(pipeline_dir)

        if name:
            names.append(name)

    logger.info(
        "Generated MkDocs pages for %d pipeline(s).",
        len(names),
    )

    _ensure_page(
        DOCS_DIR / "README.md",
        _homepage_template(names),
    )

    logger.info("MkDocs documentation completed successfully.")


if __name__ == "__main__":
    main()
