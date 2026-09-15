"""Generate Sphinx parameter-reference pages from each pipeline's nextflow_schema.json.

Run before the Sphinx build (see .github/workflows/docs-ci.yml):

    python -m docs.generate_docs
"""

import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
PIPELINES_DIR = REPO_ROOT / "pipelines"
OUTPUT_DIR = REPO_ROOT / "docs" / "source" / "pipelines"


def _format_default(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return f"`{str(value).lower()}`"
    return f"`{value}`"


def _ordered_groups(schema: dict) -> list[str]:
    defs = schema.get("$defs", {})
    ref_order = [ref["$ref"].rsplit("/", 1)[-1] for ref in schema.get("allOf", []) if "$ref" in ref]
    ordered = [name for name in ref_order if name in defs]
    ordered += [name for name in defs if name not in ordered]
    return ordered


def _render_group(group: dict) -> list[str]:
    required = group.get("required", [])
    lines = [f"## {group.get('title', '')}", ""]
    description = group.get("description")
    if description:
        lines += [description, ""]
    lines += ["| Parameter | Type | Default | Description |", "|---|---|---|---|"]
    for param, spec in group.get("properties", {}).items():
        marker = " (required)" if param in required else ""
        lines.append(
            f"| `--{param}`{marker} | {spec.get('type', '')} | "
            f"{_format_default(spec.get('default'))} | {spec.get('description', '')} |"
        )
    lines.append("")
    return lines


def render_schema(schema_path: Path) -> str:
    schema = json.loads(schema_path.read_text())
    lines = [f"# {schema.get('title', schema_path.parent.name)}", ""]
    description = schema.get("description")
    if description:
        lines += [description, ""]
    defs = schema.get("$defs", {})
    for name in _ordered_groups(schema):
        lines += _render_group(defs[name])
    return "\n".join(lines)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for schema_path in sorted(PIPELINES_DIR.glob("*/nextflow_schema.json")):
        pipeline_name = schema_path.parent.name
        out_path = OUTPUT_DIR / f"{pipeline_name}.md"
        out_path.write_text(render_schema(schema_path) + "\n")
        print(f"Wrote {out_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
