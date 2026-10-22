#!/usr/bin/env bash
# Render the nf-metro pipeline diagrams from their sources, pipelines/<pipeline>/assets/metro_map.mmd:
#   pipelines/<pipeline>/<pipeline>_animated.svg  (shown in the pipeline's README)
#   docs/img/<pipeline>.svg                       (shown on the docs site)
#
# Needs nf-metro at the pinned version on PATH (Python 3.11+):
#   pip install "nf-metro==2.1.0"
# Run from anywhere; edit the .mmd, then re-run this and commit the .mmd and both SVGs together.
#
# Usage: src/bash/render_metro_maps.sh [pipeline ...]   (default: every pipeline with a map)

set -euo pipefail

NF_METRO_VERSION="2.1.0"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if ! command -v nf-metro >/dev/null 2>&1; then
    echo "nf-metro not found; install it with: pip install \"nf-metro==${NF_METRO_VERSION}\"" >&2
    exit 1
fi
installed="$(nf-metro --version | awk '{print $NF}')"
if [[ "${installed}" != "${NF_METRO_VERSION}" ]]; then
    echo "nf-metro ${installed} found, but the diagrams are rendered with ${NF_METRO_VERSION}" >&2
    echo "(a different version can lay them out differently): pip install \"nf-metro==${NF_METRO_VERSION}\"" >&2
    exit 1
fi

if [[ $# -gt 0 ]]; then
    pipelines=("$@")
else
    pipelines=()
    for map in "${REPO_ROOT}"/pipelines/*/assets/metro_map.mmd; do
        pipelines+=("$(basename "$(dirname "$(dirname "${map}")")")")
    done
fi

for pipeline in "${pipelines[@]}"; do
    map="${REPO_ROOT}/pipelines/${pipeline}/assets/metro_map.mmd"
    svg="${REPO_ROOT}/pipelines/${pipeline}/${pipeline}_animated.svg"
    docs_svg="${REPO_ROOT}/docs/img/${pipeline}.svg"
    if [[ ! -f "${map}" ]]; then
        echo "No map for ${pipeline}: ${map}" >&2
        exit 1
    fi
    nf-metro validate "${map}"
    nf-metro render "${map}" --animate -o "${svg}"
    cp "${svg}" "${docs_svg}"
    echo "Rendered ${pipeline}: ${svg#"${REPO_ROOT}"/} and ${docs_svg#"${REPO_ROOT}"/}"
done
