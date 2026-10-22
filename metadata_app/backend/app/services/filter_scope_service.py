import logging

from fastapi import HTTPException

from metadata_app.backend.app.services.taxonomy_service import get_descendant_taxa


def as_list(value):
    if value is None:
        return []
    return [value] if isinstance(value, str) else list(value)


def expand_taxon_ids(taxon_ids):
    """Return one deterministic descendant-taxon set for all requested taxa."""
    requested_taxa = sorted({int(taxon_id) for taxon_id in as_list(taxon_ids)})
    expanded_taxa = set()

    for taxon_id in requested_taxa:
        descendants = get_descendant_taxa(taxon_id)
        if not descendants:
            logging.warning("No descendants found for taxon ID %s", taxon_id)
        expanded_taxa.update(descendants)

    if requested_taxa and not expanded_taxa:
        raise HTTPException(
            status_code=400,
            detail=(
                "No descendant taxa found for any of the provided Taxon IDs: "
                + ", ".join(map(str, requested_taxa))
            ),
        )

    return sorted(expanded_taxa)


def add_scope_conditions(
    conditions,
    parameters,
    *,
    bioproject_id=None,
    group_name=None,
    taxon_id=None,
    project_alias="b",
    group_alias="g",
    taxon_column="s.lowest_taxon_id",
):
    """Append the shared project/group/taxon SQL scope to a query."""
    bioproject_ids = as_list(bioproject_id)
    group_names = as_list(group_name)
    project_conditions = []

    if bioproject_ids:
        project_conditions.append(
            f"{project_alias}.bioproject_id IN " f"({','.join(['%s'] * len(bioproject_ids))})"
        )
        parameters.extend(bioproject_ids)

    if group_names:
        project_conditions.append(f"{group_alias}.group_name IN " f"({','.join(['%s'] * len(group_names))})")
        parameters.extend(group_names)

    if project_conditions:
        conditions.append(f"({' OR '.join(project_conditions)})")

    expanded_taxa = expand_taxon_ids(taxon_id)
    if expanded_taxa:
        conditions.append(f"{taxon_column} IN ({','.join(['%s'] * len(expanded_taxa))})")
        parameters.extend(expanded_taxa)

    return conditions, parameters
