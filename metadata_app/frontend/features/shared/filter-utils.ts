import { GROUP_NAME_VALUES } from "@/features/shared/project-options";

type ProjectSelection = {
  value: string;
};

export function splitCommaSeparated(value?: string | null): string[] {
  if (!value) {
    return [];
  }

  return value
    .split(",")
    .map((item) => item.trim())
    .filter(Boolean);
}

export function parseTaxonIds(value?: string | null): number[] | null {
  const taxonIds = splitCommaSeparated(value)
    .filter((id) => /^\d+$/.test(id))
    .map(Number);

  return taxonIds.length > 0 ? taxonIds : null;
}

export function splitProjectFilters(
  selectedProjects: ProjectSelection[],
  manualBioprojectIds?: string | null,
  groupNameValues: readonly string[] = GROUP_NAME_VALUES,
) {
  const bioprojectIds: string[] = [];
  const groupNames: string[] = [];

  selectedProjects.forEach((project) => {
    if (groupNameValues.includes(project.value)) {
      groupNames.push(project.value);
    } else {
      bioprojectIds.push(project.value);
    }
  });

  bioprojectIds.push(...splitCommaSeparated(manualBioprojectIds));

  return {
    bioprojectIds: Array.from(new Set(bioprojectIds)),
    groupNames,
  };
}

export function buildScopePayload(
  selectedProjects: ProjectSelection[],
  manualBioprojectIds?: string | null,
  taxonValue?: string | null,
) {
  const { bioprojectIds, groupNames } = splitProjectFilters(
    selectedProjects,
    manualBioprojectIds,
  );

  return {
    bioproject_id: bioprojectIds.length > 0 ? bioprojectIds : null,
    group_name: groupNames.length > 0 ? groupNames : null,
    taxon_id: parseTaxonIds(taxonValue),
  };
}

export function cleanPayload<T extends Record<string, unknown>>(payload: T) {
  return Object.fromEntries(
    Object.entries(payload).filter(([, value]) => {
      if (value === null || value === undefined) {
        return false;
      }

      return !(Array.isArray(value) && value.length === 0);
    }),
  ) as Partial<T>;
}
