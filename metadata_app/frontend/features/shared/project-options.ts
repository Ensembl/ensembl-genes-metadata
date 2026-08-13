import type { Option } from "@/components/ui/multi_select";

export const GROUP_NAME_VALUES = ["LACA", "AQUA-FAANG"] as const;

export type ProjectOption = Option & {
  slug?: string;
};

export const PROJECT_OPTIONS: ProjectOption[] = [
  { value: "PRJEB40665", label: "Darwin Tree of Life", slug: "dtol" },
  {
    value: "PRJEB61747",
    label: "European Reference Genome Atlas/Biodiversity Genomics Europe",
    slug: "erga-bge",
  },
  { value: "PRJEB43510", label: "European Reference Genome Atlas", slug: "erga" },
  { value: "PRJNA533106", label: "Earth BioGenome", slug: "ebp" },
  { value: "PRJEB47820", label: "European Reference Genome Atlas pilot", slug: "erga-pilot" },
  { value: "PRJEB43743", label: "Aquatic Symbiosis", slug: "asg" },
  { value: "PRJNA489243", label: "Vertebrate Genomes", slug: "vgp" },
  {
    value: "PRJEB80366",
    label: "Ancient Environmental Genomics Initiative for Sustainability",
    slug: "aegis",
  },
  { value: "PRJNA813333", label: "Canadian BioGenome", slug: "cbp" },
  { value: "PRJEB43745", label: "Tree of Life", slug: "tol" },
  { value: "PRJNA1399476", label: "Rodent 2K", slug: "rodent2k" },
  { value: "LACA", label: "Livestock And Companion Animals", slug: "laca" },
  { value: "AQUA-FAANG", label: "Aqua FAANG" },
  { value: "PRJEB64126", label: "ATLASea", slug: "atlasea" },
];

export const PROJECT_OPTION_BY_SLUG = Object.fromEntries(
  PROJECT_OPTIONS.filter((option) => option.slug).map((option) => [
    option.slug,
    { value: option.value, label: option.label },
  ]),
) as Record<string, Option>;
