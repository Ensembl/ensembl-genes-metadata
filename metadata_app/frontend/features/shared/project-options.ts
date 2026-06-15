import type { Option } from "@/components/ui/multi_select";

export const GROUP_NAME_VALUES = ["LACA", "AQUA-FAANG"] as const;

export const PROJECT_OPTIONS: Option[] = [
  { value: "PRJEB40665", label: "Darwin Tree of Life" },
  {
    value: "PRJEB61747",
    label: "European Reference Genome Atlas/Biodiversity Genomics Europe",
  },
  { value: "PRJEB43510", label: "European Reference Genome Atlas" },
  { value: "PRJNA533106", label: "Earth BioGenome" },
  { value: "PRJEB47820", label: "European Reference Genome Atlas pilot" },
  { value: "PRJEB43743", label: "Aquatic Symbiosis" },
  { value: "PRJNA489243", label: "Vertebrate Genomes" },
  {
    value: "PRJEB80366",
    label: "Ancient Environmental Genomics Initiative for Sustainability",
  },
  { value: "PRJNA813333", label: "Canadian BioGenome" },
  { value: "PRJEB43745", label: "Tree of Life" },
  { value: "PRJNA1399476", label: "Rodent 2K" },
  { value: "LACA", label: "Livestock And Companion Animals" },
  { value: "AQUA-FAANG", label: "Aqua FAANG" },
];
