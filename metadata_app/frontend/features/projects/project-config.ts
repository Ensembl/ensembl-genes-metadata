export type ProjectConfig = {
  slug: string;
  title: string;
  description: string;
  reportKey: string;
  reportLabel?: string;
  reportGroup: "biodiversity" | "custom";
};

export const PROJECTS: ProjectConfig[] = [
  {
    slug: "erga-bge",
    title: "ERGA-BGE*",
    description:
      "European Reference Genome Atlas Biodiversity Genomics Europe Project",
    reportKey: "ERGA/BGE",
    reportGroup: "biodiversity",
  },
  {
    slug: "erga",
    title: "ERGA*",
    description: "European Reference Genome Atlas Project",
    reportKey: "ERGA",
    reportGroup: "biodiversity",
  },
  {
    slug: "erga-pilot",
    title: "ERGA-pilot*",
    description: "European Reference Genome Atlas Pilot Project",
    reportKey: "ERGA_pilot",
    reportLabel: "ERGA Pilot",
    reportGroup: "biodiversity",
  },
  {
    slug: "cbp",
    title: "CBP*",
    description: "Canadian BioGenome Project",
    reportKey: "CBP",
    reportGroup: "biodiversity",
  },
  {
    slug: "aegis",
    title: "AEGIS*",
    description: "Ancient Environmental Genomics Initiative for Sustainability",
    reportKey: "AEGIS",
    reportGroup: "biodiversity",
  },
  {
    slug: "ebp",
    title: "EBP*",
    description: "Earth BioGenome Project",
    reportKey: "EBP",
    reportGroup: "biodiversity",
  },
  {
    slug: "vgp",
    title: "VGP*",
    description: "Vertebrate Genomes Project",
    reportKey: "VGP",
    reportGroup: "biodiversity",
  },
  {
    slug: "dtol",
    title: "DToL",
    description: "Darwin Tree of Life Project",
    reportKey: "DToL",
    reportGroup: "biodiversity",
  },
  {
    slug: "tol",
    title: "ToL",
    description: "Tree of Life Project",
    reportKey: "ToL",
    reportGroup: "biodiversity",
  },
  {
    slug: "asg",
    title: "ASG",
    description: "Aquatic Symbiosis Genomics Project",
    reportKey: "ASG",
    reportGroup: "biodiversity",
  },
    {
    slug: "atlasea",
    title: "ATLASea",
    description: "ATLASea Project",
    reportKey: "ATLASea",
    reportGroup: "biodiversity",
  },
  {
    slug: "hprc",
    title: "HPRC",
    description: "Human Pangenome Reference Consortium",
    reportKey: "HPRC",
    reportGroup: "custom",
  },
  {
    slug: "laca",
    title: "LACA",
    description: "Livestock And Companion Animals",
    reportKey: "LACA",
    reportGroup: "custom",
  },
    {
    slug: "rodent2k",
    title: "Rodent2K",
    description: "Rodent 2K",
    reportKey: "Rodent2K",
    reportGroup: "biodiversity",
  },
];

export const PROJECTS_BY_SLUG = Object.fromEntries(
  PROJECTS.map((project) => [project.slug, project]),
) as Record<string, ProjectConfig>;
