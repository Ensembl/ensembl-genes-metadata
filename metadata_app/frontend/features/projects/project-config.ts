export type ProjectConfig = {
  slug: string;
  title: string;
  description: string;
};

export const PROJECTS: ProjectConfig[] = [
  {
    slug: "erga-bge",
    title: "ERGA-BGE*",
    description:
      "European Reference Genome Atlas Biodiversity Genomics Europe Project",
  },
  {
    slug: "erga",
    title: "ERGA*",
    description: "European Reference Genome Atlas Project",
  },
  {
    slug: "erga-pilot",
    title: "ERGA-pilot*",
    description: "European Reference Genome Atlas Pilot Project",
  },
  {
    slug: "cbp",
    title: "CBP*",
    description: "Canadian BioGenome Project",
  },
  {
    slug: "aegis",
    title: "AEGIS*",
    description: "Ancient Environmental Genomics Initiative for Sustainability",
  },
  {
    slug: "ebp",
    title: "EBP*",
    description: "Earth BioGenome Project",
  },
  {
    slug: "vgp",
    title: "VGP*",
    description: "Vertebrate Genomes Project",
  },
  {
    slug: "dtol",
    title: "DToL",
    description: "Darwin Tree of Life Project",
  },
  {
    slug: "tol",
    title: "ToL",
    description: "Tree of Life Project",
  },
  {
    slug: "asg",
    title: "ASG",
    description: "Aquatic Symbiosis Genomics Project",
  },
  {
    slug: "hprc",
    title: "HPRC",
    description: "Human Pangenome Reference Consortium",
  },
  {
    slug: "laca",
    title: "LACA",
    description: "Livestock And Companion Animals",
  },
];

export const PROJECTS_BY_SLUG = Object.fromEntries(
  PROJECTS.map((project) => [project.slug, project]),
) as Record<string, ProjectConfig>;
