import { ProjectPageClient } from "@/features/projects/project-page-client";
import { PROJECTS, PROJECTS_BY_SLUG } from "@/features/projects/project-config";

export function generateStaticParams() {
  return PROJECTS.map((project) => ({
    project: project.slug,
  }));
}

type ProjectPageProps = {
  params: Promise<{
    project: string;
  }>;
};

export default async function ProjectPage({ params }: ProjectPageProps) {
  const { project: projectSlug } = await params;
  const project = PROJECTS_BY_SLUG[projectSlug];

  return (
    <ProjectPageClient
      projectSlug={projectSlug}
      projectTitle={project?.title ?? projectSlug}
      projectExists={Boolean(project)}
    />
  );
}
