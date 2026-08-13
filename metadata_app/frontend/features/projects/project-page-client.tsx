"use client";

import React, { useEffect, useState } from "react";
import { Terminal } from "lucide-react";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { columns, type ProjectGCA } from "@/features/projects/columns";
import { ProjectDataTable } from "@/features/projects/project-data-table";

type ProjectPageClientProps = {
  projectSlug: string;
  projectTitle: string;
  projectExists: boolean;
};

export function ProjectPageClient({
  projectSlug,
  projectTitle,
  projectExists,
}: ProjectPageClientProps) {
  const [assemblies, setAssemblies] = useState<ProjectGCA[]>([]);
  const [loading, setLoading] = useState(true);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);

  useEffect(() => {
    const fetchData = async () => {
      if (!projectExists) {
        setErrorMessage(`Unknown project: ${projectSlug}`);
        setLoading(false);
        return;
      }

      try {
        const res = await fetch(`/api/project/project/${projectSlug}`);

        if (!res.ok) {
          const errorText = await res.text();
          throw new Error(`Failed to fetch: ${errorText}`);
        }

        const json = await res.json();
        const formatted = json.map((item: ProjectGCA) => ({
          gca: item.gca,
          lowest_taxon_id: item.lowest_taxon_id,
          scientific_name: item.scientific_name,
          asm_level: item.asm_level,
          infra_name: item.infra_name,
          gb_status: item.gb_status,
          genebuilder: item.genebuilder,
        }));

        setAssemblies(formatted);
      } catch (err) {
        console.error("Error fetching data:", err);
        setErrorMessage(
          err instanceof Error ? err.message : "Unknown error occurred",
        );
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, [projectExists, projectSlug]);

  return (
    <div className="flex justify-center px-4 py-10 sm:px-6 lg:px-8">
      <div className="w-full max-w-6xl">
        {errorMessage && (
          <Alert variant="destructive" className="mt-8">
            <Terminal />
            <AlertTitle>Heads up!</AlertTitle>
            <AlertDescription>{errorMessage}</AlertDescription>
          </Alert>
        )}

        <div className="flex items-center justify-between py-6 sm:px-8">
          <h1 className="text-xl font-semibold sm:text-2xl">
            {projectTitle} number of GCAs: {assemblies.length}
          </h1>
        </div>

        {loading ? (
          <p className="text-muted-foreground">Loading data...</p>
        ) : (
          <div>
            <ProjectDataTable
              columns={columns}
              data={assemblies}
              projectSlug={projectSlug}
            />
          </div>
        )}
      </div>
    </div>
  );
}
