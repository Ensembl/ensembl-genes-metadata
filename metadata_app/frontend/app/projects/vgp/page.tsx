"use client";

import React, { useState, useEffect } from "react";
import { Terminal } from "lucide-react";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { DataTable } from "@/app/tables/data-table-sorting";
import { ProjectGCA, columns } from "@/app/tables/project_columns";

export default function Page() {
  const [assemblies, setAssemblies] = useState<ProjectGCA[]>([]);
  const [loading, setLoading] = useState(true);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);

  useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await fetch("http://127.0.0.1:8000/api/project/project/vgp");

        if (!res.ok) {
          const errorText = await res.text();
          throw new Error(`Failed to fetch: ${errorText}`);
        }

        const json = await res.json();
        console.log("API response:", json);

        const formatted = json.map((item: ProjectGCA) => ({
          gca: item.gca,
          lowest_taxon_id: item.lowest_taxon_id,
          scientific_name: item.scientific_name,
          asm_name: item.asm_name,
          asm_level: item.asm_level,
          gb_status: item.gb_status,
          genebuilder: item.genebuilder,
        }));

        setAssemblies(formatted);
      } catch (err) {
        console.error("Error fetching data:", err);
        setErrorMessage(
          err instanceof Error ? err.message : "Unknown error occurred"
        );
      } finally {
        setLoading(false);
      }
    };

    fetchData();
  }, []);

  return (
    <div className="flex items-center justify-center">
      <div className="container mt-4 m-16 max-w-6xl">
        {errorMessage && (
          <Alert variant="destructive" className="mt-8">
            <Terminal />
            <AlertTitle>Heads up!</AlertTitle>
            <AlertDescription>{errorMessage}</AlertDescription>
          </Alert>
        )}

        <div className="flex items-center justify-between px-8 py-6">
          <h1 className="text-2xl font-semibold">
            Darwin Tree of Life Number of GCAs ({assemblies.length})
          </h1>
        </div>

        {loading ? (
          <p className="text-muted-foreground">Loading data...</p>
        ) : (
            <div className="border-border border-2 rounded-md shadow-border">
            <DataTable columns={columns} data={assemblies} />
            </div>
        )}
      </div>
    </div>
  );
}