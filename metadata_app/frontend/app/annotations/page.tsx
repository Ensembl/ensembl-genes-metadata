"use client";

import React, { useState } from "react";
import { Loader2 } from "lucide-react";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { DataTable } from "@/app/tables/data-table";
import { Annotations, columns } from "@/app/tables/annotations_columns";
import { Checkbox } from "@/components/ui/checkbox";

type Downloadables = {
  anno_main: string;
  anno_wide: string;
};

export default function Page() {
  const baseFields = [
    { label: "BioProject ID", placeholder: "PRJNA123456" },
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Annotation date", placeholder: "2024-12-31" },
  ];

  const [baseFieldValues, setBaseFieldValues] = useState<{ [key: string]: string }>({});
  const [annotations, setAnnotations] = useState<Annotations[]>([]);
  const [downloadables, setDownloadables] = useState<Downloadables | null>(null);
  const [loading, setLoading] = useState(false);
  const [releaseSites, setReleaseSites] = useState({
    main: true,
    beta: true
  });

  const handleGetAnnotations = async () => {
    setLoading(true);
    try {
      let bioprojectArray = null;
      const bioprojectId = baseFieldValues["BioProject ID"];

      if (bioprojectId) {
        bioprojectArray = bioprojectId.includes(",")
          ? bioprojectId.split(",").map((id) => id.trim()).filter((id) => id)
          : [bioprojectId.trim()];
      }

      const taxon_id = baseFieldValues["Taxon ID"]
        ? parseInt(baseFieldValues["Taxon ID"], 10)
        : null;

     let release_type: string[] | null = [];
    if (releaseSites.main) release_type.push("main");
    if (releaseSites.beta) release_type.push("beta");

    if (release_type.length === 0 || release_type.length === 2) {
      release_type = null; // means "no filter" to the backend
    }

      const payload = {
        bioproject_id: bioprojectArray,
        annotation_date: baseFieldValues["Annotation date"] || null,
        taxon_id: taxon_id,
        release_type: release_type,
      };

      const cleanPayload = Object.fromEntries(
        Object.entries(payload).filter(([_, value]) => value !== null && value !== undefined)
      );

      const res = await fetch("http://127.0.0.1:8000/api/annotations/annotations/filter", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(cleanPayload),
      });

      if (!res.ok) {
        const errorText = await res.text();
        console.error("Failed to fetch annotations", errorText);
        alert("Failed to fetch annotations: " + errorText);
        return;
      }

      const result = await res.json();
      console.log("API response:", result);

      if (result.anno_main) {
        setAnnotations(result.anno_main);
      } else {
        console.warn("No annotation data found in response");
        alert("No annotation data found.");
      }

      if (result.downloadables_anno) {
        setDownloadables(result.downloadables_anno);
      }
    } catch (error) {
      console.error("Error fetching annotations:", error);
      alert("Error fetching annotations: " + (error instanceof Error ? error.message : String(error)));
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = (content: string | undefined, filename: string, type = "text/plain") => {
    if (!content) {
      alert(`No ${filename} data available to download.`);
      return;
    }

    try {
      const blob = new Blob([content], { type });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      URL.revokeObjectURL(url);
    } catch (error) {
      console.error(`Error downloading ${filename}:`, error);
      alert(`Error downloading ${filename}: ${error instanceof Error ? error.message : String(error)}`);
    }
  };

  return (
    <div className="min-h-screen justify-center flex flex-wrap align-items-center">
      <div className="container m-16 mt-10 max-w-6xl">
        <div className="rounded-2xl border-accent">
          <div className="rounded-2xl bg-secondary p-8 gap-10 shadow-lg">
            <div className="grid justify-center grid-cols-4 gap-4">
              {baseFields.map(({ label, placeholder }, index) => (
                <div key={index}>
                  <Label htmlFor={label.toLowerCase().replace(" ", "-")}>
                    {label}
                  </Label>
                  <Input
                    id={label.toLowerCase().replace(" ", "-")}
                    type="text"
                    placeholder={placeholder}
                    className="mt-3 gap-2 bg-filter-input-bg dark:text-background"
                    value={baseFieldValues[label] || ""}
                    onChange={(e) =>
                      setBaseFieldValues((prev) => ({
                        ...prev,
                        [label]: e.target.value,
                      }))
                    }
                  />
                </div>
              ))}

              <div>
                <Label className="mb-2 block">Release site</Label>
                <div className="flex space-x-4 mt-3">
                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="checkbox-main"
                      className= "cursor-pointer"
                      checked={releaseSites.main}
                      onCheckedChange={(checked) => {
                        setReleaseSites(prev => ({
                          ...prev,
                          main: checked === true
                        }));
                      }}
                    />
                    <Label htmlFor="checkbox-main" className="cursor-pointer">Main</Label>
                  </div>
                  <div className="flex items-center space-x-2">
                    <Checkbox
                      id="checkbox-beta"
                      className= "cursor-pointer"
                      checked={releaseSites.beta}
                      onCheckedChange={(checked) => {
                        setReleaseSites(prev => ({
                          ...prev,
                          beta: checked === true
                        }));
                      }}
                    />
                    <Label htmlFor="checkbox-beta" className="cursor-pointer">Beta</Label>
                  </div>
                </div>
              </div>
            </div>
          </div>

          {/* Get Results Button */}
          <div className="mt-8 flex justify-end">
            <Button size="lg" onClick={handleGetAnnotations} disabled={loading}>
              {loading ? (
                <>
                  <Loader2 className="animate-spin mr-2" />
                  Loading...
                </>
              ) : (
                "Get Results"
              )}
            </Button>
          </div>

          {/* Results */}
          {annotations.length > 0 && (
            <div className="mt-10 shadow-lg border border:border rounded-2xl">
              <div className="flex items-center justify-between px-8 py-6 border-b border:border">
                <h2 className="text-lg font-semibold">Filtered Annotations ({annotations.length})</h2>
                <div className="flex gap-2">
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(downloadables?.anno_main, "anno_main.csv", "text/csv")}
                  >
                    Download CSV
                  </Button>
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(downloadables?.anno_wide, "full_table_filtered_annotations.csv", "text/csv")}
                  >
                    Download Full Table
                  </Button>
                </div>
              </div>
              <DataTable columns={columns} data={annotations} />
            </div>
          )}
        </div>
      </div>
    </div>
  );
}