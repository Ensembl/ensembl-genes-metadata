"use client";

import React, { useState } from "react";
import {InfoIcon, Loader2, Terminal} from "lucide-react";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { DataTable } from "@/components/tables/data-table";
import { Annotations, columns } from "@/features/annotations/columns";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import MultipleSelector, { Option } from "@/components/ui/multi_select";
import {InputGroup, InputGroupAddon, InputGroupButton, InputGroupInput} from "@/components/ui/input-group";
import {Tooltip, TooltipContent, TooltipTrigger} from "@/components/ui/tooltip";
import { cleanPayload, parseTaxonIds, splitCommaSeparated, splitProjectFilters } from "@/features/shared/filter-utils";
import { PROJECT_OPTIONS } from "@/features/shared/project-options";


type DownloadableFile = {
  filename: string;
  csv: string;
};

type Downloadables = {
  anno_main: DownloadableFile;
  anno_wide: DownloadableFile;
};

export default function Page() {
  const baseFields = [
      { label: "BioProject ID", placeholder: "PRJNA123456" },
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Annotation date", placeholder: "2024-12-31" },
  ];

  const [baseFieldValues, setBaseFieldValues] = useState<{ [key: string]: string }>({});
  const [selectedProjects, setSelectedProjects] = useState<Option[]>([]);
  const [annotations, setAnnotations] = useState<Annotations[]>([]);
  const [downloadables, setDownloadables] = useState<Downloadables | null>(null);
  const [loading, setLoading] = useState(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);
  const [gcaInput, setGcaInput] = useState<string>("");

  const handleGetAnnotations = async (): Promise<void> => {
    setErrorMessage(null);
    setAnnotations([]);
    setLoading(true);

    try {
      const { bioprojectIds, groupNames } = splitProjectFilters(
        selectedProjects,
        baseFieldValues["BioProject ID"],
      );
      const taxonIdArray = parseTaxonIds(baseFieldValues["Taxon ID"]);
      const uniqueGCA = splitCommaSeparated(gcaInput);

      const payload = {
        bioproject_id: bioprojectIds.length > 0 ? bioprojectIds : null,
        group_name: groupNames.length > 0 ? groupNames : null,
        annotation_date: baseFieldValues["Annotation date"] || null,
        taxon_id: taxonIdArray,
        gca: uniqueGCA.length > 0 ? uniqueGCA : null,
      };

      const cleanedPayload = cleanPayload(payload);

      const res = await fetch("/api/annotations/annotations/filter", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(cleanedPayload),
      });

      if (!res.ok) {
        const errorText = await res.text();
        console.error("Failed to fetch annotations", errorText);
        setErrorMessage("Failed to fetch report: " + errorText);
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
      setErrorMessage("Error fetching data: " + (error instanceof Error ? error.message : String(error)));
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = (file: DownloadableFile | undefined, type = "text/plain"): void => {
  if (!file) {
    alert("No file available to download.");
    return;
  }

  try {
    const blob = new Blob([file.csv], { type });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.href = url;
    link.download = file.filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  } catch (error) {
    console.error(`Error downloading ${file?.filename}:`, error);
    alert(
      `Error downloading ${file?.filename}: ${
        error instanceof Error ? error.message : String(error)
      }`
    );
  }
};

  return (
    <div className="min-h-screen flex justify-center px-4 py-10 sm:px-6 lg:px-8">
      <div className="w-full max-w-6xl">
        <div className="rounded-2xl border-accent">
          <div className="rounded-2xl bg-secondary p-4 gap-10 shadow-lg sm:p-8">
            <div className="mb-6">
              <h1 className="text-2xl font-bold">Set filters</h1>
            </div>
            <div className="grid grid-cols-1 gap-4 mb-4 md:grid-cols-2">
            <div>
                <Label className="mb-3 block">Select project name</Label>
                <MultipleSelector
                  placeholder="Select projects or groups..."
                  defaultOptions={PROJECT_OPTIONS}
                  onChange={(values) => setSelectedProjects(values)}
                />
              </div>
              <div>
                <Label className="mb-3 block">GCA(s)</Label>
                <InputGroup className="mt-3 gap-2 bg-filter-input-bg dark:bg-transparent">
                <InputGroupInput
                    placeholder="GCA_12345333.1, GCA_23456782.2"
                    value={gcaInput}
                    onChange={(e) => setGcaInput(e.target.value)}/>
                <InputGroupAddon align="inline-end">
                  <Tooltip>
                    <TooltipTrigger asChild>
                      <InputGroupButton
                        variant="ghost"
                        aria-label="Info"
                        size="icon-xs"
                      >
                        <InfoIcon />
                      </InputGroupButton>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p>GCAs must be comma seperated</p>
                  </TooltipContent>
                </Tooltip>
              </InputGroupAddon>
            </InputGroup>

              </div>
              </div>
            <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">


              {baseFields.map(({ label, placeholder }, index) => (
                <div key={index}>
                  <Label htmlFor={label.toLowerCase().replace(" ", "-")}>
                    {label}
                  </Label>
                  <Input
                    id={label.toLowerCase().replace(" ", "-")}
                    type="text"
                    placeholder={placeholder}
                    className="mt-3 gap-2 bg-filter-input-bg dark:bg-transparent"
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

            </div>
          </div>

          <div className="mt-8 flex justify-end">
            <Button className="w-full sm:w-auto" size="lg" onClick={handleGetAnnotations} disabled={loading}>
              {loading ? (
                <>
                  <Loader2 className="animate-spin mr-2" />
                  Loading...
                </>
              ) : (
                "Get Annotations"
              )}
            </Button>
          </div>

          {errorMessage && (
            <Alert variant="destructive" className="mt-8">
              <Terminal />
              <AlertTitle>Heads up!</AlertTitle>
              <AlertDescription>{errorMessage}</AlertDescription>
            </Alert>
          )}

          {annotations.length > 0 && (
            <div className="mt-10 shadow-lg border border:border rounded-2xl">
              <div className="flex flex-col gap-4 px-4 py-6 border-b border:border sm:px-8 md:flex-row md:items-center md:justify-between">
                <h2 className="text-lg font-semibold">
                  Filtered Annotations ({annotations.length})
                </h2>
                <div className="grid w-full grid-cols-1 gap-2 sm:grid-cols-2 md:w-auto md:flex">
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(downloadables?.anno_main, "text/csv")}
                  >
                    Download CSV
                  </Button>

                  <Button
                    variant="outline"
                    onClick={() => handleDownload(downloadables?.anno_wide, "text/csv")}
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
