"use client";

import React, { useEffect, useState } from "react";
import {XCircle, Loader2, Terminal} from "lucide-react";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Switch } from "@/components/ui/switch";
import { PopoverWithMultiSelect } from "@/components/ui/metrics_select";
import {
  ToggleGroup,
  ToggleGroupItem,
} from "@/components/ui/assembly_toggle";
import { Assemblies, columns } from "@/app/tables/assemblies_columns";
import { DataTable } from "@/app/tables/data-table";
import {cn} from "@/lib/utils";
import {StartAnnotationDialog} from "@/components/start_anno_dialog";
import {Tooltip, TooltipContent, TooltipTrigger} from "@/components/ui/tooltip";
import {Alert, AlertDescription, AlertTitle} from "@/components/ui/alert";



export default function Page() {
  const baseFields = [
    { label: "BioProject ID", placeholder: "PRJNA123456" },
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Release date", placeholder: "2024-12-31" },
  ];

  const [selectedMetrics, setSelectedMetrics] = useState<string[]>([]);
  const [baseFieldValues, setBaseFieldValues] = useState<{ [key: string]: string }>({});
  const [toggleStates, setToggleStates] = useState<{ [key: string]: string[] }>({});
  const [metricValues, setMetricValues] = useState<{ [key: string]: string }>({});
  const [autoFillRequested, setAutoFillRequested] = useState(false);
  const [assemblies, setAssemblies] = useState<Assemblies[]>([]);
  const [checkENA, setCheckENA] = useState(false);
  const [nonAnnotated, setNonAnnotated] = useState(false);
  const [checkTranscript, setCheckTranscript] = useState(false);
  const [checkCurrent, setCheckCurrent] = useState(false);
  const [downloadables, setDownloadables] = useState<{
    gca_list: string;
    df_main: string;
    df_wide: string
  } | null>(null);
  const [loading, setLoading] = useState(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);


  const isNumeric = (value: string) => /^\d+(\.\d+)?$/.test(value);

  const handleToggleChange = (metric: string, values: string[]) => {
    setToggleStates((prev) => ({ ...prev, [metric]: values }));
  };

  const handleGetResults = async () => {
      // Reset states at the beginning
    setErrorMessage(null);
    setAssemblies([]);

    // Validate numeric inputs before proceeding
    for (const [metric, value] of Object.entries(metricValues)) {
      if (value && !isNumeric(value)) {
        alert(`Invalid value for metric "${metric}". Please enter a numeric threshold.`);
        return; // Stop execution if validation fails
      }
    }

    setLoading(true);
    try {
      // Format BioProject ID as an array
      let bioprojectArray = null;
      const bioprojectId = baseFieldValues["BioProject ID"];

      if (bioprojectId) {
        if (bioprojectId.includes(',')) {
          bioprojectArray = bioprojectId.split(',').map(id => id.trim()).filter(id => id);
        } else {
          bioprojectArray = [bioprojectId.trim()];
        }
      }

      // Format metric thresholds
      const metric_thresholds: Record<string, number> = {};      Object.entries(metricValues)
        .filter(([_, value]) => value) // Only include fields with values
        .forEach(([metric, value]) => {
          metric_thresholds[metric] = Number(value);
        });

      // Get Assembly Level and Assembly Type values from toggleStates
      const asm_level = toggleStates["Assembly level"] || null;
      const asm_type = toggleStates["Assembly type"] || null;
      const pipeline = toggleStates["Pipeline"] || null;

      // Format taxon_id as number
      let taxonIdArray = null;
      const taxonInput = baseFieldValues["Taxon ID"];

      if (taxonInput) {
        if (taxonInput.includes(',')) {
          taxonIdArray = taxonInput
            .split(',')
            .map(id => parseInt(id.trim(), 10))
            .filter(id => !isNaN(id));
        } else {
          const parsed = parseInt(taxonInput.trim(), 10);
          if (!isNaN(parsed)) {
            taxonIdArray = [parsed];
          }
        }
      }

      // Format the payload according to API expectations
      const payload = {
        bioproject_id: bioprojectArray,
        metric_thresholds: Object.keys(metric_thresholds).length > 0 ? metric_thresholds : null,
        asm_level: asm_level,
        asm_type: asm_type,
        release_date: baseFieldValues["Release date"] || null,
        taxon_id: taxonIdArray,
        current: checkCurrent,
        pipeline: pipeline,
        transc: checkTranscript,
        transc_ena: checkENA,
        non_annotated: nonAnnotated,
      };

      // Remove null fields to keep payload clean
      const cleanPayload = Object.fromEntries(
      Object.entries(payload).filter(([_, value]) => value !== null && value !== undefined)
      );

      console.log("Sending payload:", JSON.stringify(payload));

      // Update the endpoint URL to match your API
      const res = await fetch("http://127.0.0.1:8000/api/assemblies/assemblies/filter", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "accept": "application/json"
        },
        body: JSON.stringify(payload),
      });

      if (!res.ok) {
        const errorText = await res.text();
        console.error("Failed to fetch report", errorText);
        setErrorMessage("Failed to fetch report: " + errorText);
        return;
      }

      const result = await res.json();
      console.log("API response:", result);

      if (result.df_main) {
        setAssemblies(result.df_main);
      } else {
        console.error("No assemblies data in response");
        alert("No assemblies data found in response");
      }

      if (result.downloadables) {
        console.log("Downloadables from backend:", result.downloadables);
        setDownloadables(result.downloadables);
      } else {
        console.error("No downloadable data in response");
      }
    } catch (error) {
      console.error("Error fetching assemblies:", error);
      setErrorMessage("Error fetching data: " + (error instanceof Error ? error.message : String(error)));
    } finally {
      setLoading(false);
    }
  };

  const handleDownload = (content: string | undefined, filename: string, type: string = "text/plain") => {
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

  const handleAutoFillHighQuality = () => {
    setSelectedMetrics((prev) => {
      const updated = [...prev];
      if (!updated.includes("Contig N50")) updated.push("Contig N50");
      if (!updated.includes("Assembly level")) updated.push("Assembly level");
      return updated;
    });
    setAutoFillRequested(true);
  };

  useEffect(() => {
    if (
      autoFillRequested &&
      selectedMetrics.includes("Contig N50") &&
      selectedMetrics.includes("Assembly level")
    ) {
      setMetricValues((prev) => ({
        ...prev,
        "Contig N50": "100000",
      }));
      setToggleStates((prev) => ({
        ...prev,
        "Assembly level": ["Complete genome", "Chromosome"],
      }));
      setAutoFillRequested(false);
    }
  }, [autoFillRequested, selectedMetrics]);

  const toggleMetrics = [
    {
      label: "Assembly level",
      options: ["Contig", "Scaffold", "Chromosome", "Complete genome"],
    },
      {
      label: "Pipeline",
      options: ["anno", "main", "hprc"],
    },
    {
      label: "Assembly type",
      options: [
        "haploid",
        "alternate-pseudohaplotype",
        "unresolved-diploid",
        "haploid-with-alt-loci",
        "diploid",
      ],
    },
  ];

  const hasToggleGroup = selectedMetrics.some((metric) =>
    toggleMetrics.some((tm) => tm.label === metric)
  );

  const hasInputFields = selectedMetrics.some((metric) =>
    !toggleMetrics.some((tm) => tm.label === metric)
  );

  return (
    <div className="min-h-screen justify-center flex flex-wrap align-items-center">
      <div className="container m-16 mt-10 max-w-6xl">
        <div className="rounded-2xl border-accent shadow-lg">
          {/* Filter Section */}
          <div className="rounded-t-2xl bg-secondary p-8 gap-10">
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
              <div className="flex items-end">
                <PopoverWithMultiSelect
                  selectedItems={selectedMetrics}
                  setSelectedItems={setSelectedMetrics}
                  onAutoFillHighQuality={handleAutoFillHighQuality}
                />
              </div>
            </div>

            <div className="flex gap-10 mt-6">

              <div className="flex items-center space-x-2">
                <Tooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center space-x-2">
                <Switch id="ena" checked={checkENA} onCheckedChange={setCheckENA} />
                <Label htmlFor="ena">Check ENA</Label>
                      </div>
                </TooltipTrigger>
                <TooltipContent>
                  <p>This will take longer to process</p>
                </TooltipContent>
              </Tooltip>
              </div>

              <div className="flex items-center space-x-2">
                <Switch id="transcript_check" checked={checkTranscript} onCheckedChange={setCheckTranscript} />
                <Label htmlFor="transcript_check">Check transcriptomic registry</Label>
              </div>
              <div className="flex items-center space-x-2">
                <Switch id="current_check" checked={checkCurrent} onCheckedChange={setCheckCurrent} />
                <Label htmlFor="current_check">Only show current records</Label>
              </div>
              <div className="flex items-center space-x-2">
                <Switch id="non_annotated" checked={nonAnnotated} onCheckedChange={setNonAnnotated} />
                <Label htmlFor="non_annotated">Only show non-annotated assemblies</Label>
              </div>
            </div>
          </div>

          {/* Metric Toggles */}
          {selectedMetrics.length > 0 && (
          <div
            className={cn(
              "bg-color-sidebar-accent dark:border-x-2 dark:border-b-2 rounded-b-2xl flex flex-col p-8",
              hasToggleGroup && hasInputFields && "space-y-8"
            )}
          >
            {hasToggleGroup && (
              <div className="space-y-6">
                {selectedMetrics.map((metric) =>
                  toggleMetrics.some((tm) => tm.label === metric) ? (
                    <div key={metric}>
                      <div className="flex items-center justify-start mb-2">
                        <Label
                        className="gap-1"
                        htmlFor={metric.toLowerCase().replace(" ", "-")}
                      >
                        {metric}
                        <XCircle
                          className="h-5 w-5 cursor-pointer text-muted-foreground hover:text-foreground"
                          strokeWidth={2.5}
                          fill="currentColor"
                          color="background"
                          onClick={() =>
                            setSelectedMetrics((prev) => prev.filter((m) => m !== metric))
                          }
                        />
                      </Label>
                      </div>
                      <ToggleGroup
                        type="multiple"
                        className="flex flex-wrap gap-2"
                        value={toggleStates[metric] || []}
                        onValueChange={(value) =>
                          handleToggleChange(metric, value)
                        }
                      >
                        {toggleMetrics
                          .find((tm) => tm.label === metric)!
                          .options.map((option) => (
                            <ToggleGroupItem key={option} value={option} variant="outline_filter">
                              {option}
                            </ToggleGroupItem>
                          ))}
                      </ToggleGroup>
                    </div>
                  ) : null
                )}
              </div>
            )}

            {hasInputFields && (
              <div className="grid grid-cols-4 gap-4">
                {selectedMetrics.map((metric) =>
                  !toggleMetrics.some((tm) => tm.label === metric) ? (
                    <div key={metric}>
                      <Label
                        className="gap-1"
                        htmlFor={metric.toLowerCase().replace(" ", "-")}
                      >
                        {metric}
                        <XCircle
                          className="h-5 w-5 cursor-pointer text-muted-foreground hover:text-foreground"
                          strokeWidth={2.5}
                          fill="currentColor"
                          color="background"
                          onClick={() =>
                            setSelectedMetrics((prev) => prev.filter((m) => m !== metric))
                          }
                        />
                      </Label>
                      <Input
                        id={metric.toLowerCase().replace(" ", "-")}
                        type="text"
                        placeholder={`Enter threshold for ${metric}`}
                        className="mt-1 my-2 bg-transparent"
                        value={metricValues[metric] || ""}
                        onChange={(e) =>
                          setMetricValues((prev) => ({
                            ...prev,
                            [metric]: e.target.value,
                          }))
                        }
                      />
                    </div>
                  ) : null
                )}
              </div>
            )}
          </div>
        )}
        </div>

        {/* Get Results Button */}
        <div className="mt-8 flex justify-end">
          <Button size="lg" onClick={handleGetResults} disabled={loading}>
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

        {errorMessage && (
            <Alert variant="destructive" className="mt-8">
              <Terminal />
              <AlertTitle>Heads up!</AlertTitle>
              <AlertDescription>
                {errorMessage}
              </AlertDescription>
            </Alert>
        )}

        {/* Results */}
        {assemblies.length > 0 &&  (
          <div className="mt-10 shadow-lg border border:border rounded-2xl">
            <div className="flex items-center justify-between px-8 py-6 border-b border:border">
              <h2 className="text-lg font-semibold">Filtered Assemblies ({assemblies.length}) </h2>
              <div className="flex gap-2">
                {/* Download GCA List (.txt) */}
                <Button
                  variant="outline"
                  onClick={() => handleDownload(downloadables?.gca_list, "filtered_gca_list.txt")}
                >
                  Download GCA List
                </Button>

                {/* Download main Table (CSV, from backend) */}
                <Button
                  variant="outline"
                  onClick={() => handleDownload(downloadables?.df_main, "df_main.csv", "text/csv;charset=utf-8;")}
                >
                  Download CSV
                </Button>

                {/* Download Full Table from Backend */}
                <Button
                  variant="outline"
                  onClick={() => handleDownload(downloadables?.df_wide, "full_table_filtered_assemblies.csv", "text/csv")}
                >
                  Download Full Table
                </Button>

                {/* Annotation start button */}
                <StartAnnotationDialog/>

              </div>
            </div>
            <DataTable columns={columns} data={assemblies} />
          </div>
        )}
      </div>
    </div>
  );
}