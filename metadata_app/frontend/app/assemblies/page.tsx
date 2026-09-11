"use client";

import React, { useEffect, useMemo, useState } from "react";
import { type ColumnDef, type RowSelectionState } from "@tanstack/react-table";
import {XCircle, Loader2, Terminal, InfoIcon, DownloadIcon, Sparkles} from "lucide-react";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Switch } from "@/components/ui/switch";
import { Checkbox } from "@/components/ui/checkbox";
import { CopyButton } from "@/components/ui/button-copy";
import { PopoverWithMultiSelect } from "@/components/ui/metrics_select";
import MultipleSelector, { Option } from "@/components/ui/multi_select";
import {
  ToggleGroup,
  ToggleGroupItem,
} from "@/components/ui/assembly_toggle";
import { Assemblies, columns } from "@/features/assemblies/columns";
import { DataTable } from "@/components/tables/data-table";
import {cn} from "@/lib/utils";
import {Tooltip, TooltipContent, TooltipTrigger} from "@/components/ui/tooltip";
import {Alert, AlertDescription, AlertTitle} from "@/components/ui/alert";
import {InputGroup, InputGroupAddon, InputGroupButton, InputGroupInput} from "@/components/ui/input-group";
import { cleanPayload, parseTaxonIds, splitCommaSeparated, splitProjectFilters } from "@/features/shared/filter-utils";
import { PROJECT_OPTIONS } from "@/features/shared/project-options";



export default function Page() {
  const baseFields = [
    { label: "BioProject ID", placeholder: "PRJNA123456"},
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Release date", placeholder: "2024-12-31" },
  ];

  const [selectedProjects, setSelectedProjects] = useState<Option[]>([]);
  const [selectedMetrics, setSelectedMetrics] = useState<string[]>([]);
  const [baseFieldValues, setBaseFieldValues] = useState<{ [key: string]: string }>({});
  const [toggleStates, setToggleStates] = useState<{ [key: string]: string[] }>({});
  const [metricValues, setMetricValues] = useState<{ [key: string]: string }>({});
  const [autoFillRequested, setAutoFillRequested] = useState(false);
  const [assemblies, setAssemblies] = useState<Assemblies[]>([]);
  const [checkENA, setCheckENA] = useState(false);
  const [nonAnnotated, setNonAnnotated] = useState(true);
  const [checkTranscript, setCheckTranscript] = useState(false);
  const [checkCurrent, setCheckCurrent] = useState(true);
  const [gcaInput, setGcaInput] = useState<string>("");
  const [downloadables, setDownloadables] = useState<{
    gca_list: string;
    df_wide: string
  } | null>(null);
  const [loading, setLoading] = useState(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);
  const [filterShortReadPresent, setFilterShortReadPresent] = useState(true);
  const [filterLongReadPresent, setFilterLongReadPresent] = useState(false);
  const [rowSelection, setRowSelection] = useState<RowSelectionState>({});

  const isNumeric = (value: string) => /^\d+(\.\d+)?$/.test(value);

  const handleToggleChange = (metric: string, values: string[]) => {
    setToggleStates((prev) => ({ ...prev, [metric]: values }));
  };

  const handleGetResults = async () => {
    // Reset states at the beginning
    setErrorMessage(null);
    setAssemblies([]);
    setRowSelection({});

    // Validate numeric inputs before proceeding
    for (const [metric, value] of Object.entries(metricValues)) {
      if (value && !isNumeric(value)) {
        alert(`Invalid value for metric "${metric}". Please enter a numeric threshold.`);
        return; // Stop execution if validation fails
      }
    }

    setLoading(true);

    try {
      const { bioprojectIds, groupNames } = splitProjectFilters(
        selectedProjects,
        baseFieldValues["BioProject ID"],
      );

      const metric_thresholds: Record<string, number> = {};
      Object.entries(metricValues)
        .filter(([, value]) => value) // Only include fields with values
        .forEach(([metric, value]) => {
          metric_thresholds[metric] = Number(value);
        });

      // Get Assembly Level, Assembly Type, and Pipeline values from toggleStates
      const asm_level = toggleStates["Assembly level"] || null;
      const asm_type = toggleStates["Assembly type"] || null;
      const pipeline = (toggleStates["Pipeline"] || []).map((value) => {
        return value.toLowerCase();
      });

      const taxonIdArray = parseTaxonIds(baseFieldValues["Taxon ID"]);
      const uniqueGCA = splitCommaSeparated(gcaInput);

      // Format the payload according to API expectations
      const payload = {
        bioproject_id: bioprojectIds.length > 0 ? bioprojectIds : null,
        group_name: groupNames.length > 0 ? groupNames : null,
        metric_thresholds: Object.keys(metric_thresholds).length > 0 ? metric_thresholds : null,
        asm_level: asm_level,
        asm_type: asm_type,
        pipeline: pipeline.length > 0 ? pipeline : null,
        release_date: baseFieldValues["Release date"] || null,
        taxon_id: taxonIdArray,
        current: checkCurrent,
        transc: checkTranscript,
        transc_ena: checkENA,
        non_annotated: nonAnnotated,
        gca: uniqueGCA.length > 0 ? uniqueGCA : null,
      };


      console.log("Sending payload:", JSON.stringify(payload));

      // Call API
      const res = await fetch("/api/assemblies/assemblies/filter", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "accept": "application/json"
        },
        body: JSON.stringify(cleanPayload(payload)),
      });

      if (!res.ok) {
        const errorText = await res.text();
        console.error("Failed to fetch report", errorText);
        setErrorMessage("Failed to fetch report: " + errorText);
        return;
      }

      const result = await res.json();
      console.log("API response:", result);

      if (result.df_wide) {
        setAssemblies(result.df_wide);
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
      label: "Assembly type",
      options: [
        "haploid",
        "alternate-pseudohaplotype",
        "unresolved-diploid",
        "haploid-with-alt-loci",
        "diploid",
      ],
    },
    {
      label: "Pipeline",
      options: ["Main", "Annotation", "HPRC"],
    },
  ];

  const hasToggleGroup = selectedMetrics.some((metric) =>
    toggleMetrics.some((tm) => tm.label === metric)
  );

  const hasInputFields = selectedMetrics.some((metric) =>
    !toggleMetrics.some((tm) => tm.label === metric)
  );

  const filteredAssemblies = useMemo(() => {
    if (!checkENA) {
      return assemblies;
    }

    return assemblies.filter((assembly) => {
      const hasShortReadData =
        Number(assembly.short_read_paired_end_illumina_lowest ?? 0) > 0 ||
        Number(assembly.short_read_paired_end_illumina ?? 0) > 0;
      const hasLongReadData = Number(assembly.long_read_pacbio ?? 0) > 0;

      if (filterShortReadPresent && !hasShortReadData) {
        return false;
      }

      if (filterLongReadPresent && !hasLongReadData) {
        return false;
      }

      return true;
    });
  }, [assemblies, checkENA, filterLongReadPresent, filterShortReadPresent]);

  const selectableColumns = useMemo<ColumnDef<Assemblies>[]>(() => [
    {
      id: "select",
      header: ({ table }) => (
        <Checkbox
          checked={
            table.getIsAllPageRowsSelected() ||
            (table.getIsSomePageRowsSelected() && "indeterminate")
          }
          onCheckedChange={(value) => table.toggleAllPageRowsSelected(Boolean(value))}
          aria-label="Select all rows"
        />
      ),
      cell: ({ row }) => (
        <Checkbox
          checked={row.getIsSelected()}
          onCheckedChange={(value) => row.toggleSelected(Boolean(value))}
          aria-label={`Select ${row.original.gca}`}
        />
      ),
      enableSorting: false,
    },
    ...columns,
  ], []);

  const selectedAssemblies = useMemo(
    () => filteredAssemblies.filter((assembly) => rowSelection[assembly.gca]),
    [filteredAssemblies, rowSelection],
  );
  const selectedGcaText = useMemo(
    () => selectedAssemblies.map((assembly) => assembly.gca).join("\n"),
    [selectedAssemblies],
  );
  const quickFiltersActive =
    checkENA && (filterShortReadPresent || filterLongReadPresent);
  const quickFiltersEmptiedTable =
    assemblies.length > 0 && filteredAssemblies.length === 0 && quickFiltersActive;

  return (
    <div className="min-h-screen flex justify-center px-4 py-10 sm:px-6 lg:px-8">
      <div className="w-full max-w-6xl">
        <div className="rounded-2xl border-accent shadow-lg">
          {/* Filter Section */}
          <div className="rounded-t-2xl bg-secondary p-4 gap-10 sm:p-8">
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
              <div className="flex items-end">
                <PopoverWithMultiSelect
                  selectedItems={selectedMetrics}
                  setSelectedItems={setSelectedMetrics}
                  onAutoFillHighQuality={handleAutoFillHighQuality}
                />
              </div>
            </div>

            <div className="grid grid-cols-1 gap-4 mt-6 sm:grid-cols-2 lg:grid-cols-4">

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
                <Label htmlFor="current_check">Only show current assemblies</Label>
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
              "bg-color-sidebar-accent dark:border-x-2 dark:border-b-2 rounded-b-2xl flex flex-col p-4 sm:p-8",
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
                          className="h-5 w-5 cursor-pointer text-muted-foreground hover:text-primary"
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
              <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">
                {selectedMetrics.map((metric) =>
                  !toggleMetrics.some((tm) => tm.label === metric) ? (
                    <div key={metric}>
                      <Label
                        className="gap-1"
                        htmlFor={metric.toLowerCase().replace(" ", "-")}
                      >
                        {metric}
                        <XCircle
                          className="h-5 w-5 cursor-pointer text-muted-foreground hover:text-primary"
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
          <Button className="w-full cursor-pointer sm:w-auto" size="lg" onClick={handleGetResults} disabled={loading}>
            {loading ? (
              <>
                <Loader2 className="animate-spin mr-2" />
                Loading...
              </>
            ) : (
              <>
                <Sparkles className="mr-2" />
                Get Assemblies
              </>
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
            <div className="flex flex-col gap-4 px-4 py-6 border-b border:border sm:px-8 md:flex-row md:items-center md:justify-between">
              <h2 className="text-lg font-semibold">Filtered Assemblies ({filteredAssemblies.length}) </h2>
              <div className="grid w-full grid-cols-1 gap-2 sm:grid-cols-2 md:w-auto md:flex">
                {selectedAssemblies.length > 0 ? (
                  <CopyButton
                    text={selectedGcaText}
                    defaultLabel="Copy selected GCAs"
                  />
                ) : (
                  <Button
                    variant="outline"
                    onClick={() => handleDownload(downloadables?.gca_list, "filtered_gca_list.txt")}
                  >
                    Download GCA List
                  </Button>
                )}


                {/* Download Full Table from Backend */}
                <Button
                  variant="outline"
                  onClick={() => handleDownload(downloadables?.df_wide, "full_table_filtered_assemblies.csv", "text/csv")}
                >
                  <DownloadIcon></DownloadIcon>
                  Download Full Table
                </Button>

              </div>
            </div>
            {checkENA && (
              <div className="border-b border:border px-4 py-4 sm:px-8">
                <div className="flex flex-col gap-4 sm:flex-row sm:flex-wrap sm:items-center">
                  <div className="flex items-center space-x-2">
                    <Switch
                      id="quick-filter-short-read"
                      checked={filterShortReadPresent}
                      onCheckedChange={setFilterShortReadPresent}
                    />
                    <Label htmlFor="quick-filter-short-read">
                      Only show records with short read data
                    </Label>
                  </div>
                  <div className="flex items-center space-x-2">
                    <Switch
                      id="quick-filter-long-read"
                      checked={filterLongReadPresent}
                      onCheckedChange={setFilterLongReadPresent}
                    />
                    <Label htmlFor="quick-filter-long-read">
                      Only show records with long read data
                    </Label>
                  </div>
                </div>
              </div>
            )}
            {quickFiltersEmptiedTable && (
              <div className="border-b border:border px-4 py-4 sm:px-8">
                <Alert>
                  <InfoIcon />
                  <AlertTitle>Warning</AlertTitle>
                  <AlertDescription>
                    The quick filters are currently on, so all rows are hidden. Turn off one or both quick filters to see the matching assemblies again.
                  </AlertDescription>
                </Alert>
              </div>
            )}
            <DataTable
              columns={selectableColumns}
              data={filteredAssemblies}
              enableRowSelection
              getRowId={(row) => row.gca}
              onRowSelectionChange={setRowSelection}
              rowSelection={rowSelection}
            />
          </div>
        )}
      </div>
    </div>
  );
}
