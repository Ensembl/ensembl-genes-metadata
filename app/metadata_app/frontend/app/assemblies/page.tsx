"use client";

import React, { useEffect, useState } from "react";
import { XCircle } from "lucide-react";
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

// 👇 Mock data
async function getMockAssemblies(): Promise<Assemblies[]> {
  return [
    {
      id: 1,
      bioproject_id: "PRJNA123456",
      associated_project: "Bird Genomes",
      gca: "GCA_000001405.28",
      scientific_name: "Anser anser",
      release_date: "2024-12-01",
      lowest_taxon_id: 8839,
      internal_clade: "anseriformes",
      asm_level: "Chromosome",
      asm_type: "haploid",
      asm_name: "Anser_anser_v1.0",
      refseq_accession: "GCF_000001405.39",
      is_current: "yes",
      contig_n50: 15400000,
      total_sequence_length: 1200000000,
    },
    {
      id: 2,
      bioproject_id: "PRJEB654321",
      associated_project: "Waterfowl Assembly",
      gca: "GCA_000002305.32",
      scientific_name: "Anas platyrhynchos",
      release_date: "2023-05-20",
      lowest_taxon_id: 8840,
      internal_clade: "anseriformes",
      asm_level: "Scaffold",
      asm_type: "diploid",
      asm_name: "Anas_platyrhynchos_v2.1",
      refseq_accession: "GCF_000002305.36",
      is_current: "no",
      contig_n50: 10200000,
      total_sequence_length: 1120000000,
    },
  ];
}

export default function Page() {
  const baseFields = [
    { label: "BioProject ID", placeholder: "PRJNA123456" },
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Release date", placeholder: "2024-12-31" },
  ];

  const [selectedMetrics, setSelectedMetrics] = useState<string[]>([]);
  const [toggleStates, setToggleStates] = useState<{ [key: string]: string[] }>({});
  const [metricValues, setMetricValues] = useState<{ [key: string]: string }>({});
  const [autoFillRequested, setAutoFillRequested] = useState(false);
  const [assemblies, setAssemblies] = useState<Assemblies[]>([]);

  const handleToggleChange = (metric: string, values: string[]) => {
    setToggleStates((prev) => ({ ...prev, [metric]: values }));
  };

  const handleGetResults = async () => {
    const data = await getMockAssemblies();
    setAssemblies(data);
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
  ];

    const hasToggleGroup = selectedMetrics.some((metric) =>
  toggleMetrics.some((tm) => tm.label === metric)
  );

  const hasInputFields = selectedMetrics.some((metric) =>
    !toggleMetrics.some((tm) => tm.label === metric)
  );


  return (
    <div className="min-h-screen justify-center flex flex-wrap align-items-center">
      <div className="container m-16 mt-10 max-w-7xl">
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
                    className="mt-3 gap-2 bg-background"
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
                <Switch id="reference" />
                <Label htmlFor="reference">Check reference status from NCBI</Label>
              </div>
              <div className="flex items-center space-x-2">
                <Switch id="transcript_check" />
                <Label htmlFor="transcript_check">Check transcriptomic data from ENA</Label>
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
                            <ToggleGroupItem key={option} value={option}>
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
                        className="mt-1 my-2"
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
          <Button size="lg" onClick={handleGetResults}>
            Get Results
          </Button>
        </div>

        {/* Results */}
        {assemblies.length > 0 && (
          <div className="mt-10 shadow-lg border border:border rounded-2xl">
            <div className="flex items-center justify-between px-8 py-6 border-b border:border">
              <h2 className="text-lg font-semibold">Filtered Assemblies</h2>
              <div className="flex gap-2">
                {/* Download GCA List (.txt) */}
                <Button
                  variant="secondary"
                  onClick={() => {
                    const gcaList = assemblies
                      .map((row) => row.gca) // replace with correct key if different
                      .join("\n");

                    const blob = new Blob([gcaList], { type: "text/plain" });
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement("a");
                    link.href = url;
                    link.download = "filtered_gca_list.txt";
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                  }}
                >
                  Download GCA List
                </Button>

                {/* Download Full Table (CSV, from frontend) */}
                <Button
                  variant="secondary"
                  onClick={() => {
                    const csv = [
                      Object.keys(assemblies[0]).join(","),
                      ...assemblies.map((row) =>
                        Object.values(row)
                          .map((value) =>
                            typeof value === "string" && value.includes(",")
                              ? `"${value}"`
                              : value
                          )
                          .join(",")
                      ),
                    ].join("\n");

                    const blob = new Blob([csv], { type: "text/csv" });
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement("a");
                    link.href = url;
                    link.download = "filtered_assemblies.csv";
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                  }}
                >
                  Download CSV
                </Button>

                {/* Download Full Table from Backend */}
                <Button
                  variant="secondary"
                  onClick={async () => {
                    const res = await fetch("/api/download-full-table"); // update to your actual route
                    const blob = await res.blob();
                    const url = URL.createObjectURL(blob);
                    const link = document.createElement("a");
                    link.href = url;
                    link.download = "full_table_filtered_assemblies.csv"; // or .xlsx depending on backend
                    document.body.appendChild(link);
                    link.click();
                    document.body.removeChild(link);
                  }}
                >
                  Download Full Table
                </Button>
              </div>
            </div>
            <DataTable columns={columns} data={assemblies} />
          </div>
        )}
      </div>
    </div>
  );
}