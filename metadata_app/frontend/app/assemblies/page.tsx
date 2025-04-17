"use client";

import React, { useState } from "react";
import { CircleX } from "lucide-react";
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

// 👇 Mock data for now
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
  const [assemblies, setAssemblies] = useState<Assemblies[]>([]);

  const handleToggleChange = (metric: string, values: string[]) => {
    setToggleStates((prev) => ({ ...prev, [metric]: values }));
  };

  const handleGetResults = async () => {
    const data = await getMockAssemblies();
    setAssemblies(data);
  };

    const handleAutoFillHighQuality = () => {
    if (!selectedMetrics.includes("Contig N50")) {
      setSelectedMetrics((prev) => [...prev, "Contig N50"]);
    }
    if (!selectedMetrics.includes("Assembly level")) {
      setSelectedMetrics((prev) => [...prev, "Assembly level"]);
    }

    setToggleStates((prev) => ({
      ...prev,
      "Assembly level": ["Complete genome", "Chromosome"],
    }));

    // Optional: if you make the Contig N50 input controlled by state
    const input = document.getElementById("contig-n50") as HTMLInputElement;
    if (input) {
      input.value = "100000";
    }
  };

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
            <div className="bg-color-sidebar-accent dark:border-x-2 dark:border-b-2 rounded-b-2xl grid gap-6 p-8">
              <div className="space-y-6">
                {selectedMetrics.map((metric) =>
                  toggleMetrics.some((tm) => tm.label === metric) ? (
                    <div key={metric}>
                      <div className="flex items-center justify-start mb-2">
                        <Label className="block mr-2">{metric}</Label>
                        <Button
                          variant="ghost"
                          size="icon"
                          className="h-4 w-4 p-0 text-muted-foreground hover:text-foreground"
                          onClick={() =>
                            setSelectedMetrics((prev) =>
                              prev.filter((m) => m !== metric)
                            )
                          }
                        >
                          <CircleX className="h-3 w-3" strokeWidth={2.5} />
                        </Button>
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

              {/* Metric Input Fields */}
              <div className="grid grid-cols-4 gap-4">
                {selectedMetrics.map((metric) =>
                  !toggleMetrics.some((tm) => tm.label === metric) ? (
                    <div key={metric}>
                      <Label
                        className="mr-2"
                        htmlFor={metric.toLowerCase().replace(" ", "-")}
                      >
                        {metric}
                      </Label>
                      <Input
                        id={metric.toLowerCase().replace(" ", "-")}
                        type="text"
                        placeholder={`Enter threshold for ${metric}`}
                        className="mt-1 my-2"
                      />
                    </div>
                  ) : null
                )}
              </div>
            </div>
          )}
        </div>

        {/* Get Results Button */}
        <div className="mt-8 flex justify-end">
          <Button size="lg" onClick={handleGetResults}>
            Get Results
          </Button>
        </div>

        {/* Result Table */}
      {assemblies.length > 0 && (
        <div className="mt-10 shadow-lg border border:border rounded-2xl">
          <div className="flex items-center justify-between p-6 border-b border:border">
            <h2 className="text-lg font-semibold">Filtered Assemblies</h2>
            <Button
              variant="secondary"
              onClick={() => {
                // Convert data to CSV and trigger download
                const csv = [
                  Object.keys(assemblies[0]).join(","), // Header row
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
                link.download = "assemblies.csv";
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
              }}
                >
                    Download CSV
            </Button>
            </div>
                <DataTable columns={columns} data={assemblies} />
              </div>
            )}
      </div>
    </div>
  );
}