"use client";
import React, { useState } from "react";
import { CircleX } from "lucide-react"
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Switch } from "@/components/ui/switch";
import { PopoverWithMultiSelect } from "@/components/ui/metrics_select";
import {
  ToggleGroup,
  ToggleGroupItem,
} from "@/components/ui/assembly_toggle";

export default function App() {
  const baseFields = [
    { label: "BioProject ID", placeholder: "PRJNA123456" },
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Release date", placeholder: "2024-12-31" },
  ];

  const [selectedMetrics, setSelectedMetrics] = useState<string[]>([]);

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

  const [toggleStates, setToggleStates] = useState<{ [key: string]: string[] }>({});

  const handleToggleChange = (metric: string, values: string[]) => {
    setToggleStates((prev) => ({ ...prev, [metric]: values }));
  };

  return (
    <div className="min-h-screen">
      <div className="container-wrapper m-16 mt-10">
        <div className="rounded-2xl border-accent shadow-lg">
          {/* Base Fields */}
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
                />
              </div>
            </div>

            {/* Static Toggles */}
            <div className="flex gap-10 mt-6">
              <div className="flex items-center space-x-2">
                <Switch id="reference" />
                <Label htmlFor="reference">Check reference status from NCBI</Label>
              </div>
              <div className="flex items-center space-x-2">
                <Switch id="transcript_check" />
                <Label htmlFor="transcript_check">
                  Check transcriptomic data from ENA
                </Label>
              </div>
            </div>
          </div>

          {/* Metrics Section (only if any metric is selected) */}
          {selectedMetrics.length > 0 && (
            <div className="bg-color-sidebar-accent dark:border-x-2 dark:border-b-2 rounded-b-2xl grid gap-6 p-8">

              {/* Toggle Groups */}
              {selectedMetrics.some((metric) =>
                toggleMetrics.some((tm) => tm.label === metric)
              ) && (
                <div className="">
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
                              setSelectedMetrics((prev) => prev.filter((m) => m !== metric))
                            }
                          >
                            <CircleX className="h-3 w-3" strokeWidth={2.5} />
                          </Button>
                          </div>
                          <ToggleGroup
                            type="multiple"
                            className="flex flex-wrap gap-2"
                            value={toggleStates[metric] || []}
                            onValueChange={(value) => handleToggleChange(metric, value)}
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
                </div>
              )}

              {/* Input Metric Fields */}
              {selectedMetrics.some(
                (metric) => !toggleMetrics.some((tm) => tm.label === metric)
              ) && (
                <div className="">
                  <div className="grid grid-cols-4 gap-4">
                    {selectedMetrics.map((metric) =>
                      !toggleMetrics.some((tm) => tm.label === metric) ? (
                        <div key={metric}>
                          <div className="flex items-center justify-start mb-2">
                            <Label className="mr-2" htmlFor={metric.toLowerCase().replace(" ", "-")}>
                              {metric}
                            </Label>
                            <Button
                            variant="ghost"
                            size="icon"
                            className="h-4 w-4 p-0 text-muted-foreground hover:text-foreground"
                            onClick={() =>
                              setSelectedMetrics((prev) => prev.filter((m) => m !== metric))
                            }
                          >
                            <CircleX className="h-3 w-3" strokeWidth={2.5} />
                          </Button>
                          </div>
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
          )}

        </div>
        <div className="mt-8 flex justify-end">
            <Button size="lg">Get Results</Button>
          </div>
      </div>
    </div>
  );
}