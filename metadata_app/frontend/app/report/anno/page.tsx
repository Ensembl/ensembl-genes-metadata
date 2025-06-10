"use client";

import React, { useRef, useState } from "react";
import {Loader2, Terminal} from "lucide-react";
import { useReactToPrint } from 'react-to-print';
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import {
  ToggleGroup,
  ToggleGroupItem,
} from "@/components/ui/toggle-group";
import { DataTable } from "@/app/tables/data-table";
import { Report, columns } from "@/app/tables/report_columns";
import MultipleSelector, { Option } from "@/components/ui/multi_select";
import {RepStatus, StatusItem} from "@/components/ui/rep_anno_status";
import {AnnotatedBuscoCard, BuscoItem} from "@/components/ui/rep_anno_busco"
import {MethodItem, AnoMethodSummaryChart } from "@/components/ui/rep_anno_method"
import {AnnotatedTaxaCard, NumTaxaItem} from "@/components/ui/rep_anno_num_taxa"
import {RepTopTaxa, TaxaItem} from "@/components/ui/rep_anno_top_taxa"
import {ProjectItem, RepProject} from "@/components/ui/repo_anno_project"
import {Card, CardContent} from "@/components/ui/card";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert"
import {DropdownMenu, DropdownMenuContent, DropdownMenuLabel, DropdownMenuItem, DropdownMenuTrigger} from "@/components/ui/dropdown-menu";


type Downloadables = {
  main_report: string;
  anno_wide: string;
};

export default function Page() {
  const baseFields = [
    { label: "BioProject ID", placeholder: "PRJNA123456" },
    { label: "Taxon ID", placeholder: "9606" },
    { label: "Report start date", placeholder: "2024-12-31" },
    { label: "Report end date", placeholder: "2024-03-05" },
  ];

  const projectOptions: Option[] = [
    { value: "PRJEB40665", label: "Darwin Tree of Life" },
    { value: "PRJEB61747", label: "European Reference Genome Atlas/Biodiversity Genomics Europe" },
    { value: "PRJEB43510", label: "European Reference Genome Atlas" },
    { value: "PRJNA533106", label: "Earth BioGenome" },
    { value: "PRJEB47820", label: "European Reference Genome Atlas pilot" },
    { value: "PRJEB43743", label: "Aquatic Symbiosis" },
    { value: "PRJNA489243", label: "Vertebrate Genomes" },
    { value: "PRJNA813333", label: "Canadian BioGenome" },
    { value: "farmed_animals_2023", label: "Farmed animals" },
    { value: "AQUA-FAANG", label: "Aqua FAANG" },
  ];

  const [selectedProjects, setSelectedProjects] = useState<Option[]>([]);
  const [baseFieldValues, setBaseFieldValues] = useState<{ [key: string]: string }>({});
  const [annotations, setReport] = useState<Report[]>([]);
  const [downloadables, setDownloadables] = useState<Downloadables | null>(null);
  const [loading, setLoading] = useState(false);
  const [statusData, setStatus] = useState<StatusItem[]>([])
  const [methodData, setMethod] = useState<MethodItem[]>([])
  const [buscoData, setBusco] = useState<BuscoItem | null>(null);
  const [taxaData, setTaxa] = useState<NumTaxaItem | null>(null);
  const [topTaxaData, setTopTaxa] = useState<TaxaItem[]>([])
  const [projectData, setProject] = useState<ProjectItem[]>([])
  const [releaseSites, setReleaseSites] = useState<string>("");
  const [errorMessage, setErrorMessage] = useState<string | null>(null);


  const title = "Generate annotation report";
  const description =
    "Select a biodiversity project or enter a BioProject ID to generate an overview of annotations by Genebuild. Use the optional filters to further customize your report. Generate a table with annotations and download a PDF report.";

  const groupNameValues = ["farmed_animals_2023", "AQUA-FAANG"];

  const handleGetAnnotations = async () => {
    setLoading(true);
    try {
      const bioprojectArray: string[] = [];
      const groupNames: string[] = [];

      // Parse selection from dropdown
      selectedProjects.forEach((item) => {
        if (groupNameValues.includes(item.value)) {
          groupNames.push(item.value);
        } else {
          bioprojectArray.push(item.value);
        }
      });

      // Include manually entered BioProject IDs
      const manualIdInput = baseFieldValues["BioProject ID"];
      if (manualIdInput) {
        const manualIds = manualIdInput
          .split(",")
          .map((id) => id.trim())
          .filter((id) => id);
        bioprojectArray.push(...manualIds);
      }

      // Remove duplicates
      const uniqueBioprojects = Array.from(new Set(bioprojectArray));

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


      let release_type: string[] | null = null;
      if (releaseSites && releaseSites !== "both") {
        release_type = [releaseSites];
      }

      const payload = {
        bioproject_id: uniqueBioprojects.length > 0 ? uniqueBioprojects : null,
        group_name: groupNames.length > 0 ? groupNames : null,
        taxon_id: taxonIdArray,
        start_date: baseFieldValues["Report start date"] || null,
        end_date: baseFieldValues["Report end date"] || null,
        release_type: release_type,
      };

      const cleanPayload = Object.fromEntries(
        Object.entries(payload).filter(([_, value]) => value !== null && value !== undefined)
      );

      const res = await fetch("http://127.0.0.1:8000/api/report/anno/report/anno/filter", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(cleanPayload),
      });

      if (!res.ok) {
        const errorText = await res.text();
        console.error("Failed to fetch report", errorText);
        setErrorMessage("Failed to fetch report: " + errorText);
        return;
      }

      const result = await res.json();
      console.log("API response:", result);

      if (result.main_report) {
        setReport(result.main_report);
        setStatus(result.number_of_annotations);
        setMethod(result.method_report);
        setBusco(result.average_busco);
        setTaxa(result.num_unique_taxa);
        setTopTaxa(result.top_3_taxa);
        setProject(result.project_report);
      } else {
        alert("No data found.");
      }

      if (result.downloadables_report) {
        setDownloadables(result.downloadables_report);
      } else if (result.report) {
        setDownloadables(result.report);
      }
    } catch (error) {
      console.error("Error fetching data:", error);
      setErrorMessage("Error fetching data: " + (error instanceof Error ? error.message : String(error)));
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
      alert(`Error downloading ${filename}: ${error instanceof Error ? error.message : String(error)}`);
    }
  };

  const componentRef = useRef<HTMLDivElement>(null);

  const handlePrint = useReactToPrint({
    contentRef: componentRef,
    documentTitle: 'Annotation report',
    pageStyle: `
      @page {
        size: A4;
        margin: 20mm;
      }
      @media print {
        body { -webkit-print-color-adjust: exact; }
        .print-content { 
          display: block !important; 
        }
      }
    `
  });

  return (
    <div className="min-h-screen justify-center flex flex-wrap align-items-center">
      <div className="container m-16 mt-10 max-w-6xl">
        <div className="rounded-2xl border-accent">
          <div className="rounded-2xl bg-secondary p-8 gap-10 shadow-lg">
            <div className="mb-8">
              <h1 className="text-2xl font-bold mb-2">{title}</h1>
              <p className="text-muted-foreground">{description}</p>
            </div>
            <div className="grid justify-center grid-cols-2 gap-4">
              <div>
                <Label className="mb-3 block">Main projects</Label>
                <MultipleSelector
                  placeholder="Select projects or groups..."
                  defaultOptions={projectOptions}
                  onChange={(values) => setSelectedProjects(values)}
                />
              </div>

              <div>
                <Label className="mb-2 block">Release site</Label>
                <ToggleGroup
                  type="single"
                  size="lg"
                  variant="outline"
                  value={releaseSites}
                  onValueChange={setReleaseSites}
                >
                  <ToggleGroupItem value="main" aria-label="main">
                    Main
                  </ToggleGroupItem>
                  <ToggleGroupItem value="beta" aria-label="beta">
                    Beta
                  </ToggleGroupItem>
                  <ToggleGroupItem value="both" aria-label="both">
                    Both
                  </ToggleGroupItem>
                </ToggleGroup>
              </div>

              {baseFields.map(({ label, placeholder }, index) => (
                <div key={index}>
                  <Label htmlFor={label.toLowerCase().replace(" ", "-")}>{label}</Label>
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
            <div className="mt-6 flex justify-end">
              <Button onClick={handleGetAnnotations} disabled={loading}>
                {loading ? <Loader2 className="animate-spin h-4 w-4 mr-2" /> : "Generate Report"}
              </Button>
            </div>
          </div>
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

        {annotations.length > 0 && (
          <div className="p-8">
            <div className="flex items-center justify-between mt-4 mb-6">
              <h1 className="text-xl font-semibold">Annotations report</h1>
              <DropdownMenu>
                <DropdownMenuTrigger asChild>
                  <Button variant="outline">Download</Button>
                </DropdownMenuTrigger>
                <DropdownMenuContent className="w-56" align="end">
                  <DropdownMenuLabel>Download report</DropdownMenuLabel>
                  <DropdownMenuItem onClick={handlePrint}>
                    PDF</DropdownMenuItem>
                  <DropdownMenuItem onClick={() => handleDownload(downloadables?.anno_wide, "annotations_report.csv", "text/csv")}
                  >CSV</DropdownMenuItem>
                </DropdownMenuContent>
              </DropdownMenu>
            </div>
            <div className="grid grid-cols-1 gap-6">
              {/* Two-column card layout */}
              <div ref={componentRef} className="grid grid-cols-2 gap-6">
                <RepStatus data={statusData} />
                <RepProject data={projectData} />

                <AnnotatedBuscoCard data={buscoData} />
                <AnoMethodSummaryChart data={methodData} />

                <RepTopTaxa data={topTaxaData} />
                <AnnotatedTaxaCard data={taxaData} />
              </div>

              {/* Full-width annotations table */}
              <div>
                <h2 className="text-xl font-semibold my-8">Annotations table</h2>
                <Card>
                  <CardContent>
                    <DataTable columns={columns} data={annotations} />
                  </CardContent>
                </Card>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}