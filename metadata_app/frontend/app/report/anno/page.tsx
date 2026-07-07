"use client";

import React, { useRef, useState } from "react";
import {Loader2, Terminal} from "lucide-react";
import { useReactToPrint } from 'react-to-print';
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { DataTable } from "@/components/tables/data-table";
import { Report, columns } from "@/features/reports/annotation-columns";
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
import {CladeItem, RepClade} from "@/components/ui/repo_anno_clade";
import {CladeLiveItem, RepCladeLive} from "@/components/ui/repo_anno_clade_live";
import { cleanPayload, parseTaxonIds, splitProjectFilters } from "@/features/shared/filter-utils";
import { PROJECT_OPTIONS } from "@/features/shared/project-options";



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
  const [projectLiveData, setProjectLive] = useState<ProjectItem[]>([])
  const [errorMessage, setErrorMessage] = useState<string | null>(null);
  const [cladeData, setClade] = useState<CladeItem[]>([]);
  const [cladeDataLive, setCladeLive] = useState<CladeLiveItem[]>([]);



  const title = "Generate annotation report";
  const description =
    "Select a biodiversity project or enter a BioProject ID to generate an overview of annotations by Genebuild. Use the optional filters to further customize your report. Generate a table with annotations and download a PDF report.";

  const handleGetAnnotations = async () => {
    setErrorMessage(null);
    setReport([]);

    setLoading(true);
    try {
      const { bioprojectIds, groupNames } = splitProjectFilters(
        selectedProjects,
        baseFieldValues["BioProject ID"],
      );
      const taxonIdArray = parseTaxonIds(baseFieldValues["Taxon ID"]);


      const payload = {
        bioproject_id: bioprojectIds.length > 0 ? bioprojectIds : null,
        group_name: groupNames.length > 0 ? groupNames : null,
        taxon_id: taxonIdArray,
        start_date: baseFieldValues["Report start date"] || null,
        end_date: baseFieldValues["Report end date"] || null
      };

      const cleanedPayload = cleanPayload(payload);

      const res = await fetch("/api/report/anno/report/anno/filter", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(cleanedPayload),
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
        setProjectLive(result.project_report_live ?? []);
        setClade(result.clade_group);
        setCladeLive(result.clade_group_live);
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
    <div className="min-h-screen flex justify-center px-4 py-10 sm:px-6 lg:px-8">
      <div className="w-full max-w-6xl">
        <div className="rounded-2xl border-accent">
          <div className="rounded-2xl bg-secondary p-4 gap-10 shadow-lg sm:p-8">
            <div className="mb-8">
              <h1 className="text-2xl font-bold mb-2">{title}</h1>
              <p className="text-muted-foreground">{description}</p>
            </div>
            <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
              <div className="sm:col-span-2">
                <Label className="mb-3 block">Main projects</Label>
                <MultipleSelector
                  placeholder="Select projects or groups..."
                  defaultOptions={PROJECT_OPTIONS}
                  onChange={(values) => setSelectedProjects(values)}
                />
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
              <Button className="w-full cursor-pointer sm:w-auto" onClick={handleGetAnnotations} disabled={loading}>
                {loading ? <Loader2 className="animate-spin h-4 w-4 mr-2" /> : "Generate annotation report"}
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
          <div className="py-8 sm:p-8">
            <div className="flex flex-col gap-4 mt-4 mb-6 sm:flex-row sm:items-center sm:justify-between">
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
              <div ref={componentRef} className="grid grid-cols-1 gap-6 lg:grid-cols-2">
                <RepStatus data={statusData} />
                <RepProject data={projectData} />
                <RepProject
                  data={projectLiveData}
                  title="Live associated biodiversity projects"
                  description="Number of live annotations per project"
                />

                <AnnotatedBuscoCard data={buscoData} />
                <AnoMethodSummaryChart data={methodData} />

                <RepTopTaxa data={topTaxaData} />
                <AnnotatedTaxaCard data={taxaData} />

                  <div className="lg:col-span-2">
                <RepClade data={cladeData} />
                </div>
                  <div className="lg:col-span-2">
                <RepCladeLive data={cladeDataLive} />
                </div>
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
