"use client";

import React, { useRef, useState } from "react";
import {Loader2, Terminal} from "lucide-react";
import { useReactToPrint } from 'react-to-print';
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { DataTable } from "@/components/tables/data-table";
import { Report, columns } from "@/features/reports/assembly-columns";
import MultipleSelector, { Option } from "@/components/ui/multi_select";
import {AsmTaxaCard, NumTaxaItem} from "@/components/ui/rep_asm_num_taxa"
import {RepTopTaxa, TaxaItem} from "@/components/ui/rep_anno_top_taxa"
import {ProjectItem, RepProject} from "@/components/ui/repo_anno_project"
import {Card, CardContent} from "@/components/ui/card";
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuLabel,
  DropdownMenuItem,
  DropdownMenuTrigger,
  DropdownMenuSeparator
} from "@/components/ui/dropdown-menu";
import {Switch} from "@/components/ui/switch";
import {CladeItem, RepClade} from "@/components/ui/repo_asm_clade";
import {AsmTypeItem, RepAsmType} from "@/components/ui/repo_asm_type";
import {AsmLevelItem, RepAsmLevel} from "@/components/ui/repo_asm_level";
import {TranscItem, TranscCard} from "@/components/ui/rep_asm_transc";
import {LengthItem, LengthChart} from "@/components/ui/rep_asm_length";
import {RepTranscENA, TranscENAItem} from "@/components/ui/repo_asm_transc_ena";
import {Tooltip, TooltipContent, TooltipTrigger} from "@/components/ui/tooltip";
import {Alert, AlertDescription, AlertTitle} from "@/components/ui/alert";
import { cleanPayload, parseTaxonIds, splitProjectFilters } from "@/features/shared/filter-utils";
import { GROUP_NAME_VALUES, PROJECT_OPTIONS } from "@/features/shared/project-options";



type Downloadables = {
  rep_asm_main: string;
  rep_asm_wide: string;
  gca_list: string;
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
  const [assemblies, setReport] = useState<Report[]>([]);
  const [downloadables, setDownloadables] = useState<Downloadables | null>(null);
  const [loading, setLoading] = useState(false);
  const [candidate, setCandidate] = useState<boolean>(false);
  const [nonAnnotated, setNonAnnotated] = useState<boolean>(true);
  const [asmtypeData, setAsmType] = useState<AsmTypeItem[]>([]);
  const [asmlevelData, setAsmLevel] = useState<AsmLevelItem[]>([]);
  const [transcData, setTransc] = useState<TranscItem| null>(null);
  const [taxaData, setTaxa] = useState<NumTaxaItem | null>(null);
  const [topTaxaData, setTopTaxa] = useState<TaxaItem[]>([]);
  const [transcenaData, setTranscENA] = useState<TranscENAItem[]>([]);
  const [projectData, setProject] = useState<ProjectItem[]>([]);
  const [cladeData, setClade] = useState<CladeItem[]>([]);
  const [lengthData, setLength] = useState<LengthItem[]>([]);
  const [transc_ena, setENA] = useState<boolean>(false);
  const [transc, setTransc_check_reg] = useState<boolean>(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);

  const title = "Generate assembly report";
  const description =
    "Select a biodiversity project or enter a BioProject ID to generate an overview of assemblies. Use the optional filters to further customize your report. Generate a table with assemblies and download a PDF report.";

  const hasBioprojectInput =
  selectedProjects.some(item => !GROUP_NAME_VALUES.some((group) => group === item.value)) ||
  (baseFieldValues["BioProject ID"]?.trim() ?? "") !== "";

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
        end_date: baseFieldValues["Report end date"] || null,
        candidate: candidate,
        non_annotated: nonAnnotated,
        transc_ena: transc_ena,
        transc: transc,
      };

      const cleanedPayload = cleanPayload(payload);

      const res = await fetch("/api/report/asm/report/asm/filter", {
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

      if (result.rep_asm_main) {
        setReport(result.rep_asm_main);
        setAsmType(result.asm_type_group);
        setAsmLevel(result.asm_level_group);
        setTransc(result.transc_reg_count);
        setTaxa(result.num_unique_taxa);
        setTopTaxa(result.top_3_taxa);
        setProject(result.project_report);
        setClade(result.clade_group);
        setLength(result.asm_length);
        setTranscENA(result.transc_data);
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
    documentTitle: 'Assembly report',
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
              <div className="col-span-2">
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

              <div className="col-span-2 grid grid-cols-4 gap-10 mt-4 items-center">
                <Tooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center gap-2">
                      <Switch
                        checked={transc_ena}
                        onCheckedChange={setENA}
                      />
                      <Label>Check ENA for RNASeq</Label>
                    </div>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p>This will take longer to process</p>
                  </TooltipContent>
                </Tooltip>

                <div className="flex items-center gap-2">
                  <Switch
                    checked={transc}
                    onCheckedChange={setTransc_check_reg}
                  />
                  <Label>Check transcriptomic registry</Label>
                </div>

                <Tooltip>
                  <TooltipTrigger asChild>
                    <div className="flex items-center gap-2">
                      <Switch
                        checked={candidate}
                        onCheckedChange={setCandidate}
                      />
                      <Label>Annotation candidates</Label>
                    </div>
                  </TooltipTrigger>
                  <TooltipContent>
                    <p>Contig N50 min. 100.000, chromosome and complete genome.</p>
                  </TooltipContent>
                </Tooltip>

                <div className="flex items-center gap-2">
                  <Switch
                    checked={nonAnnotated}
                    onCheckedChange={setNonAnnotated}
                  />
                  <Label className="leading-tight">Only show non-annotated assemblies</Label>
                </div>
              </div>

            </div>
            <div className="mt-6 flex justify-end">
              <Button className="cursor-pointer" onClick={handleGetAnnotations} disabled={loading}>
                {loading ? <Loader2 className="animate-spin h-4 w-4 mr-2" /> : "Generate assembly report"}
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

        {assemblies.length > 0 && (
          <div className="p-8">
            <div className="flex items-center justify-between mt-4 mb-6">
              <h1 className="text-xl font-semibold">Assembly report</h1>
              <DropdownMenu>
                <DropdownMenuTrigger asChild>
                  <Button variant="outline">Download</Button>
                </DropdownMenuTrigger>
                <DropdownMenuContent className="w-56" align="end">
                  <DropdownMenuLabel>Download report</DropdownMenuLabel>
                  <DropdownMenuSeparator />
                  <DropdownMenuItem onClick={handlePrint}>
                    PDF</DropdownMenuItem>
                  <DropdownMenuItem onClick={() => handleDownload(downloadables?.rep_asm_wide, "assembly_report.csv", "text/csv")}
                  >CSV</DropdownMenuItem>
                  <DropdownMenuItem onClick={() => handleDownload(downloadables?.gca_list, "gca_report.csv", "text/csv")}
                  >GCA list</DropdownMenuItem>
                </DropdownMenuContent>
              </DropdownMenu>
            </div>
            <div className="grid grid-cols-1 gap-6">
              <div ref={componentRef} className="grid grid-cols-2 gap-6">
                {!hasBioprojectInput &&
                    <RepProject data={projectData} />}
                <RepAsmType data={asmtypeData} />
                <RepAsmLevel data={asmlevelData} />

                {baseFieldValues["Taxon ID"] && (
                  <>
                    <RepTopTaxa data={topTaxaData} />
                    <AsmTaxaCard data={taxaData} />
                  </>
                )}

                {transc &&  (
                  <TranscCard data={transcData} />
                  )}
                {transc_ena && (
                <RepTranscENA data={transcenaData} />
                )}

                <div className="col-span-2">
                <RepClade data={cladeData} />
                </div>

                <div className="col-span-2">
                <LengthChart data={lengthData} />
                </div>

              </div>



              {/* Full-width annotations table */}
              <div>
                <h2 className="text-xl font-semibold my-8">Assemblies table</h2>
                <Card>
                  <CardContent>
                    <DataTable columns={columns} data={assemblies} />
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
