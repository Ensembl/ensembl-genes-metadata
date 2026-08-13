"use client";

import React, { useEffect, useMemo, useState } from "react";
import { DownloadIcon, Loader2, Sparkles, Terminal } from "lucide-react";

import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Button } from "@/components/ui/button";
import { DataTable } from "@/components/tables/data-table";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { cleanPayload, parseTaxonIds } from "@/features/shared/filter-utils";
import {
  columns,
  type TranscriptomicRegistryRecord,
} from "@/features/transc-assess/columns";

export default function Page() {
  const [taxonInput, setTaxonInput] = useState("");
  const [records, setRecords] = useState<TranscriptomicRegistryRecord[]>([]);
  const [loading, setLoading] = useState(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);

  const csvContent = useMemo(() => {
    if (!records.length) {
      return "";
    }

    const headers = [
      "taxon_id",
      "transc_assess_date",
      "run_accession",
      "sample_tissue",
      "tissue_prediction",
      "qc_status",
      "uniquely_mapped_reads_percentage",
      "percentage_reads_mapped_to_multiple_loci",
      "percentage_reads_unmapped_too_short",
    ];

    const escapeCsvValue = (value: string | number | null | undefined) => {
      if (value === null || value === undefined) {
        return "";
      }

      const normalizedValue = String(value).replace(/"/g, '""');
      return /[",\n]/.test(normalizedValue)
        ? `"${normalizedValue}"`
        : normalizedValue;
    };

    const rows = records.map((record) =>
      headers
        .map((header) =>
          escapeCsvValue(record[header as keyof TranscriptomicRegistryRecord]),
        )
        .join(","),
    );

    return [headers.join(","), ...rows].join("\n");
  }, [records]);

  const handleDownload = (
    content: string,
    filename: string,
    type: string = "text/plain",
  ) => {
    if (!content) {
      setErrorMessage(`No ${filename} data available to download.`);
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
      setErrorMessage(
        `Error downloading ${filename}: ${
          error instanceof Error ? error.message : String(error)
        }`,
      );
    }
  };

  const handleGetRegistryData = async (inputValue: string = taxonInput) => {
    setErrorMessage(null);
    setRecords([]);

    const taxonIds = parseTaxonIds(inputValue);
    if (!taxonIds?.length) {
      setErrorMessage("Enter one or more taxonomy IDs separated by commas.");
      return;
    }

    setLoading(true);

    try {
      const res = await fetch("/api/transcriptomics/registry", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(cleanPayload({ taxon_ids: taxonIds })),
      });

      if (!res.ok) {
        const errorText = await res.text();
        setErrorMessage(`Failed to fetch transcriptomic registry data: ${errorText}`);
        return;
      }

      const result = await res.json();
      setRecords(result.records ?? []);
    } catch (error) {
      console.error("Error fetching transcriptomic registry data:", error);
      setErrorMessage(
        "Error fetching data: " + (error instanceof Error ? error.message : String(error)),
      );
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    const searchParams = new URLSearchParams(window.location.search);
    const taxonId = searchParams.get("taxon_id");
    if (!taxonId) {
      return;
    }

    setTaxonInput(taxonId);
    void handleGetRegistryData(taxonId);
  }, []);

  return (
    <div className="min-h-screen flex justify-center px-4 py-10 sm:px-6 lg:px-8">
      <div className="w-full max-w-6xl">
        <div className="rounded-2xl border-accent shadow-lg">
          <div className="rounded-2xl bg-secondary p-4 shadow-lg sm:p-8">
            <div className="mb-6">
              <h1 className="text-2xl font-bold">Transcriptomic registry</h1>
              <p className="mt-2 text-sm text-muted-foreground">
                Enter one or more taxonomy IDs to inspect transcriptomic registry runs and alignment metrics.
              </p>
            </div>

            <div className="grid grid-cols-1 gap-4 sm:grid-cols-[minmax(0,1fr)_auto] sm:items-end">
              <div>
                <Label htmlFor="taxon-id-input">Taxonomy ID(s)</Label>
                <Input
                  id="taxon-id-input"
                  type="text"
                  placeholder="9606, 10090"
                  className="mt-3 bg-filter-input-bg dark:bg-transparent"
                  value={taxonInput}
                  onChange={(event) => setTaxonInput(event.target.value)}
                />
              </div>
              <Button
                className="w-full sm:w-auto"
                size="lg"
                onClick={handleGetRegistryData}
                disabled={loading}
              >
                {loading ? (
                  <>
                    <Loader2 className="mr-2 animate-spin" />
                    Loading...
                  </>
                ) : (
                  <>
                    <Sparkles className="mr-2" />
                    Get registry data
                  </>
                )}
              </Button>
            </div>
          </div>
        </div>

        {errorMessage && (
          <Alert variant="destructive" className="mt-8">
            <Terminal />
            <AlertTitle>Heads up!</AlertTitle>
            <AlertDescription>{errorMessage}</AlertDescription>
          </Alert>
        )}

        {records.length > 0 && (
          <div className="mt-10 overflow-hidden rounded-2xl border shadow-lg">
            <div className="border-b px-4 py-6 sm:px-8">
              <div className="flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
              <h2 className="text-lg font-semibold">
                Registry records ({records.length})
              </h2>
                <Button
                  variant="outline"
                  onClick={() =>
                    handleDownload(
                      csvContent,
                      "transcriptomic_registry_records.csv",
                      "text/csv",
                    )
                  }
                >
                  <DownloadIcon className="mr-2 h-4 w-4" />
                  Download CSV
                </Button>
              </div>
            </div>
            <div className="overflow-x-auto">
              <div className="min-w-[1100px]">
                <DataTable columns={columns} data={records} />
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
