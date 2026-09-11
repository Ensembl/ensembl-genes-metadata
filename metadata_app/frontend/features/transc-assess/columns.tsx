"use client";

import { type ColumnDef, type HeaderContext } from "@tanstack/react-table";
import { ArrowUpDown } from "lucide-react";

import { Button } from "@/components/ui/button";

export type TranscriptomicRegistryRecord = {
  taxon_id: number;
  transc_assess_date?: string | null;
  run_accession?: string | null;
  sample_tissue?: string | null;
  tissue_prediction?: string | null;
  qc_status?: string | null;
  uniquely_mapped_reads_percentage?: number | null;
  percentage_reads_mapped_to_multiple_loci?: number | null;
  percentage_reads_unmapped_too_short?: number | null;
};

const stringValue = (value: string | null | undefined, fallback = "Unknown") =>
  value && value.trim() ? value : fallback;

function sortableHeader(label: string) {
  return function SortableHeader({
    column,
  }: HeaderContext<TranscriptomicRegistryRecord, unknown>) {
    return (
      <Button
        variant="ghost"
        className="cursor-pointer hover:bg-transparent hover:text-inherit"
        onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
      >
        {label}
        <ArrowUpDown className="ml-2 h-4 w-4" />
      </Button>
    );
  };
}

export const columns: ColumnDef<TranscriptomicRegistryRecord>[] = [
  {
    accessorKey: "taxon_id",
    header: sortableHeader("Taxon ID"),
  },
  {
    accessorKey: "transc_assess_date",
    header: sortableHeader("Last check"),
    cell: ({ row }) => row.getValue("transc_assess_date") ?? "Unknown",
  },
  {
    accessorKey: "run_accession",
    header: sortableHeader("Run accession"),
    cell: ({ row }) => stringValue(row.getValue("run_accession")),
  },
  {
    accessorKey: "sample_tissue",
    header: sortableHeader("Sample tissue"),
    cell: ({ row }) => stringValue(row.getValue("sample_tissue")),
  },
  {
    accessorKey: "tissue_prediction",
    header: sortableHeader("Tissue prediction"),
    cell: ({ row }) => stringValue(row.getValue("tissue_prediction")),
  },
  {
    accessorKey: "qc_status",
    header: sortableHeader("QC status"),
    cell: ({ row }) => stringValue(row.getValue("qc_status")),
  },
  {
    accessorKey: "uniquely_mapped_reads_percentage",
    header: sortableHeader("Unique mapped %"),
    cell: ({ row }) =>
      row.getValue("uniquely_mapped_reads_percentage") ?? "Unknown",
  },
  {
    accessorKey: "percentage_reads_mapped_to_multiple_loci",
    header: sortableHeader("Multiple loci %"),
    cell: ({ row }) =>
      row.getValue("percentage_reads_mapped_to_multiple_loci") ?? "Unknown",
  },
  {
    accessorKey: "percentage_reads_unmapped_too_short",
    header: sortableHeader("Unmapped too short %"),
    cell: ({ row }) =>
      row.getValue("percentage_reads_unmapped_too_short") ?? "Unknown",
  },
];
