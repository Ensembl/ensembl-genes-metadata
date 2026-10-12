"use client";

import { type ColumnDef, type HeaderContext } from "@tanstack/react-table";
import { ArrowUpDown, Badge, BadgeCheck } from "lucide-react";
import { Button } from "@/components/ui/button";

export type Report = {
  id: number;
  associated_project: string;
  gca: string;
  scientific_name: string;
  release_date: string;
  lowest_taxon_id: number;
  genus_taxon_id: number;
  transcriptomic_evidence: string;
  internal_clade: string;
  asm_type: string;
  asm_name: string;
  refseq_accession: string;
  asm_level: string;
  contig_n50: number;
  total_sequence_length: number;
};

function sortableHeader(label: string) {
  return function SortableHeader({ column }: HeaderContext<Report, unknown>) {
    return (
      <Button
        variant="ghost"
        className="hover:bg-transparent hover:text-inherit cursor-pointer"
        onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
      >
        {label}
        <ArrowUpDown className="ml-2 h-4 w-4" />
      </Button>
    );
  };
}

function dateCell(value: unknown) {
  return typeof value === "string" ? value.split("T")[0] : "";
}

export const columns: ColumnDef<Report>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA"),
  },
  {
    accessorKey: "scientific_name",
    header: sortableHeader("Species"),
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Lowest Taxon ID"),
  },
  {
    accessorKey: "genus_taxon_id",
    header: sortableHeader("Genus Taxon ID"),
  },
  {
    accessorKey: "asm_name",
    header: sortableHeader("Assembly name"),
  },
  {
    accessorKey: "internal_clade",
    header: sortableHeader("Clade"),
  },
  {
    accessorKey: "release_date",
    header: sortableHeader("Release Date"),
    cell: ({ row }) => dateCell(row.getValue("release_date")),
  },
  {
    accessorKey: "asm_type",
    header: sortableHeader("Type"),
  },
  {
    accessorKey: "refseq_accession",
    header: sortableHeader("RefSeq"),
  },
  {
    accessorKey: "asm_level",
    header: sortableHeader("Level"),
  },
  {
    accessorKey: "contig_n50",
    header: sortableHeader("Contig N50"),
  },
  {
    accessorKey: "total_sequence_length",
    header: sortableHeader("Length"),
  },
  {
    accessorKey: "transcriptomic_evidence",
    header: sortableHeader("RNA"),
    cell: ({ row }) =>
      row.getValue("transcriptomic_evidence") === "yes" ? (
        <BadgeCheck className="text-foreground w-5 h-5" />
      ) : (
        <Badge className="text-foreground w-5 h-5" />
      ),
  },
];
