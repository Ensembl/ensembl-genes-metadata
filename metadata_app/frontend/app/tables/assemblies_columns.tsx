"use client"

import { ColumnDef } from "@tanstack/react-table"
import { ArrowUpDown } from "lucide-react"
import { Button } from "@/components/ui/button"

export type Assemblies = {
  id: number
  bioproject_id: string
  associated_project: string
  gca: string
  scientific_name: string
  release_date: string
  lowest_taxon_id: number
  internal_clade: string
  asm_level: string
  asm_type: string
  asm_name: string
  refseq_accession: string
  is_current: string
  contig_n50: number
  total_sequence_length: number
}

function sortableHeader(label: string, accessor: string) {
  return ({ column }: { column: any }) => (
    <Button
      variant="ghost"
      className="p-0"
      onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
    >
      {label}
      <ArrowUpDown className="ml-2 h-4 w-4" />
    </Button>
  )
}

export const columns: ColumnDef<Assemblies>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA", "gca"),
  },
  {
    accessorKey: "bioproject_id",
    header: sortableHeader("BioProject ID", "bioproject_id"),
  },
  {
    accessorKey: "associated_project",
    header: sortableHeader("Associated Project", "associated_project"),
  },
  {
    accessorKey: "scientific_name",
    header: sortableHeader("Scientific Name", "scientific_name"),
  },
  {
    accessorKey: "release_date",
    header: sortableHeader("Release Date", "release_date"),
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Lowest Taxon ID", "lowest_taxon_id"),
  },
  {
    accessorKey: "internal_clade",
    header: sortableHeader("Internal Clade", "internal_clade"),
  },
  {
    accessorKey: "asm_level",
    header: sortableHeader("Assembly Level", "asm_level"),
  },
  {
    accessorKey: "asm_type",
    header: sortableHeader("Assembly Type", "asm_type"),
  },
  {
    accessorKey: "asm_name",
    header: sortableHeader("Assembly Name", "asm_name"),
  },
  {
    accessorKey: "refseq_accession",
    header: sortableHeader("RefSeq Accession", "refseq_accession"),
  },
  {
    accessorKey: "is_current",
    header: sortableHeader("Current", "is_current"),
  },
  {
    accessorKey: "contig_n50",
    header: sortableHeader("Contig N50", "contig_n50"),
  },
  {
    accessorKey: "total_sequence_length",
    header: sortableHeader("Sequence Length", "total_sequence_length"),
  },
]