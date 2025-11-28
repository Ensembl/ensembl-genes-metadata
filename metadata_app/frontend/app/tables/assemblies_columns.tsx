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
  is_current: string,
  short_read_paired_end_illumina_lowest: number,
  short_read_paired_end_illumina:number,
  lowest_aligned_count: number,
}

function sortableHeader(label: string, accessor: string) {
  return ({ column }: { column: any }) => (
    <Button
      variant="ghost"
      className="hover:bg-transparent hover:text-inherit cursor-pointer"
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
    cell: ({ row }) => {
      const fullDate = row.getValue("release_date") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
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
    accessorKey: "is_current",
    header: sortableHeader("Latest GCA", "is_current"),
  },
    {
    accessorKey: "lowest_aligned_count",
    header: sortableHeader("Transcr. reg. lowest", "lowest_aligned_count"),
  },
    {
    accessorKey: "short_read_paired_end_illumina_lowest",
    header: sortableHeader("RNA lowest", "short_read_paired_end_illumina_lowest"),
  },
    {
    accessorKey: "short_read_paired_end_illumina",
    header: sortableHeader("RNA genus", "short_read_paired_end_illumina"),
  },
]