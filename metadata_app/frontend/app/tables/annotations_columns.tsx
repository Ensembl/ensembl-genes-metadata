"use client"

import { ColumnDef } from "@tanstack/react-table"
import { ArrowUpDown } from "lucide-react"
import { Button } from "@/components/ui/button"

export type Annotations = {
  id: number
  bioproject_id: string
  associated_project: string
  gca: string
  scientific_name: string
  annotation_date: string
  lowest_taxon_id: number
  gb_status: string
  release_site: string
  latest_annotated: string
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

export const columns: ColumnDef<Annotations>[] = [
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
  accessorKey: "annotation_date",
  header: sortableHeader("Annotation Date", "annotation_date"),
    cell: ({ row }) => {
      const fullDate = row.getValue("annotation_date") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Lowest Taxon ID", "lowest_taxon_id"),
  },
  {
    accessorKey: "release_site",
    header: sortableHeader("Release site", "release_site"),
  },
  {
    accessorKey: "latest_annotated",
    header: sortableHeader("Latest GCA annotated", "latest_annotated"),
  }
]