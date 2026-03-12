"use client"

import { ColumnDef } from "@tanstack/react-table"
import { ArrowUpDown } from "lucide-react"
import { Button } from "@/components/ui/button"
import { BadgeCheck, Badge } from "lucide-react";


export type Annotations = {
  id: number
  bioproject_id: string
  associated_project: string
  gca: string
  scientific_name: string
  lowest_taxon_id: number
  gb_status: string
  latest_annotated: string
  date_status_update: string
  release_date: string
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
    accessorKey: "gb_status",
    header: sortableHeader("Annotation status", "gb_status"),
  },
  {
  accessorKey: "release_date",
  header: sortableHeader("Release Date", "release_date"),
    cell: ({ row }) => {
  const fullDate = row.getValue("release_date");
  return typeof fullDate === "string"
    ? fullDate.split("T")[0]
    : "";
},
  },
    {
  accessorKey: "date_status_update",
  header: sortableHeader("Status Update", "date_status_update"),
    cell: ({ row }) => {
  const fullDate = row.getValue("date_status_update");
  return typeof fullDate === "string"
    ? fullDate.split("T")[0]
    : "";
}, },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Lowest Taxon ID", "lowest_taxon_id"),
  },
{
  accessorKey: "latest_annotated",
  header: sortableHeader("Latest GCA Annotated", "latest_annotated"),
  cell: ({ row }) => {
    const isLatest = row.getValue("latest_annotated") === "Yes";
    return isLatest ? (
      <BadgeCheck className="text-foreground w-5 h-5"  />
    ) : (
      <Badge className="text-foreground w-5 h-5"  />
    );
  },
}
]