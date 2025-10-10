"use client"

import { ColumnDef } from "@tanstack/react-table"
import { ArrowUpDown } from "lucide-react"
import { Button } from "@/components/ui/button"
import { BadgeCheck, Badge } from "lucide-react";


export type Report = {
  id: number
  associated_project: string
  gca: string
  genebuilder: string
  gb_status: string
  ftp: string
  latest_annotated: string
  protein_busco: string
  release_date: string
  last_genebuild_update: string
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

export const columns: ColumnDef<Report>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA", "gca"),
  },
  {
    accessorKey: "associated_project",
    header: sortableHeader("Project", "associated_project"),
  },
  {
    accessorKey: "genebuilder",
    header: sortableHeader("Genebuilder", "genebuilder"),
  },
  {
    accessorKey: "gb_status",
    header: sortableHeader("Staus", "gb_status"),
  },

  {
    accessorKey: "ftp",
    header: sortableHeader("FTP", "ftp"),
  },
  {
  accessorKey: "last_genebuild_update",
  header: sortableHeader("Annotation Date", "last_genebuild_update"),
    cell: ({ row }) => {
      const fullDate = row.getValue("last_genebuild_update") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
  },
  {
  accessorKey: "release_date",
  header: sortableHeader("Release Date Beta", "release_date"),
    cell: ({ row }) => {
      const fullDate = row.getValue("release_date") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
  },
  {
    accessorKey: "protein_busco",
    header: sortableHeader("BUSCO", "protein_busco"),
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