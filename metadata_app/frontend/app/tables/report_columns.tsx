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
  release_type: string
  ftp: string
  latest_annotated: string
  busco_protein: string
  date_completed_beta: string
  release_date_beta: string
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
    accessorKey: "release_type",
    header: sortableHeader("Site", "release_type"),
  },
  {
    accessorKey: "ftp",
    header: sortableHeader("FTP", "ftp"),
  },
  {
  accessorKey: "date_completed_beta",
  header: sortableHeader("Annotation Date", "date_completed_beta"),
    cell: ({ row }) => {
      const fullDate = row.getValue("date_completed_beta") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
  },
  {
  accessorKey: "release_date_beta",
  header: sortableHeader("Release Date Beta", "release_date_beta"),
    cell: ({ row }) => {
      const fullDate = row.getValue("release_date_beta") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
  },
  {
    accessorKey: "busco_protein",
    header: sortableHeader("BUSCO", "busco_protein"),
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