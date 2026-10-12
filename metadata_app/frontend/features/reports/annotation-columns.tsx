"use client";

import { type ColumnDef, type HeaderContext } from "@tanstack/react-table";
import { ArrowUpDown, Badge, BadgeCheck } from "lucide-react";
import { Button } from "@/components/ui/button";

export type Report = {
  id: number;
  associated_project: string;
  infra_name: string;
  gca: string;
  genebuilder: string;
  gb_status: string;
  ftp: string;
  latest_annotated: string;
  protein_busco: string;
  release_date: string;
  last_genebuild_update: string;
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
    accessorKey: "associated_project",
    header: sortableHeader("Project"),
  },
  {
    accessorKey: "infra_name",
    header: sortableHeader("Infra Name"),
  },
  {
    accessorKey: "genebuilder",
    header: sortableHeader("Genebuilder"),
  },
  {
    accessorKey: "gb_status",
    header: sortableHeader("Status"),
  },

  {
    accessorKey: "last_genebuild_update",
    header: sortableHeader("Annotation Date"),
    cell: ({ row }) => dateCell(row.getValue("last_genebuild_update")),
  },
  {
    accessorKey: "release_date",
    header: sortableHeader("Release Date"),
    cell: ({ row }) => dateCell(row.getValue("release_date")),
  },
  {
    accessorKey: "protein_busco",
    header: sortableHeader("BUSCO"),
  },
  {
    accessorKey: "latest_annotated",
    header: sortableHeader("Latest Annotated"),
    cell: ({ row }) =>
      row.getValue("latest_annotated") === "Yes" ? (
        <BadgeCheck className="text-foreground w-5 h-5" />
      ) : (
        <Badge className="text-foreground w-5 h-5" />
      ),
  },
  {
    accessorKey: "ftp",
    header: sortableHeader("FTP"),
  },
];
