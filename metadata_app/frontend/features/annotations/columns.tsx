"use client";

import { type ColumnDef, type HeaderContext } from "@tanstack/react-table";
import { ArrowUpDown, Badge, BadgeCheck } from "lucide-react";
import { Button } from "@/components/ui/button";

export type Annotations = {
  id: number;
  bioproject_id: string;
  associated_project: string;
  gca: string;
  scientific_name: string;
  lowest_taxon_id: number;
  gb_status: string;
  latest_annotated: string;
  date_status_update: string;
  release_date: string;
  assembly_busco: string;
  assembly_busco_lineage: string;
  protein_busco: string;
  protein_busco_lineage: string;
};

function sortableHeader(label: string) {
  return function SortableHeader({ column }: HeaderContext<Annotations, unknown>) {
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

export const columns: ColumnDef<Annotations>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA"),
  },
  {
    accessorKey: "bioproject_id",
    header: sortableHeader("BioProject ID"),
  },
  {
    accessorKey: "associated_project",
    header: sortableHeader("Associated Project"),
  },
  {
    accessorKey: "scientific_name",
    header: sortableHeader("Scientific Name"),
  },
  {
    accessorKey: "gb_status",
    header: sortableHeader("Annotation status"),
  },
  {
    accessorKey: "release_date",
    header: sortableHeader("Release Date"),
    cell: ({ row }) => dateCell(row.getValue("release_date")),
  },
  {
    accessorKey: "date_status_update",
    header: sortableHeader("Status Update"),
    cell: ({ row }) => dateCell(row.getValue("date_status_update")),
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Lowest Taxon ID"),
  },
  {
    accessorKey: "latest_annotated",
    header: sortableHeader("Latest GCA Annotated"),
    cell: ({ row }) =>
      row.getValue("latest_annotated") === "Yes" ? (
        <BadgeCheck className="text-foreground w-5 h-5" />
      ) : (
        <Badge className="text-foreground w-5 h-5" />
      ),
  },
  {
    accessorKey: "assembly_busco",
    header: sortableHeader("Assembly BUSCO"),
  },
  {
    accessorKey: "assembly_busco_lineage",
    header: sortableHeader("Assembly BUSCO lineage"),
  },
  {
    accessorKey: "protein_busco",
    header: sortableHeader("Protein BUSCO"),
  },
  {
    accessorKey: "protein_busco_lineage",
    header: sortableHeader("Protein BUSCO lineage"),
  },
];
