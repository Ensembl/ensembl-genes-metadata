"use client";

import { type ColumnDef, type HeaderContext } from "@tanstack/react-table";
import { ArrowUpDown } from "lucide-react";
import { Button } from "@/components/ui/button";

export type ProjectGCA = {
  gca: string;
  lowest_taxon_id: number;
  scientific_name: string;
  asm_level: string;
  infra_name: string;
  gb_status: string;
  genebuilder: string;
};

function sortableHeader(label: string) {
  return function SortableHeader({ column }: HeaderContext<ProjectGCA, unknown>) {
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

export const columns: ColumnDef<ProjectGCA>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA"),
    filterFn: "includesString",
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Taxon ID"),
    filterFn: (row, columnId, filterValue) =>
      String(row.getValue(columnId)).includes(String(filterValue).trim()),
  },
  {
    accessorKey: "scientific_name",
    header: sortableHeader("Scientific Name"),
    filterFn: "includesString",
  },
  {
    accessorKey: "asm_level",
    header: sortableHeader("ASM Level"),
    filterFn: "includesString",
  },
    {
    accessorKey: "infra_name",
    header: sortableHeader("Infra name"),
    filterFn: "includesString",
  },
  {
    accessorKey: "gb_status",
    header: sortableHeader("Status"),
    filterFn: "includesString",
  },
  {
    accessorKey: "genebuilder",
    header: sortableHeader("Genebuilder"),
    filterFn: "includesString",
  },
];
