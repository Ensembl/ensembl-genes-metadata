"use client";

import { ColumnDef } from "@tanstack/react-table";
import { ArrowUpDown } from "lucide-react";
import { Button } from "@/components/ui/button";

export type ProjectGCA = {
  gca: string;
  lowest_taxon_id: number;
  scientific_name: string;
  asm_level: string;
  gb_status: string;
  genebuilder: string;
};


function sortableHeader(label: string) {
  return ({ column }: { column: any }) => (
    <Button
      variant="ghost"
      className="hover:bg-transparent hover:text-inherit cursor-pointer"
      onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
    >
      {label}
      <ArrowUpDown className="ml-2 h-4 w-4" />
    </Button>
  );
}

export const columns: ColumnDef<ProjectGCA>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA"),
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Taxon ID"),
  },
  {
    accessorKey: "scientific_name",
    header: sortableHeader("Scientific Name"),
  },

  {
    accessorKey: "asm_level",
    header: sortableHeader("ASM Level"),
  },
  {
    accessorKey: "gb_status",
    header: sortableHeader("Status"),
  },
  {
    accessorKey: "genebuilder",
    header: sortableHeader("Genebuilder"),
  },
];