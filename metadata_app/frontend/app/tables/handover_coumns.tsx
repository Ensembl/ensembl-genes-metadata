"use client";

import { ColumnDef } from "@tanstack/react-table";
import { ArrowUpDown } from "lucide-react";
import { Button } from "@/components/ui/button";

export type Handover = {
  gb_status: string;
  date_status_update: string;
  gca: string;
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

export const columns: ColumnDef<Handover>[] = [
  {
    accessorKey: "gca",
    header: sortableHeader("GCA"),
  },
  {
    accessorKey: "date_status_update",
    header: sortableHeader("Date Status Update"),
  },
    {
    accessorKey: "bioproject_name",
    header: sortableHeader("Bioproject"),
  },
  {
    accessorKey: "gb_status",
    header: sortableHeader("GB Status"),
  },

];