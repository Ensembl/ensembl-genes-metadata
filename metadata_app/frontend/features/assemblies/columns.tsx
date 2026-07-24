"use client";

import { type ColumnDef, type HeaderContext } from "@tanstack/react-table";
import {ArrowUpDown, Badge, BadgeCheck} from "lucide-react";
import { Button } from "@/components/ui/button";

export type Assemblies = {
  id: number;
  bioproject_id: string;
  associated_project: string;
  gca: string;
  scientific_name: string;
  release_date: string;
  lowest_taxon_id: number;
  internal_clade: string;
  is_current: string;
  gb_status: string;
  other_version_live: string;
  "assembly.busco": string;
  "assembly.busco_dataset": string;
  short_read_paired_end_illumina_lowest: number;
  short_read_paired_end_illumina: number;
  lowest_aligned_count: number;
  long_read_pacbio: number;
};

function sortableHeader(label: string) {
  return function SortableHeader({ column }: HeaderContext<Assemblies, unknown>) {
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

export const columns: ColumnDef<Assemblies>[] = [
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
    accessorKey: "release_date",
    header: sortableHeader("Release Date"),
    cell: ({ row }) => dateCell(row.getValue("release_date")),
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Taxon ID"),
  },
  {
    accessorKey: "internal_clade",
    header: sortableHeader("Internal Clade"),
  },
  {
    accessorKey: "is_current",
    header: sortableHeader("Latest GCA"),
  },
    {
    accessorKey: "lowest_aligned_count",
    header: sortableHeader("Transcr. reg. lowest"),
  },
  {
    accessorKey: "short_read_paired_end_illumina_lowest",
    header: sortableHeader("RNA lowest"),
  },
  {
    accessorKey: "short_read_paired_end_illumina",
    header: sortableHeader("RNA genus"),
  },
    {
    accessorKey: "long_read_pacbio",
    header: sortableHeader("PacBio genus"),
  },
    {
    accessorKey: "gb_status",
    header: sortableHeader("Status"),
  },
    {
  accessorKey: "other_version_live",
  header: "Other version live",
  cell: ({ row }) => (
    <div className="flex justify-center items-center w-full">
      {row.getValue("other_version_live") === "yes" ? (
        <BadgeCheck className="h-5 w-5 text-foreground" />
      ) : (
        <Badge className="h-5 w-5 text-foreground" />
      )}
    </div>
  ),
},
  {
    id: "assembly_busco",
    accessorFn: (row) => row["assembly.busco"],
    header: sortableHeader("gBUSCO"),
  },
  {
    id: "assembly_busco_dataset",
    accessorFn: (row) => row["assembly.busco_dataset"],
    header: sortableHeader(" gBUSCO lineage"),
  },

];
