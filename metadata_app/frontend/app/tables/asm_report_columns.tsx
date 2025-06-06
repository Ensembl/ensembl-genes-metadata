"use client"

import { ColumnDef } from "@tanstack/react-table"
import { ArrowUpDown } from "lucide-react"
import { Button } from "@/components/ui/button"
import { BadgeCheck, Badge } from "lucide-react";


export type Report = {
  id: number
  associated_project: string
  gca: string
  scientific_name: string
  release_date: string
  lowest_taxon_id: number
  genus_taxon_id: number
  transcriptomic_evidence: string
  internal_clade: string
  asm_type: string
  asm_name: string
  refseq_accession: string
  asm_level:string
  contig_n50: number
  total_sequence_length: number
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
    accessorKey: "scientific_name",
    header: sortableHeader("Species", "scientific_name"),
  },
  {
    accessorKey: "lowest_taxon_id",
    header: sortableHeader("Lowest Taxon ID", "lowest_taxon_id"),
  },
  {
    accessorKey: "genus_taxon_id",
    header: sortableHeader("Genus Taxon ID", "genus_taxon_id"),
  },
  {
    accessorKey: "asm_name",
    header: sortableHeader("Assembly name", "asm_name"),
  },
  {
    accessorKey: "internal_clade",
    header: sortableHeader("Clade", "internal_clade"),
  },
  {
  accessorKey: "release_date",
  header: sortableHeader("Release Date", "release_date"),
    cell: ({ row }) => {
      const fullDate = row.getValue("release_date") as string;
      const dateOnly = fullDate.split("T")[0]; // or use new Date(fullDate).toISOString().split("T")[0]
      return dateOnly;
    },
  },
  {
    accessorKey: "asm_type",
    header: sortableHeader("Type", "asm_type"),
  },
    {
    accessorKey: "refseq_accession",
    header: sortableHeader("RefSeq", "refseq_accession"),
  },
{
    accessorKey: "asm_level",
    header: sortableHeader("Level", "asm_level"),
  },
{
    accessorKey: "contig_n50",
    header: sortableHeader("Contig N50", "contig_n50"),
  },
{
    accessorKey: "total_sequence_length",
    header: sortableHeader("Length", "total_sequence_length"),
  },

{
  accessorKey: "transcriptomic_evidence",
  header: sortableHeader("RNA", "transcriptomic_evidence"),
  cell: ({ row }) => {
    const isLatest = row.getValue("transcriptomic_evidence") === "yes";
    return isLatest ? (
      <BadgeCheck className="text-foreground w-5 h-5"  />
    ) : (
      <Badge className="text-foreground w-5 h-5"  />
    );
  },
}
]