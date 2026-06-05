"use client"

import * as React from "react"
import {
  ColumnDef,
  ColumnFiltersState,
  SortingState,
  VisibilityState,
  flexRender,
  getCoreRowModel,
  getFilteredRowModel,
  getSortedRowModel,
  useReactTable,
} from "@tanstack/react-table"
import { ArrowUpDown } from "lucide-react"

import { Button } from "@/components/ui/button"
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"

export type Report_Table = {
  id: string
  project_name: string
  total_assemblies: number
  annotation_candidates: number
  unannotated: number
  unannotated_main: number
  unannotated_anno: number
  in_progress: number
  live: number
}

export const columns: ColumnDef<Report_Table>[] = [
  {
    accessorKey: "project_name",
    header: ({ column }) => (
      <Button
        variant="ghost"
        onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
      >
        Project Name
        <ArrowUpDown />
      </Button>
    ),
    cell: ({ row }) => <div>{row.getValue("project_name")}</div>,
  },

  {
    accessorKey: "total_assemblies",
    header: "Total assemblies",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("total_assemblies")}</div>
    ),
  },
    {
    accessorKey: "annotation_candidates",
    header: () => ( <> Annotation candidates <br /> total </>),
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("annotation_candidates")}</div>
    ),
  },
  {
    accessorKey: "unannotated",
    header: () => ( <> Annotation candidates <br /> unannotated total</>),
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("unannotated")}</div>
    ),
  },
  {
    accessorKey: "unannotated_main",
    header: () => ( <> Annotation candidates <br /> unannotated main</>),
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("unannotated_main")}</div>
    ),
  },
  {
    accessorKey: "unannotated_anno",
    header: () => ( <> Annotation candidates <br /> unannotated anno</>),
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("unannotated_anno")}</div>
    ),
  },
  {
    accessorKey: "in_progress",
    header: "In progress",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("in_progress")}</div>
    ),
  },
    {
    accessorKey: "live",
    header: "Live",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("live")}</div>
    ),
  },
]

export function CardsReportTable() {
  const [data, setData] = React.useState<Report_Table[]>([])
  const [loading, setLoading] = React.useState(true)
  const [sorting, setSorting] = React.useState<SortingState>([])
  const [columnFilters, setColumnFilters] = React.useState<ColumnFiltersState>([])
  const [columnVisibility, setColumnVisibility] = React.useState<VisibilityState>({})
  const [rowSelection, setRowSelection] = React.useState({})

  React.useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await fetch("/api/report/asm/report/asm/main")
        const json = await res.json()
        console.log("API response:", json)
        const formatted = json.map((item: Report_Table) => ({
          id: item.project_name,
          project_name: item.project_name || "Unknown",
          total_assemblies: item.total_assemblies,
          unannotated: item.unannotated,
          unannotated_main: item.unannotated_main ?? 0,
          unannotated_anno: item.unannotated_anno ?? 0,
          annotation_candidates: item.annotation_candidates,
          in_progress: item.in_progress,
          live: item.live,

        }))

        setData(formatted)
      } catch (err) {
        console.error("Error fetching data:", err)
      } finally {
        setLoading(false)
      }
    }

    fetchData()
  }, [])

  const table = useReactTable({
    data,
    columns,
    onSortingChange: setSorting,
    onColumnFiltersChange: setColumnFilters,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
    onColumnVisibilityChange: setColumnVisibility,
    onRowSelectionChange: setRowSelection,
    state: {
      sorting,
      columnFilters,
      columnVisibility,
      rowSelection,
    },
  })


  return (
    <Card className="dark:bg-secondary">
      <CardHeader className="mb-4">
        <CardTitle className="text-xl">Biodiversity projects overview</CardTitle>
        <CardDescription>Number of annotations per project. Qualified assemblies shows the number of chromosome level, primary, current assemblies with transciptomic data.
        In progress shows in progress, completed and handed over annotations. This takes a couple of seconds to load.</CardDescription>
      </CardHeader>
      <CardContent className="overflow-visible mb-4">
        <div className="rounded-md">
          {loading ? (
            <p className="text-muted-foreground">Loading data...</p>
          ) : (
            <Table>
              <TableHeader>
                {table.getHeaderGroups().map((headerGroup) => (
                  <TableRow key={headerGroup.id}>
                    {headerGroup.headers.map((header) => (
                      <TableHead
                        key={header.id}
                        className="[&:has([role=checkbox])]:pl-3 text-center"
                      >
                        {header.isPlaceholder
                          ? null
                          : flexRender(header.column.columnDef.header, header.getContext())}
                      </TableHead>
                    ))}
                  </TableRow>
                ))}
              </TableHeader>
              <TableBody>
                {table.getRowModel().rows.length ? (
                  table.getRowModel().rows.map((row) => (
                    <TableRow
                      key={row.id}
                      data-state={row.getIsSelected() && "selected"}
                    >
                      {row.getVisibleCells().map((cell) => (
                        <TableCell
                          key={cell.id}
                          className="[&:has([role=checkbox])]:pl-3 text-center"
                        >
                          {flexRender(cell.column.columnDef.cell, cell.getContext())}
                        </TableCell>
                      ))}
                    </TableRow>
                  ))
                ) : (
                  <TableRow>
                    <TableCell
                      colSpan={columns.length}
                      className="h-24 text-center"
                    >
                      No results.
                    </TableCell>
                  </TableRow>
                )}
              </TableBody>
            </Table>
          )}
        </div>
      </CardContent>
    </Card>
  )
}
