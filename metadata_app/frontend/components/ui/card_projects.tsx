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
  getPaginationRowModel,
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
import { Input } from "@/components/ui/input"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"

export type Project = {
  id: string
  bioproject_id: string
  bioproject_name: string
  annotation_count: number
  qualified_assembly_count: number
  in_progress: number
}

export const columns: ColumnDef<Project>[] = [
  {
    accessorKey: "bioproject_name",
    header: ({ column }) => (
      <Button
        variant="ghost"
        onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
      >
        Name
        <ArrowUpDown />
      </Button>
    ),
    cell: ({ row }) => <div>{row.getValue("bioproject_name")}</div>,
  },
  {
    accessorKey: "annotation_count",
    header: "Live annotations",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("annotation_count")}</div>
    ),
  },
    {
    accessorKey: "in_progress",
    header: "In-progress",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("in_progress")}</div>
    ),
  },
  {
    accessorKey: "qualified_assembly_count",
    header: "Unannotated",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("qualified_assembly_count")}</div>
    ),
  },
]

export function CardsDataTable() {
  const [data, setData] = React.useState<Project[]>([])
  const [loading, setLoading] = React.useState(true)
  const [sorting, setSorting] = React.useState<SortingState>([])
  const [columnFilters, setColumnFilters] = React.useState<ColumnFiltersState>([])
  const [columnVisibility, setColumnVisibility] = React.useState<VisibilityState>({})
  const [rowSelection, setRowSelection] = React.useState({})

  React.useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await fetch("http://127.0.0.1:8000/api/home_page/home/bioproject")
        const json = await res.json()
        console.log("API response:", json)
        const formatted = json.map((item: Project) => ({
          id: item.bioproject_id,
          annotation_count: item.annotation_count,
          bioproject_id: item.bioproject_id,
          bioproject_name: item.bioproject_name || "Unknown",
          qualified_assembly_count: item.qualified_assembly_count,
          in_progress: item.in_progress,
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
    getPaginationRowModel: getPaginationRowModel(),
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
        <CardTitle className="text-xl -mb-2">Biodiversity projects</CardTitle>
        <CardDescription>Number of annotations per project. Unannotated shows the number of chromosome level, primary, current assemblies.</CardDescription>
      </CardHeader>
      <CardContent>
        <div className="mb-4 flex items-center gap-4">
          <Input
            placeholder="Type project name"
            value={(table.getColumn("bioproject_name")?.getFilterValue() as string) ?? ""}
            onChange={(event) =>
              table.getColumn("bioproject_name")?.setFilterValue(event.target.value)
            }
            className="max-w-sm"
          />
        </div>
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
                        className="[&:has([role=checkbox])]:pl-3"
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
                          className="[&:has([role=checkbox])]:pl-3"
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