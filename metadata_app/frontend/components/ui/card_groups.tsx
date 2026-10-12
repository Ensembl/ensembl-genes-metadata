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
  group_id: string
  group_name: string
  annotation_count: number
  qualified_assembly_count: number
  in_progress: number
}

export const columns: ColumnDef<Project>[] = [
  {
    accessorKey: "group_name",
    header: ({ column }) => (
      <Button
        variant="ghost"
        onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
      >
        Name
        <ArrowUpDown />
      </Button>
    ),
    cell: ({ row }) => <div>{row.getValue("group_name")}</div>,
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

export function CardsDataTableGroup() {
  const [data, setData] = React.useState<Project[]>([])
  const [loading, setLoading] = React.useState(true)
  const [sorting, setSorting] = React.useState<SortingState>([])
  const [columnFilters, setColumnFilters] = React.useState<ColumnFiltersState>([])
  const [columnVisibility, setColumnVisibility] = React.useState<VisibilityState>({})
  const [rowSelection, setRowSelection] = React.useState({})

  React.useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await fetch("/api/home_page/home/group")
        const json = await res.json()
        console.log("API response:", json)
        const formatted = json.map((item: Project) => ({
          id: item.group_id,
          annotation_count: item.annotation_count,
          group_id: item.group_id,
          group_name: item.group_name || "Unknown",
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
        <CardTitle className="text-lg -mb-2">Other projects</CardTitle>
        <CardDescription>Number of annotations per project. Unannotated shows the number of chromosome level, primary, current assemblies.</CardDescription>
      </CardHeader>
      <CardContent className="overflow-visible">
        <div className="mb-4 flex items-center gap-4">
          <Input
            placeholder="Type group name"
            value={(table.getColumn("group_name")?.getFilterValue() as string) ?? ""}
            onChange={(event) =>
              table.getColumn("group_name")?.setFilterValue(event.target.value)
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