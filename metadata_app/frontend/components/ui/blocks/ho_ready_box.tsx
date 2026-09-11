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
import { Checkbox } from '@/components/ui/checkbox'
import {
  Dialog,
  DialogClose,
  DialogContent,
  DialogDescription,
  DialogFooter,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
} from "@/components/ui/dialog"
import { Button } from "@/components/ui/button"
import {
    Card, CardAction,
    CardContent,
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
import {
  Select,
  SelectContent,
  SelectGroup,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select"

export type Handover = {
  id: string
  bioproject_name: string
  gca: string
    scientific_name:string
  date_status_update: string
  annotation_method: string
}

export const columns: ColumnDef<Handover>[] = [
  {
    id: 'select',
    header: ({ table }) => (
      <Checkbox
        checked={table.getIsAllPageRowsSelected() || (table.getIsSomePageRowsSelected() && 'indeterminate')}
        onCheckedChange={value => table.toggleAllPageRowsSelected(!!value)}
        aria-label='Select all'
      />
    ),
    cell: ({ row }) => (
      <Checkbox
        checked={row.getIsSelected()}
        onCheckedChange={value => row.toggleSelected(!!value)}
        aria-label='Select row'
      />
    ),
    enableSorting: false,
    enableHiding: false
  },
    {
    accessorKey: "bioproject_name",
    header: ({ column }) => (

        <Button
        variant="ghost"
        onClick={() => column.toggleSorting(column.getIsSorted() === "asc")}
      >
        BioProject
        <ArrowUpDown />
      </Button>
    ),
    cell: ({ row }) => <div>{row.getValue("bioproject_name")}</div>,
  },

  {
    accessorKey: "gca",
    header: "Accession",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("gca")}</div>
    ),
  },
    {
    accessorKey: "scientific_name",
    header: "Scientific name",
    cell: ({ row }) => (
      <div>{row.getValue("scientific_name")}</div>
    ),
  },
    {
    accessorKey: "date_status_update",
    header: "Last status update",
    cell: ({ row }) => (
      <div className="capitalize">{row.getValue("date_status_update")}</div>
    ),
  },
    {
    accessorKey: "annotation_method",
    header: () => ( <> Method </>),
    cell: ({ row }) => (
      <div>{row.getValue("annotation_method")}</div>
    ),
  },
]

type HoReadyBoxProps = {
  data: Handover[];
  loading: boolean;
  genebuilder: string | null;
};

export function HoReadyBox({ data, loading, genebuilder }: HoReadyBoxProps) {
  const [sorting, setSorting] = React.useState<SortingState>([])
  const [columnFilters, setColumnFilters] = React.useState<ColumnFiltersState>([])
  const [columnVisibility, setColumnVisibility] = React.useState<VisibilityState>({})
  const [rowSelection, setRowSelection] = React.useState({})
    const [tableData, setTableData] = React.useState<Handover[]>(data)
    const [selectedStatus, setSelectedStatus] = React.useState<string | undefined>();

    React.useEffect(() => {
  setTableData(data)
}, [data])

  const table = useReactTable({
    data: tableData,
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
    const [isDialogOpen, setIsDialogOpen] = React.useState(false)

    const handleChangeStatus = () => {
  const selectedRows = table.getSelectedRowModel().rows
  if (!selectedRows.length) return

  const selectedData = selectedRows.map(row => row.original)
  console.log("Selected rows:", selectedData)

}

const reloadData = async () => {
  if (!genebuilder) return
  try {
    const response = await fetch(`/api/handover/handover/genebuilder`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ genebuilder }),
    })
    if (!response.ok) throw new Error("Failed to fetch data")
    const json = await response.json()
    setTableData(json.df_ready)
  } catch (error) {
    console.error("Error reloading table data:", error)
  }
}


  return (
    <Card id="ho-ready-table" className="min-w-0 dark:bg-secondary">
      <CardHeader className="mb-4">
        <CardTitle className="text-xl">Handover ready cores</CardTitle>
          <CardAction>
              <Dialog open={isDialogOpen} onOpenChange={setIsDialogOpen}>
              <DialogTrigger asChild>
                      <Button
                        onClick={() => setIsDialogOpen(true)}
                        disabled={!Object.keys(rowSelection).length}
                        variant="default"
                      >
                        Change status
                      </Button>
                    </DialogTrigger>
                  <DialogContent className="sm:max-w-[425px]">
                  <DialogHeader>
                    <DialogTitle>Change status in the registry</DialogTitle>
                    <DialogDescription>
                      Click apply if you are certain
                    </DialogDescription>
                  </DialogHeader>

                      <div>
                        <p className="text-sm text-foreground">
                            Select new status:
                           </p>
                          <Select onValueChange={setSelectedStatus}>
                              <SelectTrigger className="w-fit">
                                <SelectValue placeholder="gb_status" />
                              </SelectTrigger>
                              <SelectContent>
                                <SelectGroup>
                                  <SelectItem value="abandoned">Abandoned</SelectItem>
                                  <SelectItem value="poor_genome_busco">Low genome BUSCO</SelectItem>
                                  <SelectItem value="insufficient_data">Insufficient data</SelectItem>
                                    <SelectItem value="check_busco">Low protein BUSCO</SelectItem>
                                </SelectGroup>
                              </SelectContent>
                            </Select>
                          {/* Add selected count info */}
                        <p className="mt-4 text-sm text-foreground">
                          {Object.keys(rowSelection).length} record
                          {Object.keys(rowSelection).length !== 1 ? 's' : ''} selected. These will be changed in the registry.
                        </p>
                      </div>
                      <DialogFooter>
                          <DialogClose asChild>
                            <Button variant="outline">Cancel</Button>
                          </DialogClose>
                          <Button
                            type="button"
                            onClick={async () => {
                              const selectedRows = table.getSelectedRowModel().rows;
                              if (!genebuilder) {
                                  alert("No genebuilder specified");
                                  return;
                                }
                              if (!selectedRows.length) return;
                              if (!selectedStatus) {
                              alert("Please select a new status");
                              return;
                            }


                              // Extract the gca info
                              const selectedItems = selectedRows.map(row => ({
                                  gca: row.original.gca,
                                  annotation_method: row.original.annotation_method,
                                }));
                              console.log("Sending GCAs:", selectedItems);

                              try {
                                const response = await fetch("/api/handover/handover/change_status", {
                                  method: "POST",
                                  headers: {
                                    "Content-Type": "application/json",
                                  },
                                  body: JSON.stringify({ genebuilder, items: selectedItems, new_status: selectedStatus }),
                                });

                                if (!response.ok) {
                                  throw new Error("Failed to update status");
                                }

                                // Handle success
                                alert(`Updated ${selectedItems.length} records successfully`);
                                setIsDialogOpen(false);
                                table.resetRowSelection();
                                // Reload table
                                await reloadData()
                              } catch (error) {
                                console.error(error);
                                alert("Error updating records");
                              }
                            }}
                          >
                            Apply
                          </Button>
                        </DialogFooter>
                  </DialogContent>
                  </Dialog>
          </CardAction>

      </CardHeader>
      <CardContent className="mb-4 min-w-0 overflow-x-auto">
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
