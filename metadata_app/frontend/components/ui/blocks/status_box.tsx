"use client"

import * as React from "react"
import Link from 'next/link';
import {ArrowRight, ChartNoAxesCombined, Ban, Users} from 'lucide-react';
import { BackgroundGradient } from "@/components/ui/backround-gradient"
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table"
import {
  Drawer,
  DrawerContent,
  DrawerDescription,
  DrawerHeader,
  DrawerTitle,
  DrawerTrigger,
} from "@/components/ui/drawer"
import {
    Card,
    CardContent,
} from "@/components/ui/card"
import {ColumnDef, getCoreRowModel, getFilteredRowModel, getSortedRowModel, useReactTable, flexRender} from "@tanstack/react-table";



export type NumHOItem = {
  value: number
}
export type NumDataItem = {
  value: number
}
export type NumPENItem = {
  value: number
}

export type DataItem = {
  gca: string;
  scientific_name: string;
  gb_status: string;
}

export type PendingItem = {
  gca: string;
  scientific_name: string;
  gb_status: string;
}

type StatusBoxProps = {
  horeadyCount: number;
  dataCount: number;
  pendingCount: number;
  listData: DataItem[];
  listPending: PendingItem[];
}

export const columns: ColumnDef<DataItem>[] = [
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
    accessorKey: "gb_status",
    header: "Status",
    cell: ({ row }) => (
      <div >{row.getValue("gb_status")}</div>
    ),
  },
]

export const columns_pending: ColumnDef<PendingItem>[] = [
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
    accessorKey: "gb_status",
    header: "Status",
    cell: ({ row }) => (
      <div >{row.getValue("gb_status")}</div>
    ),
  },
]

export function StatusBox({ horeadyCount, dataCount, pendingCount, listData, listPending}: StatusBoxProps) {
  const table = useReactTable({
    data: listData,
    columns,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  });

    const table_pending = useReactTable({
    data: listPending,
    columns: columns_pending,
    getCoreRowModel: getCoreRowModel(),
    getSortedRowModel: getSortedRowModel(),
    getFilteredRowModel: getFilteredRowModel(),
  });

    return (
      <div className="flex items-center justify-center p-8">
          <div className="@container grow w-full">

      <div className="grid grid-cols-1 @3xl:grid-cols-3 gap-8 w-full">
          <BackgroundGradient>
          <Card className="rounded-xl overflow-hidden shadow-lg p-0 border-0">
          <CardContent className="relative overflow-hidden flex flex-col justify-end py-6 px-0 pb-0">
                {/* Icon */}
                <div className="px-6 mb-3.5">
                  <ChartNoAxesCombined className="size-8 text-foreground" />
                </div>
                {/* Main content */}
              <div className="flex-1 flex flex-col justify-center items-start px-6">
                  <div className="text-foreground text-4xl font-bold mb-6">{horeadyCount.toLocaleString()}</div>
                  <div className="text-foreground text-lg font-semibold mb-1">Handover ready</div>
                  <div className="text-foreground text-sm mb-2">Cores ready to be handed over</div>
                </div>
              <Link
                  href="#ho-ready-table"
                  className="group/card w-full bg-muted-foreground px-6 py-4 flex items-center justify-between mt-6"
                >
                  <span className="text-background text-sm font-medium"> Jump to table</span>
                  <ArrowRight className="group-hover/card:translate-x-1 transition-transform duration-300 w-5 h-5 text-background" />
                </Link>
          </CardContent>
        </Card>
              </BackgroundGradient>

          <BackgroundGradient>
              <Card className="rounded-xl overflow-hidden shadow-lg p-0 border-0">
                <CardContent className="relative overflow-hidden flex flex-col justify-end py-6 px-0 pb-0">
                  {/* Icon */}
                  <div className="px-6 mb-3.5">
                    <Users className="size-8 text-foreground" />
                  </div>

                  {/* Main content */}
                  <div className="flex-1 flex flex-col justify-center items-start px-6">
                    <div className="text-foreground text-4xl font-bold mb-6">
                      {dataCount.toLocaleString()}
                    </div>
                    <div className="text-foreground text-lg font-semibold mb-1">
                      Check for more data
                    </div>
                    <div className="text-foreground text-sm mb-2">
                      Low BUSCO or insufficient data
                    </div>
                  </div>

                  {/* Drawer trigger */}
                  <Drawer>
                    <DrawerTrigger asChild>
                      <button className="group/card w-full bg-muted-foreground px-6 py-4 flex items-center justify-between mt-6">
                        <span className="text-background text-sm font-medium">Show list</span>
                        <ArrowRight className="group-hover/card:translate-x-1 transition-transform duration-300 w-5 h-5 text-background" />
                      </button>
                    </DrawerTrigger>

                    <DrawerContent>
                      <DrawerHeader>
                        <DrawerTitle>Data List</DrawerTitle>
                        <DrawerDescription>
                          Here is the list of annotations that have been stuck for more then 6 months due to low data.
                        </DrawerDescription>
                      </DrawerHeader>

                        <Card className="mx-16 py-2">
                        <CardContent className="max-h-[60vh] overflow-y-auto mb-4">

                      <div className={
                            listData && listData.length
                              ? "mt-4 rounded-md overflow-y-auto"
                              : "mt-4 rounded-md flex items-center justify-center h-[60vh]"
                          }>
                        {listData && listData.length ? (
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
                        ) : (
                          <p>No data available.</p>
                        )}
                      </div>
                        </CardContent>
                            </Card>
                    </DrawerContent>
                  </Drawer>
                </CardContent>
              </Card>
            </BackgroundGradient>

          <BackgroundGradient>
          <Card className="rounded-xl overflow-hidden shadow-lg p-0 border-0">
          <CardContent className="relative overflow-hidden flex flex-col justify-end py-6 px-0 pb-0">
                {/* Icon */}
                <div className="px-6 mb-3.5">
                  <Ban className="size-8 text-foreground" />
                </div>
                {/* Main content */}
              <div className="flex-1 flex flex-col justify-center items-start px-6">
                  <div className="text-foreground text-4xl font-bold mb-6">{pendingCount.toLocaleString()}</div>
                  <div className="text-foreground text-lg font-semibold mb-1">Pending</div>
                  <div className="text-foreground text-sm mb-2">In progress for more than 6 months</div>
                </div>
              {/* Drawer trigger */}
                  <Drawer>
                    <DrawerTrigger asChild>
                      <button className="group/card w-full bg-muted-foreground px-6 py-4 flex items-center justify-between mt-6">
                        <span className="text-background text-sm font-medium">Show list</span>
                        <ArrowRight className="group-hover/card:translate-x-1 transition-transform duration-300 w-5 h-5 text-background" />
                      </button>
                    </DrawerTrigger>

                    <DrawerContent>
                      <DrawerHeader>
                        <DrawerTitle>Data List</DrawerTitle>
                        <DrawerDescription>
                          Here is the list of annotations that have been stuck for more then 6 months as in progress.
                        </DrawerDescription>
                      </DrawerHeader>
                    <Card className="mx-16 py-2">
                        <CardContent className="max-h-[60vh] overflow-y-auto mb-4">
                      <div className={
                            listData && listData.length
                              ? "mt-4 rounded-md overflow-y-auto"
                              : "mt-4 rounded-md flex items-center justify-center h-[60vh]"
                          }>
                        {listPending && listPending.length ? (
                          <Table>
                              <TableHeader>
                                {table_pending.getHeaderGroups().map((headerGroup) => (
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
                                {table_pending.getRowModel().rows.length ? (
                                  table_pending.getRowModel().rows.map((row) => (
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
                                      colSpan={columns_pending.length}
                                      className="h-24 text-center"
                                    >
                                      No results.
                                    </TableCell>
                                  </TableRow>
                                )}
                              </TableBody>
                            </Table>
                        ) : (
                          <p>No data available.</p>
                        )}
                      </div>
                        </CardContent>
                    </Card>
                    </DrawerContent>
                  </Drawer>
          </CardContent>
        </Card>
              </BackgroundGradient>

          </div>
          </div>
           </div>
  )
}