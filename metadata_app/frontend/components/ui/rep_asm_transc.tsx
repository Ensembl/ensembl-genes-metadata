"use client"

import * as React from "react"
import {
  Card,
  CardHeader,
  CardTitle,
  CardContent, CardDescription,
} from "@/components/ui/card"

export type TranscItem = {
  value: number
}

type Props = {
  data: TranscItem | null
}

export function TranscCard({ data }: Props) {
  if (!data) return null

  return (
    <Card>
  <CardHeader>
    <CardTitle>Transcriptomic registry</CardTitle>
    <CardDescription>Percentage of assemblies that have been assessed in the transcriptomic registry species or genus level</CardDescription>
  </CardHeader>
  <CardContent className="flex items-center justify-center h-60">
    <div className="text-7xl font-bold text-primary text-center">
      {data.value.toLocaleString()}
    </div>
  </CardContent>
</Card>
  )
}