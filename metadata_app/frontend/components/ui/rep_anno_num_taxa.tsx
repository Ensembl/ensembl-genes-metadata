"use client"

import * as React from "react"
import {
  Card,
  CardHeader,
  CardTitle,
  CardContent,
} from "@/components/ui/card"

export type NumTaxaItem = {
  value: number
}

type Props = {
  data: NumTaxaItem | null
}

export function AnnotatedTaxaCard({ data }: Props) {
  if (!data) return null

  return (
    <Card>
  <CardHeader>
    <CardTitle>Number of Annotated Taxa</CardTitle>
  </CardHeader>
  <CardContent className="flex items-center justify-center h-60">
    <div className="text-7xl font-bold text-primary text-center">
      {data.value.toLocaleString()}
    </div>
  </CardContent>
</Card>
  )
}