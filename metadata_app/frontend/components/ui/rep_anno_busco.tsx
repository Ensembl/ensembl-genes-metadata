"use client"

import * as React from "react"
import {
  Card,
  CardHeader,
  CardTitle,
  CardContent,
} from "@/components/ui/card"

export type BuscoItem = {
  value: number
}

type Props = {
  data: BuscoItem | null
}

export function AnnotatedBuscoCard({ data }: Props) {
  if (!data) return null

  return (
    <Card>
  <CardHeader>
    <CardTitle>Average BUSCO score</CardTitle>
  </CardHeader>
  <CardContent className="flex items-center justify-center h-60">
    <div className="text-7xl font-bold text-primary text-center">
      {data.value.toLocaleString()}
    </div>
  </CardContent>
</Card>
  )
}