"use client"
import * as React from "react"
import {
  Card,
  CardContent, CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import {
  ChartConfig,
  ChartContainer,
  ChartTooltip,
  ChartTooltipContent,
} from "@/components/ui/chart"

import {
  ScatterChart,
  Scatter,
  CartesianGrid,
  XAxis,
  YAxis,
  ResponsiveContainer,
} from "recharts"

const chartConfig: Record<string, { label: string; color?: string }> = {
  total_sequence_length_Gb: { label: "Genome size", color: "var(--chart-1)" },
} satisfies ChartConfig

export type LengthItem = {
  total_sequence_length_Gb: number
  gca: string
}

type Props = {
  data: LengthItem[]
}

export function LengthChart({ data }: Props) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>Genome size</CardTitle>
        <CardDescription>Assembly ungapped length in Gb</CardDescription>
      </CardHeader>

      <CardContent>
        <div style={{ height: 300 }}>
        <ChartContainer config={chartConfig}>
          <div style={{ height: 300 }}>
          <ResponsiveContainer width="100%" height={300}>
            <ScatterChart
              data={data}
              margin={{ top: 30, right: 30 }}
            >
              <CartesianGrid vertical={false} />
              <XAxis
                type="category"
                dataKey="gca"
                name="Assembly"
                tick={false}
                tickLine={false}
                axisLine={false}
              />
              <YAxis
                type="number"
                dataKey="total_sequence_length_Gb"
                unit="Gb"
                name="Genome size (Gb)"
                tickLine={false}
                axisLine={false}
              />
              <ChartTooltip content={<ChartTooltipContent />} />
              <Scatter
                name="Genome size"
                data={data}
                fill="var(--chart-1)"
              />
            </ScatterChart>
          </ResponsiveContainer>
          </div>
        </ChartContainer>
           </div>
      </CardContent>
    </Card>
  )
}