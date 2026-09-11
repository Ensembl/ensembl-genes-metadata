"use client"

import {Bar, BarChart, CartesianGrid, LabelList, XAxis} from "recharts"
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"

import {
  ChartConfig,
  ChartContainer,
  ChartTooltip,
  ChartTooltipContent,
} from "@/components/ui/chart"
import * as React from "react";

export type AsmLevelItem = {
  asm_level: string
  count: number
}
type Props = {
  data: AsmLevelItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  count: { label: "Assemblies" },
  "Contig": { label: "Contig", color: "var(--chart-1)" },
  "Chromosome": { label: "Chromosome", color: "var(--chart-1)" },
  "Scaffold": { label: "Scaffold", color: "var(--chart-1)" },
  "Complete genome": { label: "Complete genome", color: "var(--chart-1)" },
} satisfies ChartConfig


export function RepAsmLevel({ data }: Props) {
const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.asm_level]?.label || item.asm_level,
      fill: chartConfig[item.asm_level]?.color,
    }))
  }, [data])

  const totalAnnotations = React.useMemo(() => {
    return transformedData.reduce((sum, item) => sum + (item.count || 0), 0)
  }, [transformedData])


  return (
    <Card>
      <CardHeader>
        <CardTitle>Assembly level</CardTitle>
        <CardDescription>Number of assemblies per level</CardDescription>
      </CardHeader>
      <CardContent>
        <ChartContainer config={chartConfig}>
        <BarChart accessibilityLayer data={transformedData} margin={{
              top: 30,
            }} >
          <CartesianGrid vertical={false} />
          <XAxis
              dataKey="displayName"
              tickLine={false}
              axisLine={false}
          />
          <ChartTooltip
              cursor={false}
              content={<ChartTooltipContent />}
            />
          <Bar dataKey="count" fill="var(--color-chart-1)" radius={8}>
              <LabelList
                position="top"
                offset={12}
                className="fill-foreground"
                fontSize={12}
              />
            </Bar>
        </BarChart>
      </ChartContainer>
      </CardContent>
    </Card>
  )
}