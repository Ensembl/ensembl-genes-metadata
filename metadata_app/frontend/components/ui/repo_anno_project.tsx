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
import {MethodItem} from "@/components/ui/rep_anno_method";
import * as React from "react";

export type ProjectItem = {
  associated_project: string
  number_of_annotations: number
}
type Props = {
  data: ProjectItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  number_of_annotations: { label: "Annotations" },
  DToL: { label: "DToL", color: "var(--chart-1)" },
  "ERGA/BGE": { label: "ERGA/BG", color: "var(--chart-2)" },
  ERGA: { label: "ERGA", color: "var(--chart-3)" },
  EBP: { label: "EBP", color: "var(--chart-4)" },
  ERGA_pilot: { label: "ERGA_pilot", color: "var(--chart-5)" },
  ASG: { label: "ASG", color: "var(--chart-6)" },
  VGP: { label: "VGP", color: "var(--chart-7)" },
  CBP: { label: "CBP", color: "var(--chart-7)" },
} satisfies ChartConfig



export function RepProject({ data }: Props) {
const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.associated_project]?.label || item.associated_project,
      fill: chartConfig[item.associated_project]?.color,
    }))
  }, [data])

  const totalAnnotations = React.useMemo(() => {
    return transformedData.reduce((sum, item) => sum + (item.number_of_annotations || 0), 0)
  }, [transformedData])


  return (
    <Card>
      <CardHeader>
        <CardTitle>Associated biodiversity projects</CardTitle>
        <CardDescription>Number of annotations per project</CardDescription>
      </CardHeader>
      <CardContent>
        <ChartContainer config={chartConfig}>
        <BarChart accessibilityLayer data={transformedData} margin={{
              top: 20,
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
          <Bar dataKey="number_of_annotations" fill="var(--color-chart-1)" radius={8}>
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