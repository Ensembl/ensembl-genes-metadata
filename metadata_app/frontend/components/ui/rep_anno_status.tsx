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

export type StatusItem = {
  gb_status: string
  count: number
};

type Props = {
  data: StatusItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  count: { label: "Annotations" },
  pre_released: { label: "Pre-released", color: "var(--chart-3)" },
    archive: { label: "Archived", color: "var(--chart-3)" },
  live: { label: "Live", color: "var(--chart-3)" },
  handed_over: { label: "Handed over", color: "var(--chart-3)" },
  in_progress: { label: "In progress", color: "var(--chart-3)" },
  completed: { label: "Completed", color: "var(--chart-3)" },
    check_busco: { label: "Low pBUSCO", color: "var(--chart-3)"},
  insufficient_data: { label: "Low evidence", color: "var(--chart-3)" },
    poor_genome_busco: { label: "Low gBUSCO", color: "var(--chart-3)" },
    abandoned: { label: "Abandoned", color: "var(--chart-3)" }
} satisfies ChartConfig



export function RepStatus({ data }: Props) {
  const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.gb_status]?.label || item.gb_status,
      fill: chartConfig[item.gb_status]?.color,
    }))
  }, [data])

  const totalAnnotations = React.useMemo(() => {
    return transformedData.reduce((sum, item) => sum + (item.count || 0), 0)
  }, [transformedData])


  return (
    <Card>
      <CardHeader>
        <CardTitle>Annotation status</CardTitle>
        <CardDescription>Current status of annotations</CardDescription>
      </CardHeader>
      <CardContent>
        <ChartContainer config={chartConfig}>
        <BarChart accessibilityLayer data={transformedData} margin={{
              top: 20,
            bottom: 80,
            }}>
          <CartesianGrid vertical={false} />
          <XAxis
              dataKey="displayName"
              tickLine={false}
              axisLine={false}
              angle={-90}
              textAnchor="end"
              alignmentBaseline="middle"
          />
          <ChartTooltip
              cursor={false}
              content={<ChartTooltipContent />}
            />
          <Bar dataKey="count" fill="var(--color-chart-1)" radius={8} >
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