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

export type TranscENAItem = {
  transcriptomic_evidence: string
  count: number
}
type Props = {
  data: TranscENAItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  count: { label: "Assemblies" },
  yes: { label: "Yes", color: "var(--chart-1)" },
  no: { label: "No", color: "var(--chart-2)" },
  "not checked": { label: "Not checked", color: "var(--chart-3)" },
} satisfies ChartConfig



export function RepTranscENA({ data }: Props) {
const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.transcriptomic_evidence]?.label || item.transcriptomic_evidence,
      fill: chartConfig[item.transcriptomic_evidence]?.color,
    }))
  }, [data])

  return (
    <Card>
      <CardHeader>
        <CardTitle>Transcriptomic evidence check</CardTitle>
        <CardDescription>Number of assemblies that have species or genus level transcriptomic evidence in ENA</CardDescription>
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