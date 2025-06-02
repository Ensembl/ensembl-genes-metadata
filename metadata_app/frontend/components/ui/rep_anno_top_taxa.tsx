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

export type TaxaItem = {
  scientific_name: string
  count: number
};

type Props = {
  data: TaxaItem[]
}

const chartConfig: Record<string, { label: string }> = {
  count: { label: "Annotations" },
} satisfies ChartConfig



export function RepTopTaxa({ data }: Props) {
  const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.scientific_name]?.label || item.scientific_name,
    }))
  }, [data])


  return (
    <Card>
      <CardHeader>
        <CardTitle>Top 3 taxa</CardTitle>
        <CardDescription>The most frequently annotated taxa</CardDescription>
      </CardHeader>
      <CardContent>
        <ChartContainer config={chartConfig}>
        <BarChart accessibilityLayer data={transformedData} margin={{
              top: 30,
            }}>
          <CartesianGrid vertical={false} />
          <XAxis
              dataKey="scientific_name"
              tickLine={false}
              axisLine={false}
          />
          <ChartTooltip
              cursor={false}
              content={<ChartTooltipContent hideLabel />}
            />
          <Bar dataKey="count" fill="var(--chart-1)" radius={8} >
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