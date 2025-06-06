"use client"

import {Bar, BarChart, CartesianGrid, LabelList, ResponsiveContainer, XAxis} from "recharts"
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

export type CladeItem = {
  internal_clade: string
  count: number
}
type Props = {
  data: CladeItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  count: { label: "Annotations" },
} satisfies ChartConfig



export function RepClade({ data }: Props) {
const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
    }))
  }, [data])

  const totalAnnotations = React.useMemo(() => {
    return transformedData.reduce((sum, item) => sum + (item.count || 0), 0)
  }, [transformedData])


  return (
    <Card>
      <CardHeader>
        <CardTitle>Associated internal clades</CardTitle>
        <CardDescription>Number of assemblies per clade</CardDescription>
      </CardHeader>
      <CardContent>
          <div style={{ height: 300 }}>
        <ChartContainer config={chartConfig}>
            <div style={{ height: 300 }}>
                <ResponsiveContainer width="100%" height={300}>
        <BarChart accessibilityLayer data={transformedData} margin={{
              top: 30,
            }} >
          <CartesianGrid vertical={false} />
          <XAxis
              dataKey="internal_clade"
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
                </ResponsiveContainer>
            </div>
      </ChartContainer>
              </div>
      </CardContent>
    </Card>
  )
}