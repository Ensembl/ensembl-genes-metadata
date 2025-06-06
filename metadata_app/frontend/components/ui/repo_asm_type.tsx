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

export type AsmTypeItem = {
  asm_type: string
  count: number
}
type Props = {
  data: AsmTypeItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  count: { label: "Assemblies" },
  haploid: { label: "haploid", color: "var(--chart-1)" },
  "alternate-pseudohaplotype": { label: "alternate-pseudohaplotype", color: "var(--chart-2)" },
  "unresolved-diploid": { label: "unresolved-diploid", color: "var(--chart-3)" },
  "haploid-with-alt-loci": { label: "haploid-with-alt-loci", color: "var(--chart-4)" },
  diploid: { label: "diploid", color: "var(--chart-5)" },
} satisfies ChartConfig



export function RepAsmType({ data }: Props) {
const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.asm_type]?.label || item.asm_type,
      fill: chartConfig[item.asm_type]?.color,
    }))
  }, [data])

  const totalAnnotations = React.useMemo(() => {
    return transformedData.reduce((sum, item) => sum + (item.count || 0), 0)
  }, [transformedData])


  return (
    <Card>
      <CardHeader>
        <CardTitle>Assembly types</CardTitle>
        <CardDescription>Number of assemblies per type</CardDescription>
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