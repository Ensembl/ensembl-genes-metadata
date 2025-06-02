"use client"

import * as React from "react"
import { Label, Pie, PieChart } from "recharts"
import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import {
  ChartConfig,
  ChartContainer,
  ChartTooltip,
  ChartTooltipContent,
} from "@/components/ui/chart"

const chartConfig: Record<string, { label: string; color?: string }> = {
  number_of_annotations: { label: "Annotations" },
  anno: { label: "Anno", color: "var(--chart-1)" },
  braker: { label: "Braker", color: "var(--chart-2)" },
  external_annotation_import: { label: "External Annotation Import", color: "var(--chart-3)" },
  full_genebuild: { label: "Full Genebuild", color: "var(--chart-4)" },
  import: { label: "Import", color: "var(--chart-5)" },
  mixed_strategy_build: { label: "Mixed Strategy Build", color: "var(--chart-6)" },
  projection_build: { label: "Projection Build", color: "var(--chart-7)" },
} satisfies ChartConfig

export type MethodItem = {
  annotation_method: string
  number_of_annotations: number
}

type Props = {
  data: MethodItem[]
}

export function AnoMethodSummaryChart({ data }: Props) {
  const transformedData = React.useMemo(() => {
    return data.map((item) => ({
      ...item,
      displayName: chartConfig[item.annotation_method]?.label || item.annotation_method,
      fill: chartConfig[item.annotation_method]?.color,
    }))
  }, [data])

  const totalAnnotations = React.useMemo(() => {
    return transformedData.reduce((sum, item) => sum + (item.number_of_annotations || 0), 0)
  }, [transformedData])

  return (
    <Card className="flex flex-col">
      <CardHeader className="pb-0">
        <CardTitle>Annotation method summary</CardTitle>
      </CardHeader>

      <CardContent className="flex-1 pb-0">
        {transformedData.length > 0 ? (
          <ChartContainer
            config={chartConfig}
            className="mx-auto aspect-square max-h-[250px]"
          >
            <PieChart>
              <ChartTooltip cursor={false} content={<ChartTooltipContent hideLabel />} />
              <Pie
                data={transformedData}
                dataKey="number_of_annotations"
                nameKey="displayName"
                innerRadius={60}
                strokeWidth={5}
              >
                <Label
                  content={({ viewBox }) => {
                    if (viewBox && "cx" in viewBox && "cy" in viewBox) {
                      return (
                        <text
                          x={viewBox.cx}
                          y={viewBox.cy}
                          textAnchor="middle"
                          dominantBaseline="middle"
                        >
                          <tspan
                            x={viewBox.cx}
                            y={viewBox.cy}
                            className="fill-foreground text-3xl font-bold"
                          >
                            {totalAnnotations.toLocaleString()}
                          </tspan>
                          <tspan
                            x={viewBox.cx}
                            y={(viewBox.cy || 0) + 24}
                            className="fill-muted-foreground"
                          >
                            Annotations
                          </tspan>
                        </text>
                      )
                    }
                    return null
                  }}
                />
              </Pie>
            </PieChart>
          </ChartContainer>
        ) : (
          <p className="text-muted-foreground text-center">No data available</p>
        )}
      </CardContent>
    </Card>
  )
}