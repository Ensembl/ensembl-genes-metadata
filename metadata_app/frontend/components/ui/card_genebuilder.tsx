"use client"

import * as React from "react"
import { Label, Pie, PieChart } from "recharts"

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
const chartData = [
  { genebuilder: "Anna", annotations: 275, fill: "var(--color-anna)" },
  { genebuilder: "Jose", annotations: 200, fill: "var(--color-jose)" },
  { genebuilder: "Francesca", annotations: 287, fill: "var(--color-francesca)" },
  { genebuilder: "Vianey", annotations: 173, fill: "var(--color-vianey)" },
  { genebuilder: "Swati", annotations: 190, fill: "var(--color-swati)" },
  { genebuilder: "Jack", annotations: 190, fill: "var(--color-jack)" },

]

const chartConfig = {
  annotations: {
    label: "Annotations",
  },
  anna: {
    label: "Anna",
    color: "var(--chart-1)",
  },
  jose: {
    label: "Jose",
    color: "var(--chart-2)",
  },
  francesca: {
    label: "Francesca",
    color: "var(--chart-3)",
  },
  vianey: {
    label: "Vianey",
    color: "var(--chart-4)",
  },
  swati: {
    label: "Swati",
    color: "var(--chart-5)",
  },
  jack: {
    label: "Jack",
    color: "var(--chart-6)",
  },

} satisfies ChartConfig

export function CardGeneBuilder() {
  const totalAnnotations = React.useMemo(() => {
    return chartData.reduce((acc, curr) => acc + curr.annotations, 0)
  }, [])

  return (
    <Card className="flex flex-col">
      <CardHeader className="items-center pb-0">
        <CardTitle>Top Genebuilders</CardTitle>
        <CardDescription>Last 30 days</CardDescription>
      </CardHeader>
      <CardContent className="flex-1 pb-0">
        <ChartContainer
          config={chartConfig}
          className="mx-auto aspect-square max-h-[250px]"
        >
          <PieChart>
            <ChartTooltip
              cursor={false}
              content={<ChartTooltipContent hideLabel />}
            />
            <Pie
              data={chartData}
              dataKey="annotations"
              nameKey="genebuilder"
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
                }}
              />
            </Pie>
          </PieChart>
        </ChartContainer>
      </CardContent>
    </Card>
  )
}
