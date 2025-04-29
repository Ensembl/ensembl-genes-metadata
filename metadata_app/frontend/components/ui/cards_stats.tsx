"use client"

import React from "react"
import { XAxis, Line, LineChart, CartesianGrid } from "recharts"

import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
  CardDescription,
} from "@/components/ui/card"
import { ChartConfig, ChartContainer, ChartTooltip, ChartTooltipContent } from "@/components/ui/chart"

const chartConfig = {
  assembly_count: {
    label: "Assemblies",
    color: "hsl(var(--primary))",
  },
  annotation_count: {
    label: "Annotations",
    color: "hsl(var(--primary))",
  },
} satisfies ChartConfig

export function CardsStats() {
  const [assemblyData, setAssemblyData] = React.useState([])
  const [annotationData, setAnnotationData] = React.useState([])

  React.useEffect(() => {
    const fetchAssemblies = async () => {
      const res = await fetch("/api/assemblies")
      const json = await res.json()
      setAssemblyData(json)
    }

    const fetchAnnotations = async () => {
      const res = await fetch("/api/annotations")
      const json = await res.json()
      setAnnotationData(json)
    }

    fetchAssemblies()
    fetchAnnotations()
  }, [])

  return (
    <div className="grid grid-cols-2 gap-4">
      <Card>
        <CardHeader>
          <CardTitle className="text-lg font-bold">Assemblies</CardTitle>
          <CardDescription className="text-sm font-light mb-2">
            Number of assemblies by year
          </CardDescription>
        </CardHeader>
        <CardContent className="pb-0">
          <ChartContainer config={chartConfig} className="h-[80px] w-full">
            <LineChart
              data={assemblyData}
              margin={{ top: 5, right: 10, left: 10, bottom: 4 }}
            >
              <ChartTooltip
              cursor={false}
              content={
                <ChartTooltipContent
                  labelFormatter={(label, payload) =>
                    payload?.[0]?.payload?.year ?? label
                  }
                />
              }
            />
              <Line
                type="natural"
                strokeWidth={2}
                dataKey="assembly_count"
                stroke="var(--chart-1)"
                activeDot={{ r: 6 }}
              />
            </LineChart>
          </ChartContainer>
        </CardContent>
      </Card>

      <Card>
        <CardHeader>
          <CardTitle className="text-lg font-bold">Annotations</CardTitle>
          <CardDescription className="text-sm font-light mb-2">
            Number of annotations by year
          </CardDescription>
        </CardHeader>
        <CardContent className="pb-0">
          <ChartContainer config={chartConfig} className="h-[80px] w-full">
            <LineChart
              data={annotationData}
              margin={{ top: 5, right: 10, left: 10, bottom: 4}}
            >
              <ChartTooltip
              cursor={false}
              content={
                <ChartTooltipContent
                  labelFormatter={(label, payload) =>
                    payload?.[0]?.payload?.year ?? label
                  }
                />
              }
            />
              <Line
                type="monotone"
                strokeWidth={2}
                dataKey="annotation_count"
                stroke="var(--chart-1)"
                activeDot={{ r: 6 }}
              />
            </LineChart>
          </ChartContainer>
        </CardContent>
      </Card>
    </div>
  )
}