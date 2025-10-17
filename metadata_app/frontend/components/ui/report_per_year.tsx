"use client"
import { Bar, BarChart, CartesianGrid, XAxis } from "recharts"

import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import {
  ChartTooltip,
  ChartTooltipContent,
  ChartContainer,
} from "@/components/ui/chart"
import React from "react"

const chartConfig = {
  AEGIS: {
    label: "AEGIS",
  },
  "AQUA-FAANG": {
    label: "AQUA-FAANG",
  },
  ASG: {
    label: "ASG",
  },
  CBP: {
    label: "CBP",
  },
  DToL: {
    label: "DToL",
  },
  EBP: {
    label: "EBP",
  },
  ERGA: {
    label: "ERGA",
  },
  "ERGA/BGE": {
    label: "ERGA/BGE",
  },
  ERGA_pilot: {
    label: "ERGA Pilot",
  },
  HPRC: {
    label: "HPRC",
  },
  LACA: {
    label: "LACA",
  },
  VGP: {
    label: "VGP",
  },
}

export function ProjectsBar() {
  const [assemblyData, setAssemblyData] = React.useState<any[]>([])
  const [annotationData, setAnnotationData] = React.useState<any[]>([])

  React.useEffect(() => {
    const fetchData = async () => {
      try {
        const res = await fetch("/api/report/asm/report/asm/bar")
        const json = await res.json()

        setAssemblyData(json.df_assembly)
        setAnnotationData(json.df_annotations)
      } catch (err) {
        console.error("Error fetching chart data:", err)
      }
    }
    fetchData()
  }, [])

  return (
      <div className="grid grid-cols-2 gap-4">
        <Card>
          <CardHeader>
            <CardTitle>Assemblies per year per project</CardTitle>
            <CardDescription>Biodiversity projects</CardDescription>
          </CardHeader>
          <CardContent>
              <ChartContainer config={chartConfig} className="h-[400px] w-full">
                <BarChart accessibilityLayer data={assemblyData}>
                  <CartesianGrid vertical={false} />
                  <XAxis
                    dataKey="release_year"
                    tickLine={false}
                    tickMargin={10}
                    axisLine={false}
                    tickFormatter={(value) => value.toString()}
                  />
                  <ChartTooltip cursor={false} content={<ChartTooltipContent indicator="dashed" />} />
                  <Bar dataKey="AEGIS" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ASG" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="CBP" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="DToL" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="EBP" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ERGA" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ERGA/BGE" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ERGA_pilot" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="VGP" fill="var(--chart-7)" radius={4} />
                </BarChart>
              </ChartContainer>
          </CardContent>
          </Card>

        <Card>
          <CardHeader>
            <CardTitle>Assemblies per year per project</CardTitle>
            <CardDescription>Custom groups</CardDescription>
          </CardHeader>
          <CardContent>
              <ChartContainer config={chartConfig} className="h-[400px] w-full">
                <BarChart accessibilityLayer data={assemblyData}>
                  <CartesianGrid vertical={false} />
                  <XAxis
                    dataKey="release_year"
                    tickLine={false}
                    tickMargin={10}
                    axisLine={false}
                    tickFormatter={(value) => value.toString()}
                  />
                  <ChartTooltip cursor={false} content={<ChartTooltipContent indicator="dashed" />} />
                  <Bar dataKey="AQUA-FAANG" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="HPRC" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="LACA" fill="var(--chart-7)" radius={4} />
                </BarChart>
              </ChartContainer>
          </CardContent>
          </Card>

          <Card>
            <CardHeader>
            <CardTitle>Annotations in Beta per year per project</CardTitle>
            <CardDescription>Biodiversity projects</CardDescription>
          </CardHeader>
          <CardContent>
              <ChartContainer config={chartConfig} className="h-[400px] w-full">
                <BarChart accessibilityLayer data={annotationData}>
                  <CartesianGrid vertical={false} />
                  <XAxis
                    dataKey="release_year"
                    tickLine={false}
                    tickMargin={10}
                    axisLine={false}
                    tickFormatter={(value) => value.toString()}
                  />
                  <ChartTooltip cursor={false} content={<ChartTooltipContent indicator="dashed" />} />
                  <Bar dataKey="AEGIS" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ASG" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="CBP" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="DToL" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="EBP" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ERGA" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ERGA/BGE" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="ERGA_pilot" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="VGP" fill="var(--chart-7)" radius={4} />
                </BarChart>
              </ChartContainer>
          </CardContent>
            </Card>

        <Card>
          <CardHeader>
            <CardTitle>Annotations in Beta per year per project</CardTitle>
            <CardDescription>Custom groups</CardDescription>
          </CardHeader>
          <CardContent>
              <ChartContainer config={chartConfig} className="h-[400px] w-full">
                <BarChart accessibilityLayer data={annotationData}>
                  <CartesianGrid vertical={false} />
                  <XAxis
                    dataKey="release_year"
                    tickLine={false}
                    tickMargin={10}
                    axisLine={false}
                    tickFormatter={(value) => value.toString()}
                  />
                  <ChartTooltip cursor={false} content={<ChartTooltipContent indicator="dashed" />} />
                  <Bar dataKey="AQUA-FAANG" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="HPRC" fill="var(--chart-7)" radius={4} />
                  <Bar dataKey="LACA" fill="var(--chart-7)" radius={4} />
                </BarChart>
              </ChartContainer>
          </CardContent>
        </Card>
      </div>
  )
}