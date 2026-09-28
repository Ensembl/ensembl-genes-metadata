"use client"
import { Bar, BarChart, CartesianGrid, XAxis } from "recharts"
import { PROJECTS } from "@/features/projects/project-config"

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
  type ChartConfig,
} from "@/components/ui/chart"
import React from "react"

const reportProjects = PROJECTS.filter((project) => project.reportKey)
const biodiversityProjects = reportProjects.filter(
  (project) => project.reportGroup === "biodiversity",
)
const customProjects = reportProjects.filter(
  (project) => project.reportGroup === "custom",
)
const chartConfig = Object.fromEntries(
  reportProjects.map((project) => [
    project.reportKey,
    {
      label: project.reportLabel ?? project.reportKey,
    },
  ]),
) satisfies ChartConfig

const renderProjectBars = (projects: typeof reportProjects) =>
  projects.map((project) => (
    <Bar key={project.reportKey} dataKey={project.reportKey} fill="var(--chart-3)" radius={4} />
  ))

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
      <div className="grid min-w-0 grid-cols-1 gap-4 lg:grid-cols-2">
        <Card className="min-w-0">
          <CardHeader>
            <CardTitle>Assemblies per year per project</CardTitle>
            <CardDescription>Biodiversity projects</CardDescription>
          </CardHeader>
          <CardContent className="min-w-0">
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
                  {renderProjectBars(biodiversityProjects)}
                </BarChart>
              </ChartContainer>
          </CardContent>
          </Card>

        <Card className="min-w-0">
          <CardHeader>
            <CardTitle>Assemblies per year per project</CardTitle>
            <CardDescription>Custom groups</CardDescription>
          </CardHeader>
          <CardContent className="min-w-0">
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
                  {renderProjectBars(customProjects)}
                </BarChart>
              </ChartContainer>
          </CardContent>
          </Card>

          <Card className="min-w-0">
            <CardHeader>
            <CardTitle>Annotations in Beta per year per project</CardTitle>
            <CardDescription>Biodiversity projects</CardDescription>
          </CardHeader>
          <CardContent className="min-w-0">
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
                  {renderProjectBars(biodiversityProjects)}
                </BarChart>
              </ChartContainer>
          </CardContent>
            </Card>

        <Card className="min-w-0">
          <CardHeader>
            <CardTitle>Annotations in Beta per year per project</CardTitle>
            <CardDescription>Custom groups</CardDescription>
          </CardHeader>
          <CardContent className="min-w-0">
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
                  {renderProjectBars(customProjects)}
                </BarChart>
              </ChartContainer>
          </CardContent>
        </Card>
      </div>
  )
}
