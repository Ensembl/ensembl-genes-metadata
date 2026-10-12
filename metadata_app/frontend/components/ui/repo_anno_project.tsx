"use client"

import { Bar, BarChart, CartesianGrid, LabelList, XAxis } from "recharts"
import { PROJECTS } from "@/features/projects/project-config"
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
import * as React from "react"

export type ProjectItem = {
  associated_project: string
  count: number
}
type Props = {
  data: ProjectItem[]
  title?: string
  description?: string
}

const biodiversityProjects = PROJECTS.filter(
  (project) => project.reportGroup === "biodiversity",
)

const chartConfig = {
  count: { label: "Annotations" },
  ...Object.fromEntries(
    biodiversityProjects.map((project) => [
      project.reportKey,
      {
        label: project.reportLabel ?? project.reportKey,
        color: "var(--chart-1)",
      },
    ]),
  ),
} satisfies ChartConfig

export function RepProject({
  data,
  title = "Associated biodiversity projects",
  description = "Number of annotations per project",
}: Props) {
  const transformedData = React.useMemo(() => {
    const countByProject = new Map(
      data
        .filter((item) => item.associated_project && item.associated_project !== "count")
        .map((item) => [item.associated_project.trim(), item.count]),
    )
    const knownProjectRows = biodiversityProjects.map((project) => ({
      associated_project: project.reportKey,
      count: countByProject.get(project.reportKey) ?? 0,
      displayName: project.reportLabel ?? project.reportKey,
    }))
    const unknownProjectRows = data
      .filter(
        (item) =>
          item.associated_project &&
          item.associated_project.trim() !== "count" &&
          !chartConfig[item.associated_project.trim()],
      )
      .map((item) => ({
        ...item,
        associated_project: item.associated_project.trim(),
        displayName: item.associated_project.trim(),
      }))

    return [...knownProjectRows, ...unknownProjectRows].filter((item) => item.count > 0)
  }, [data])

  return (
    <Card>
      <CardHeader>
        <CardTitle>{title}</CardTitle>
        <CardDescription>{description}</CardDescription>
      </CardHeader>
      <CardContent>
        <ChartContainer config={chartConfig}>
          <BarChart
            accessibilityLayer
            data={transformedData}
            margin={{ top: 20, bottom: 80 }}
          >
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
            <Bar dataKey="count" name="Annotations" fill="var(--color-chart-1)" radius={8}>
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
