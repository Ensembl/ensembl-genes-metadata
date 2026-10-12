"use client"

import * as React from "react"
import { Bar, BarChart, CartesianGrid, LabelList, XAxis } from "recharts"
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
import {
  DropdownMenu,
  DropdownMenuCheckboxItem,
  DropdownMenuContent,
  DropdownMenuTrigger,
} from "@/components/ui/dropdown-menu"
import { Button } from "@/components/ui/button"

export type StatusItem = {
  gb_status: string
  count: number
}

type Props = {
  data: StatusItem[]
}

const chartConfig: Record<string, { label: string; color?: string }> = {
  count: { label: "Annotations" },
  in_progress: { label: "In progress", color: "var(--chart-3)" },
  completed: { label: "Completed", color: "var(--chart-3)" },
  pre_released: { label: "Pre-released", color: "var(--chart-3)" },
  handed_over: { label: "Handed over", color: "var(--chart-3)" },
  coming_soon: { label: "Almost live", color: "var(--chart-3)" },
  live: { label: "Live", color: "var(--chart-3)" },
  archive: { label: "Archived", color: "var(--chart-3)" },
  check_busco: { label: "Low pBUSCO", color: "var(--chart-3)" },
  insufficient_data: { label: "Low evidence", color: "var(--chart-3)" },
  poor_genome_busco: { label: "Low gBUSCO", color: "var(--chart-3)" },
  abandoned: { label: "Abandoned", color: "var(--chart-3)" },
} satisfies ChartConfig

const statusOrder = Object.keys(chartConfig).filter((key) => key !== "count")

export function RepStatus({ data }: Props) {
  const [visibleStatuses, setVisibleStatuses] = React.useState<string[]>(statusOrder)

  const transformedData = React.useMemo(() => {
    const orderMap = new Map(statusOrder.map((key, index) => [key, index]))

    return [...data]
      .filter((item) => visibleStatuses.includes(item.gb_status))
      .sort(
        (a, b) =>
          (orderMap.get(a.gb_status) ?? Infinity) -
          (orderMap.get(b.gb_status) ?? Infinity)
      )
      .map((item) => ({
        ...item,
        displayName: chartConfig[item.gb_status]?.label || item.gb_status,
        fill: chartConfig[item.gb_status]?.color,
      }))
  }, [data, visibleStatuses])

  return (
    <Card>
      <CardHeader className="flex flex-row items-start justify-between gap-4">
        <div>
          <CardTitle>Annotation status</CardTitle>
          <CardDescription>Current status of annotations</CardDescription>
        </div>

        <DropdownMenu>
          <DropdownMenuTrigger asChild>
            <Button variant="outline">Show bars</Button>
          </DropdownMenuTrigger>
          <DropdownMenuContent align="end" className="max-h-80 overflow-auto">
            {statusOrder.map((statusKey) => (
              <DropdownMenuCheckboxItem
                key={statusKey}
                checked={visibleStatuses.includes(statusKey)}
                onCheckedChange={(checked) => {
                  setVisibleStatuses((current) =>
                    checked
                      ? [...current, statusKey]
                      : current.filter((s) => s !== statusKey)
                  )
                }}
              >
                {chartConfig[statusKey].label}
              </DropdownMenuCheckboxItem>
            ))}
          </DropdownMenuContent>
        </DropdownMenu>
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
            <ChartTooltip cursor={false} content={<ChartTooltipContent />} />
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