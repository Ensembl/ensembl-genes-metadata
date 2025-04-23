"use client"

import { Tooltip, Line, LineChart } from "recharts"


import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
  CardDescription,
} from "@/components/ui/card"
import { ChartConfig, ChartContainer } from "@/components/ui/chart"

const data = [
  {
    revenue: 10400,
    subscription: 240,
  },
  {
    revenue: 14405,
    subscription: 300,
  },
  {
    revenue: 9400,
    subscription: 200,
  },
  {
    revenue: 8200,
    subscription: 278,
  },
  {
    revenue: 7000,
    subscription: 189,
  },
  {
    revenue: 9600,
    subscription: 239,
  },
  {
    revenue: 11244,
    subscription: 278,
  },
  {
    revenue: 26475,
    subscription: 189,
  },
]

const chartConfig = {
  revenue: {
    label: "Revenue",
    color: "hsl(var(--primary))",
  },
  subscription: {
    label: "Subscriptions",
    color: "hsl(var(--primary))",
  },
} satisfies ChartConfig

export function CardsStats() {
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
              data={data}
              margin={{
                top: 5,
                right: 10,
                left: 10,
                bottom: 0,
              }}
            >
              <Tooltip
                cursor={{ fill: 'rgba(0, 0, 0, 0.1)' }}
                contentStyle={{ backgroundColor: 'var(--background)',
                  borderRadius: '0.5rem',
                  border: '1px solid var(--border)',
                  boxShadow: '0 4px 12px rgba(0, 0, 0, 0.15)'
                }}
                labelStyle={{ fontWeight: 'bold', color: 'var(--chart-1)' }}
                itemStyle={{ color: 'var(--chart-1)' }}
              />
              <Line
                type="monotone"
                strokeWidth={2}
                dataKey="revenue"
                stroke="var(--chart-1)"
                activeDot={{
                  r: 6,
                }}
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
              data={data}
              margin={{
                top: 5,
                right: 10,
                left: 10,
                bottom: 0,
              }}
            >
              <Tooltip
                cursor={{ fill: 'rgba(0, 0, 0, 0.1)' }}
                contentStyle={{ backgroundColor: 'var(--background)',
                  borderRadius: '0.5rem',
                  border: '1px solid var(--border)',
                  boxShadow: '0 4px 12px rgba(0, 0, 0, 0.15)'
                }}
                labelStyle={{ fontWeight: 'bold', color: 'var(--chart-1)' }}
                itemStyle={{ color: 'var(--chart-1)' }}
              />
              <Line
                type="monotone"
                strokeWidth={2}
                dataKey="revenue"
                stroke="var(--chart-1)"
                activeDot={{
                  r: 6,
                }}
              />
            </LineChart>
          </ChartContainer>
        </CardContent>
      </Card>
    </div>
  )
}