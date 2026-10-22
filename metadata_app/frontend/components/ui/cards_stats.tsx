"use client"

import React from "react"
import {LabelList, Legend, Line, LineChart, XAxis} from "recharts"

import {Card, CardContent, CardDescription, CardHeader, CardTitle,} from "@/components/ui/card"
import {ChartConfig, ChartContainer, ChartTooltip, ChartTooltipContent} from "@/components/ui/chart"

const chartConfig = {
    assembly_count: {
        label: "Assemblies",
        color: "var(--chart-3)",
    },
    annotation_count: {
        label: "Annotations",
        color: "var(--chart-3)",
    },
} satisfies ChartConfig

type ChartValueKey = "assembly_count" | "annotation_count"

type ChartRow = {
    year: number
    assembly_count?: number
    annotation_count?: number
    trend?: number | null
}

function addTrendLine(data: ChartRow[], valueKey: ChartValueKey) {
    const points = data.filter((row) => typeof row[valueKey] === "number")
    if (points.length < 2) {
        return data.map((row) => ({...row, trend: null}))
    }

    const meanX = points.reduce((sum, row) => sum + row.year, 0) / points.length
    const meanY = points.reduce((sum, row) => sum + Number(row[valueKey]), 0) / points.length
    const denominator = points.reduce((sum, row) => sum + (row.year - meanX) ** 2, 0)
    const slope = denominator === 0 ? 0 : points.reduce(
        (sum, row) => sum + (row.year - meanX) * (Number(row[valueKey]) - meanY),
        0,
    ) / denominator
    const intercept = meanY - slope * meanX

    return data.map((row) => ({...row, trend: slope * row.year + intercept}))
}

export function CardsStats() {
    const [assemblyData, setAssemblyData] = React.useState<ChartRow[]>([])
    const [annotationData, setAnnotationData] = React.useState<ChartRow[]>([])

    React.useEffect(() => {
        const fetchAssemblies = async () => {
            const res = await fetch("/api/home_page/home/assemblies")
            const json = await res.json()
            setAssemblyData(addTrendLine(json, "assembly_count"))
        }

        const fetchAnnotations = async () => {
            const res = await fetch("/api/home_page/home/annotations")
            const json = await res.json()
            setAnnotationData(addTrendLine(json, "annotation_count"))
        }

        fetchAssemblies()
        fetchAnnotations()
    }, [])

    return (
        <div className="grid grid-cols-1 gap-4 sm:grid-cols-2">
            <Card className="dark:bg-secondary">
                <CardHeader>
                    <CardTitle className="text-lg font-bold">Assemblies</CardTitle>
                    <CardDescription className="text-sm font-light mb-2">
                        Cumulative assemblies released by year based on NCBI release date
                    </CardDescription>
                </CardHeader>
                <CardContent className="pb-0">
                    <ChartContainer config={chartConfig} className="h-[160px] w-full">
                        <LineChart
                            data={assemblyData}
                            margin={{top: 20, right: 28, left: 20, bottom: 20}}
                        >
                            <XAxis
                                dataKey="year"
                                tickLine={false}
                                axisLine={false}
                                tickMargin={8}
                                interval={0}
                                padding={{left: 16, right: 16}}
                            />
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
                                dataKey="assembly_count"
                                stroke="var(--chart-3)"
                                dot={{r: 4, fill: "var(--chart-3)", stroke: "none"}}
                                activeDot={{r: 6}}
                                name="Assemblies"
                            >
                                <LabelList position="top" offset={8} className="fill-foreground" fontSize={12}/>
                            </Line>
                            <Line
                                type="monotone"
                                dataKey="trend"
                                stroke="var(--muted-foreground)"
                                strokeDasharray="5 5"
                                dot={false}
                                name="Trend"
                            />
                            <Legend verticalAlign="bottom" height={20}/>
                        </LineChart>
                    </ChartContainer>
                </CardContent>
            </Card>

            <Card className="dark:bg-secondary">
                <CardHeader>
                    <CardTitle className="text-lg font-bold">Annotations</CardTitle>
                    <CardDescription className="text-sm font-light mb-2">
                        Cumulative annotations completed by year based on last status update
                    </CardDescription>
                </CardHeader>
                <CardContent className="pb-0">
                    <ChartContainer config={chartConfig} className="h-[160px] w-full">
                        <LineChart
                            data={annotationData}
                            margin={{top: 20, right: 28, left: 20, bottom: 20}}
                        >
                            <XAxis
                                dataKey="year"
                                tickLine={false}
                                axisLine={false}
                                tickMargin={8}
                                interval={0}
                                padding={{left: 16, right: 16}}
                            />
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
                                stroke="var(--chart-3)"
                                dot={{r: 4, fill: "var(--chart-3)", stroke: "none"}}
                                activeDot={{r: 6}}
                                name="Annotations"
                            >
                                <LabelList position="top" offset={8} className="fill-foreground" fontSize={12}/>
                            </Line>
                            <Line
                                type="monotone"
                                dataKey="trend"
                                stroke="var(--muted-foreground)"
                                strokeDasharray="5 5"
                                dot={false}
                                name="Trend"
                            />
                            <Legend verticalAlign="bottom" height={20}/>
                        </LineChart>
                    </ChartContainer>
                </CardContent>
            </Card>
        </div>
    )
}
