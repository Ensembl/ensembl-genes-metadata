"use client"

import * as React from "react"
import {Download} from "lucide-react"
import {CartesianGrid, Line, LineChart, Tooltip, XAxis, YAxis} from "recharts"
import {Button} from "@/components/ui/button"
import {Card, CardContent, CardDescription, CardHeader, CardTitle} from "@/components/ui/card"
import {ChartContainer} from "@/components/ui/chart"

type SummaryRow = {
    annotation_method: string
    sample_size: number
    average_days: number | null
    median_days: number | null
    minimum_days: number | null
    maximum_days: number | null
}

type TrendRow = { month: string; [key: string]: string | number | null }

const methods = [
    ["anno", "Anno", "var(--chart-6)"],
    ["full_genebuild", "Full genebuild", "var(--chart-5)"],
] as const

function formatDays(value: number | null | undefined) {
    return value == null ? "No data" : `${value.toLocaleString()} days`
}

function downloadCsv(trend: TrendRow[]) {
    const headers = [
        "month",
        ...methods.flatMap(([method]) => [
            `average_days_${method}`,
            `median_days_${method}`,
            `lower_quartile_days_${method}`,
            `upper_quartile_days_${method}`,
            `sample_size_${method}`,
        ]),
    ]
    const csv = [headers, ...trend.map((row) => headers.map((header) => row[header] ?? ""))]
        .map((row) => row.join(","))
        .join("\n")
    const url = URL.createObjectURL(new Blob([csv], {type: "text/csv;charset=utf-8"}))
    const link = document.createElement("a")
    link.href = url
    link.download = "annotation_duration_history.csv"
    document.body.appendChild(link)
    link.click()
    document.body.removeChild(link)
    URL.revokeObjectURL(url)
}

export function AnnotationDurationChart() {
    const [summary, setSummary] = React.useState<SummaryRow[]>([])
    const [trend, setTrend] = React.useState<TrendRow[]>([])
    const [windowLabel, setWindowLabel] = React.useState("latest available 3 months")

    React.useEffect(() => {
        fetch("/api/home_page/home/annotation_time")
            .then((response) => {
                if (!response.ok) throw new Error(`Server error: ${response.status}`)
                return response.json()
            })
            .then((result) => {
                setSummary(result.summary ?? [])
                setTrend(result.trend ?? [])
                if (result.window_start && result.window_end) {
                    setWindowLabel(`${result.window_start} to ${result.window_end}`)
                }
            })
            .catch((error) => console.error("Error fetching annotation durations:", error))
    }, [])

    return (
        <Card className="dark:bg-secondary">
            <CardHeader>
                <div className="flex items-center justify-between gap-3">
                    <CardTitle>Annotation time</CardTitle>
                    <Button variant="default" size="sm" onClick={() => downloadCsv(trend)} disabled={!trend.length}>
                        <Download/> Download history
                    </Button>
                </div>
                <CardDescription>
                    Start to last geneset update. Cards show the {windowLabel}; the chart shows monthly
                    medians across all available history. Tooltips include the mean, middle 50% range, and
                    sample size.
                </CardDescription>
            </CardHeader>
            <CardContent>
                <div className="grid gap-4 md:grid-cols-2">
                    {methods.map(([method, label, color]) => {
                        const item = summary.find((row) => row.annotation_method === method)
                        return (
                            <div key={method} className="rounded-lg border border-t-4 p-4"
                                 style={{borderTopColor: color}}>
                                <p className="text-sm font-medium text-muted-foreground">{label}</p>
                                <p className="mt-2 text-xs text-muted-foreground">Median duration</p>
                                <p className="text-2xl font-bold">{formatDays(item?.median_days)}</p>
                                <p className="mt-1 text-xs text-muted-foreground">
                                    {item?.sample_size ?? 0} completed genome{item?.sample_size === 1 ? "" : "s"}
                                </p>
                                <p className="mt-3 text-xs text-muted-foreground">
                                    Mean: {formatDays(item?.average_days)} ·
                                    Range: {formatDays(item?.minimum_days)}–{formatDays(item?.maximum_days)}
                                </p>
                            </div>
                        )
                    })}
                </div>

                {trend.length > 0 ? (
                    <>
                        <ChartContainer config={{}} className="mt-8 h-[360px] w-full">
                            <LineChart accessibilityLayer data={trend}
                                       margin={{top: 16, right: 16, left: 8, bottom: 8}}>
                                <CartesianGrid vertical={false}/>
                                <XAxis dataKey="month" tickLine={false} axisLine={false}/>
                                <YAxis allowDecimals={false} tickLine={false} axisLine={false} width={60}/>
                                <Tooltip
                                    content={({active, label, payload}) => {
                                        if (!active || !payload?.length) return null
                                        const row = payload[0].payload as TrendRow

                                        return (
                                            <div className="rounded-lg border bg-background p-3 text-sm shadow-sm">
                                                <p className="font-medium">{label}</p>
                                                {methods.map(([method, methodLabel]) => {
                                                    const median = row[`median_days_${method}`]
                                                    const average = row[`average_days_${method}`]
                                                    const lower = row[`lower_quartile_days_${method}`]
                                                    const upper = row[`upper_quartile_days_${method}`]
                                                    const sampleSize = row[`sample_size_${method}`]

                                                    if (sampleSize === 0 || median == null) return null

                                                    return (
                                                        <div key={method} className="mt-2">
                                                            <p className="font-medium">{methodLabel}</p>
                                                            <p>Median: {Number(median).toLocaleString()} days</p>
                                                            <p>Mean: {Number(average).toLocaleString()} days</p>
                                                            <p>
                                                                Middle 50%: {Number(lower).toLocaleString()}–
                                                                {Number(upper).toLocaleString()} days
                                                            </p>
                                                            <p>Completed: {sampleSize}</p>
                                                        </div>
                                                    )
                                                })}
                                            </div>
                                        )
                                    }}
                                />
                                {methods.map(([method, label, color]) => (
                                    <Line key={method} type="monotone" dataKey={`median_days_${method}`} name={label}
                                          stroke={color} connectNulls dot={false}/>
                                ))}
                            </LineChart>
                        </ChartContainer>
                    </>
                ) : (
                    <p className="mt-8 text-center text-muted-foreground">No completed annotation durations are
                        available.</p>
                )}
            </CardContent>
        </Card>
    )
}
