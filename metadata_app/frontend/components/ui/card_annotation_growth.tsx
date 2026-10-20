"use client"

import * as React from "react"
import {Bar, BarChart, CartesianGrid, Legend, XAxis, YAxis} from "recharts"
import {Download} from "lucide-react"
import {Button} from "@/components/ui/button"
import {Checkbox} from "@/components/ui/checkbox"
import {Card, CardContent, CardDescription, CardHeader, CardTitle} from "@/components/ui/card"
import {type ChartConfig, ChartContainer, ChartTooltip, ChartTooltipContent} from "@/components/ui/chart"

type GrowthRow = {
    quarter: string
    total_genomes: number
    full_genebuild: number
    anno: number
    braker: number
    helixer: number
}

const methods = [
    ["full_genebuild", "Full Genebuild", "var(--chart-5)"],
    ["braker", "Braker", "var(--chart-7)"],
    ["helixer", "Helixer", "var(--chart-8)"],
    ["anno", "Anno", "var(--chart-6)"]
] as const

const chartConfig = {
    total_genomes: {label: "All methods", color: "var(--chart-5)"},
    full_genebuild: {label: "Full Genebuild", color: "var(--chart-5)"},
    braker: {label: "Braker", color: "var(--chart-7)"},
    helixer: {label: "Helixer", color: "var(--chart-8)"},
    anno: {label: "Anno", color: "var(--chart-6)"},
} satisfies ChartConfig

const renderMethodBars = (showMethods: boolean) =>
    showMethods ? (
        methods.map(([key, label, color], index) => (
            <Bar
                key={key}
                dataKey={key}
                name={label}
                fill={color}
                radius={
                    index === 0
                        ? [0, 0, 4, 4]
                        : index === methods.length - 1
                            ? [4, 4, 0, 0]
                            : [0, 0, 0, 0]
                }
                stackId="methods"
            />
        ))
    ) : (
        <Bar dataKey="total_genomes" name="All methods" fill="var(--chart-5)" radius={4}/>
    )

function downloadCsv(data: GrowthRow[]) {
    const headers = ["quarter", "total_genomes", "full_genebuild", "anno", "braker", "helixer"]
    const csv = [headers, ...data.map((row) => headers.map((header) => row[header as keyof GrowthRow]))]
        .map((row) => row.join(","))
        .join("\n")
    const url = URL.createObjectURL(new Blob([csv], {type: "text/csv;charset=utf-8"}))
    const link = document.createElement("a")
    link.href = url
    link.download = "ensembl_annotation_growth_last_4_years.csv"
    link.click()
    URL.revokeObjectURL(url)
}

export function AnnotationGrowthChart() {
    const [data, setData] = React.useState<GrowthRow[]>([])
    const [groupByMethod, setGroupByMethod] = React.useState(false)

    React.useEffect(() => {
        fetch("/api/home_page/home/annotation_growth")
            .then((response) => response.json())
            .then(setData)
            .catch((error) => console.error("Error fetching annotation growth:", error))
    }, [])

    return (
        <Card className="dark:bg-secondary">
            <CardHeader>
                <div className="flex flex-wrap items-start justify-between gap-3">
                    <div>
                        <CardTitle>Ensembl annotation growth</CardTitle>
                        <CardDescription>
                            Live annotated genomes to date by release date, including the pre-window baseline
                            and quarterly additions
                        </CardDescription>
                    </div>
                    <div className="flex flex-wrap items-center gap-2">
                        <label className="flex items-center gap-2 text-sm">
                            <Checkbox
                                checked={groupByMethod}
                                onCheckedChange={(checked) => setGroupByMethod(checked === true)}
                                aria-label="Show annotation methods"
                            />
                            Show methods
                        </label>
                        <Button variant="default" size="sm" onClick={() => downloadCsv(data)} disabled={!data.length}>
                            <Download/> Download CSV
                        </Button>
                    </div>
                </div>
            </CardHeader>
            <CardContent>
                <ChartContainer config={chartConfig} className="h-[360px] w-full">
                    <BarChart accessibilityLayer data={data} margin={{top: 16, right: 16, left: 8, bottom: 8}}>
                        <CartesianGrid vertical={false}/>
                        <XAxis dataKey="quarter" tickLine={false} axisLine={false} angle={-35} textAnchor="end"
                               height={70}/>
                        <YAxis allowDecimals={false} tickLine={false} axisLine={false} width={55}/>
                        <ChartTooltip content={<ChartTooltipContent/>}/>
                        <Legend/>
                        {renderMethodBars(groupByMethod)}
                    </BarChart>
                </ChartContainer>

                <div className="mt-4 overflow-x-auto">
                    <table className="w-full min-w-[620px] text-sm">
                        <caption className="sr-only">Annotated genomes per quarter by method</caption>
                        <thead>
                        <tr className="border-b text-left">
                            <th className="px-2 py-2">Quarter</th>
                            <th className="px-2 py-2 text-right">All methods</th>
                            {methods.map(([key, label]) => <th key={key} className="px-2 py-2 text-right">{label}</th>)}
                        </tr>
                        </thead>
                        <tbody>
                        {data.map((row) => (
                            <tr key={row.quarter} className="border-b last:border-0">
                                <td className="px-2 py-1.5">{row.quarter}</td>
                                <td className="px-2 py-1.5 text-right font-medium">{row.total_genomes.toLocaleString()}</td>
                                {methods.map(([key]) => <td key={key}
                                                            className="px-2 py-1.5 text-right">{row[key as keyof GrowthRow].toLocaleString()}</td>)}
                            </tr>
                        ))}
                        </tbody>
                    </table>
                </div>
            </CardContent>
        </Card>
    )
}
