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
import { Input } from "@/components/ui/input"
import { Button } from "@/components/ui/button"
import { Label as FormLabel } from "@/components/ui/label"

// Mapping API names to display names
const nameMapping: Record<string, string> = {
  lazar: "Anna",
  jose: "Jose",
  ftricomi: "Francesca",
  vianey: "Vianey",
  swati: "Swati",
  jackt: "Jack",
}

// Define chart config entries
const chartConfig: Record<string, { label: string; color?: string }> = {
  annotations: { label: "Annotations" },
  anno: { label: "Anno", color: "var(--chart-1)" },
  braker: { label: "Braker", color: "var(--chart-2)" },
  external_annotation_import: { label: "External Annotation Import", color: "var(--chart-3)" },
  full_genebuild: { label: "Full Genebuild", color: "var(--chart-4)" },
  import: { label: "Import", color: "var(--chart-5)" },
  mixed_strategy_build: { label: "Mixed Strategy Build", color: "var(--chart-6)" },
  projection_build: { label: "Projection Build", color: "var(--chart-7)" },
}

interface MethodItem {
  annotation_method: string
  number_of_annotations: number
  [key: string]: any
}

interface TransformedItem extends MethodItem {
  displayName: string
  fill?: string
}

interface FilterData {
  bioproject_id?: string[]
  release_type?: string[]
  release_date?: string
  taxon_id?: number
}

export function AnoMethodSummaryChart() {
  const [methodData, setMethodData] = React.useState<TransformedItem[]>([])
  const [totalAnnotations, setTotalAnnotations] = React.useState(0)
  const [filters, setFilters] = React.useState({
    bioproject_id: [""],
    release_date: "",
    taxon_id: 0,
    release_type: [""],
  })

const fetchMethodData = async (filtersData: FilterData = {}) => {
    const isEmpty = Object.keys(filtersData).length === 0

  const cleanedFilters = isEmpty
    ? {}
    : {
        bioproject_id: filtersData.bioproject_id?.filter((id: string) => id.trim() !== "") ?? [],
        release_type: filtersData.release_type?.filter((r: string) => r.trim() !== "") ?? [],
        release_date: filtersData.release_date || undefined,
        taxon_id: filtersData.taxon_id || 0,
      }

  try {
    const res = await fetch("http://127.0.0.1:8000/api/home_page/home/method_summary", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(cleanedFilters),
    })

    if (!res.ok) throw new Error(`Server error: ${res.status}`)

    const json = await res.json()

    if (!Array.isArray(json)) throw new Error("API did not return an array")

    const transformedData = json.map((item: MethodItem): TransformedItem => {
      const displayName = nameMapping[item.annotation_method] || item.annotation_method
      const configKey = displayName.toLowerCase()
      const configEntry = chartConfig[configKey as keyof typeof chartConfig]
      const color = configEntry?.color ?? ""

      return {
        ...item,
        displayName,
        fill: color,
      }
    })

    setMethodData(transformedData)

    const total = transformedData.reduce(
      (sum, item) => sum + (item.number_of_annotations || 0),
      0
    )
    setTotalAnnotations(total)
  } catch (error) {
    console.error("Error fetching data:", error)
    setMethodData([])
    setTotalAnnotations(0)
  }
}

  React.useEffect(() => {
    fetchMethodData() // initial load with empty filters
  }, [])

  const handleChange = (field: string, value: string) => {
    setFilters((prev) => ({
      ...prev,
      [field]: field === "taxon_id" ? Number(value) : field === "bioproject_id" || field === "release_type" ? [value] : value,
    }))
  }

  return (
    <Card className="flex flex-col">
      <CardHeader className="pb-0">
        <CardTitle>Annotation method summary</CardTitle>
        <CardDescription>Use the filters below to refine the data</CardDescription>

        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mt-4">
          <div>
            <FormLabel htmlFor="bioproject_id" className="mb-2 block">BioProject ID</FormLabel>
            <Input
              id="bioproject_id"
              placeholder="PRJNA123456"
              value={filters.bioproject_id[0]}
              onChange={(e) => handleChange("bioproject_id", e.target.value)}
            />
          </div>

          <div>
            <FormLabel htmlFor="release_date" className="mb-2 block">Release Date</FormLabel>
            <Input
              id="release_date"
              placeholder="2024-03-05"
              value={filters.release_date}
              onChange={(e) => handleChange("release_date", e.target.value)}
            />
          </div>

          <div>
            <FormLabel htmlFor="taxon_id" className="mb-2 block">Taxon ID</FormLabel>
            <Input
              id="taxon_id"
              type="number"
              placeholder="9606"
              value={filters.taxon_id || ""}
              onChange={(e) => handleChange("taxon_id", e.target.value)}
            />
          </div>

          <div>
            <FormLabel htmlFor="release_type" className="mb-2 block">Release Site</FormLabel>
            <Input
              id="release_type"
              placeholder="main or beta"
              value={filters.release_type[0]}
              onChange={(e) => handleChange("release_type", e.target.value)}
            />
          </div>
        </div>
        <Button className="mt-2" onClick={() => fetchMethodData(filters)}>
          Apply Filters
        </Button>
      </CardHeader>

      <CardContent className="flex-1 pb-0">
        {methodData.length > 0 ? (
          <ChartContainer
            config={chartConfig}
            className="mx-auto aspect-square max-h-[250px]"
          >
            <PieChart>
              <ChartTooltip cursor={false} content={<ChartTooltipContent hideLabel />} />
              <Pie
                data={methodData}
                dataKey="number_of_annotations"
                nameKey="displayName"
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
                    return null
                  }}
                />
              </Pie>
            </PieChart>
          </ChartContainer>
        ) : (
          <div className="text-center text-muted-foreground py-12">
            No annotations match the filters.
          </div>
        )}
      </CardContent>
    </Card>
  )
}