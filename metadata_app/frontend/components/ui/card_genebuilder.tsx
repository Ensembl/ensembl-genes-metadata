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

// Map API genebuilder names to display names
const nameMapping: Record<string, string> = {
  lazar: "Anna",
  jose: "Jose",
  ftricomi: "Francesca",
  vianey: "Vianey",
  swati: "Swati",
  jackt: "Jack",
}

// Define type for chart config entries
interface ChartConfigEntry {
  label: string;
  color?: string;
}

const chartConfig: Record<string, ChartConfigEntry> = {
  annotations: {
    label: "Annotations",
  },
  anna: {
    label: "Anna",
    color: "var(--chart-1)",
  },
  jose: {
    label: "Jose",
    color: "var(--chart-2)",
  },
  francesca: {
    label: "Francesca",
    color: "var(--chart-3)",
  },
  vianey: {
    label: "Vianey",
    color: "var(--chart-4)",
  },
  swati: {
    label: "Swati",
    color: "var(--chart-5)",
  },
  jack: {
    label: "Jack",
    color: "var(--chart-6)",
  },
} satisfies ChartConfig

  // Define types for our data
  interface GenebuilderItem {
    genebuilder: string;
    annotations: number;
    [key: string]: any;
  }

  interface TransformedItem extends GenebuilderItem {
    displayName: string;
    fill?: string; // Add fill property for color
  }

export function CardGeneBuilder() {
  const [genebuilderData, setGenebuilderData] = React.useState<GenebuilderItem[]>([])
  const [totalAnnotations, setTotalAnnotations] = React.useState(0)

  React.useEffect(() => {
    const fetchGenebuilderData = async () => {
      try {
        const res = await fetch("/api/genebuilder")
        const json = await res.json()

        // Define the structure of the genebuilder data
        interface GenebuilderItem {
          genebuilder: string;
          annotations: number;
          [key: string]: any;
        }

        interface TransformedItem extends GenebuilderItem {
          displayName: string;
        }

        // Transform the data to use display names instead of API names
        // Only include genebuilders who are in our nameMapping list
        const transformedData = json
          .filter((item: GenebuilderItem) => nameMapping[item.genebuilder] !== undefined)
          .map((item: GenebuilderItem): TransformedItem => {
            const displayName = nameMapping[item.genebuilder];
            // Get the corresponding chart config key (lowercase name)
            const configKey = displayName.toLowerCase();
            // Get the color from chartConfig if available
            const configEntry = chartConfig[configKey as keyof typeof chartConfig];
            const color = configEntry && 'color' in configEntry ? configEntry.color : "";

            return {
              ...item,
              displayName,
              fill: color // Add the color directly to each data item
            };
          });

        setGenebuilderData(transformedData)

        // Calculate total annotations by summing all annotation counts
        const total = transformedData.reduce((sum: number, item: TransformedItem) => sum + (item.annotations || 0), 0)
        setTotalAnnotations(total)
      } catch (error) {
        console.error("Error fetching genebuilder data:", error)
        setGenebuilderData([])
        setTotalAnnotations(0)
      }
    }

    fetchGenebuilderData()
  }, [])

  return (
    <Card className="flex flex-col">
      <CardHeader className="items-center pb-0">
        <CardTitle>Top Genebuilders</CardTitle>
        <CardDescription>Last 30 days</CardDescription>
      </CardHeader>
      <CardContent className="flex-1 pb-0">
        <ChartContainer
          config={chartConfig}
          className="mx-auto aspect-square max-h-[250px]"
        >
          <PieChart>
            <ChartTooltip
              cursor={false}
              content={<ChartTooltipContent hideLabel />}
            />
            <Pie
              data={genebuilderData}
              dataKey="annotations"
              nameKey="displayName"
              innerRadius={60}
              strokeWidth={5}
              // No need to specify cell colors - they come from each data point's fill property
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
                  return null;
                }}
              />
            </Pie>
          </PieChart>
        </ChartContainer>
      </CardContent>
    </Card>
  )
}