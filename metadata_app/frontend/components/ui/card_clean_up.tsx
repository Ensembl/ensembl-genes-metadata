"use client"

import * as React from "react"
import { Download, Loader2 } from "lucide-react"
import { Button } from "@/components/ui/button"
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"

import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select"

export function DatabaseCleanup() {
  const [loading, setLoading] = React.useState(false)
  const [genebuilder, setGenebuilder] = React.useState<string>("")
  const [downloaded, setDownloaded] = React.useState(false)

  const genebuilderOptions = [
    { value: "lazar", label: "Anna" },
    { value: "leanne", label: "Leanne" },
    { value: "jackt", label: "Jack" },
    { value: "vianey", label: "Vianey" },
    { value: "swati", label: "Swati" },
    { value: "ereboperezsilva", label: "Jose" },
    { value: "ftricomi", label: "Francesca" },
  ]

const fetchCleanupData = async () => {
  if (!genebuilder) return

  setLoading(true)
  setDownloaded(false)

  try {
    const formData = new FormData()
    formData.append("genebuilder", genebuilder)

    const response = await fetch("/api/clean/db_clean/genebuilder", {
      method: "POST",
      body: formData,
    })

    if (!response.ok) {
      throw new Error(`Failed to generate SQL script: ${response.statusText}`)
    }

    const blob = await response.blob()
    const url = window.URL.createObjectURL(blob)
    const a = document.createElement("a")
    a.href = url
    a.download = `cleanup_${genebuilder}.sql`
    document.body.appendChild(a)
    a.click()
    a.remove()
    window.URL.revokeObjectURL(url)

    setDownloaded(true)
  } catch (err) {
    console.error("Download error:", err)
    alert("Failed to download SQL script. See console for details.")
  } finally {
    setLoading(false)
  }
}

  return (
    <Card className="dark:bg-secondary">
      <CardHeader>
        <CardTitle className="text-lg">Database cleanup</CardTitle>
        <CardDescription>
          Select a genebuilder to list databases eligible for cleanup and
          download SQL script. Anno pipe DBs need to be checked manually as they are not listed here.
        </CardDescription>
      </CardHeader>

      <CardContent className="space-y-6">
        <div className="grid grid-cols-2 items-center gap-4">
          <Select onValueChange={setGenebuilder}>
            <SelectTrigger className="w-full">
              <SelectValue placeholder="Select genebuilder" />
            </SelectTrigger>
            <SelectContent>
              {genebuilderOptions.map((gb) => (
                <SelectItem key={gb.value} value={gb.value}>
                  {gb.label}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>

          <Button
            onClick={fetchCleanupData}
            disabled={!genebuilder || loading}
            variant={downloaded ? "ghost" : "default"}
          >
            {loading ? (
              <>
                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                Generating SQL script...
              </>
            ) : downloaded ? (
              <>
                <Download className="mr-2 h-4 w-4" />
                Downloaded
              </>
            ) : (
              "Generate SQL script"
            )}
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}