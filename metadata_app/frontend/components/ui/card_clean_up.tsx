"use client"

import * as React from "react"
import { ClipboardCheck, Loader2 } from "lucide-react"
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

async function copyTextToClipboard(text: string) {
  try {
    if (navigator.clipboard && window.isSecureContext) {
      await navigator.clipboard.writeText(text)
      return
    }
  } catch (err) {
    console.warn("navigator.clipboard failed, trying fallback copy", err)
  }

  const textArea = document.createElement("textarea")
  textArea.value = text
  textArea.style.position = "fixed"
  textArea.style.left = "-9999px"
  textArea.style.top = "0"
  textArea.setAttribute("readonly", "")

  document.body.appendChild(textArea)
  textArea.focus()
  textArea.select()

  try {
    const copied = document.execCommand("copy")
    if (!copied) {
      throw new Error("Fallback copy command was rejected")
    }
  } finally {
    document.body.removeChild(textArea)
  }
}

export function DatabaseCleanup() {
  const [loading, setLoading] = React.useState(false)
  const [genebuilder, setGenebuilder] = React.useState<string>("")
  const [copied, setCopied] = React.useState(false)
  const [sqlScript, setSqlScript] = React.useState("")

  const genebuilderOptions = [
    { value: "lazar", label: "Anna" },
    { value: "leanne", label: "Leanne" },
    { value: "jackt", label: "Jack" },
    { value: "vianey", label: "Vianey" },
    { value: "swati", label: "Swati" },
    { value: "ereboperezsilva", label: "Jose" },
    { value: "ftricomi", label: "Francesca" },
  ]

  const copyGeneratedSql = async (script = sqlScript, showAlert = true) => {
    if (!script) return false

    try {
      await copyTextToClipboard(script)
      setCopied(true)
      return true
    } catch (err) {
      console.error("Clipboard error:", err)
      if (showAlert) {
        alert("Failed to copy SQL script to clipboard. See console for details.")
      }
      return false
    }
  }

  const fetchCleanupData = async () => {
    if (!genebuilder) return
    if (sqlScript) {
      await copyGeneratedSql()
      return
    }

    setLoading(true)
    setCopied(false)
    setSqlScript("")

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

      const generatedSqlScript = await response.text()
      setSqlScript(generatedSqlScript)

      await copyGeneratedSql(generatedSqlScript, false)
    } catch (err) {
      console.error("SQL generation error:", err)
      alert("Failed to generate SQL script. See console for details.")
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
          copy the SQL script to your clipboard. Anno pipe DBs need to be checked manually as they are not listed here.
        </CardDescription>
      </CardHeader>

      <CardContent className="space-y-6">
        <div className="grid grid-cols-2 items-center gap-4">
          <Select
            onValueChange={(value) => {
              setGenebuilder(value)
              setCopied(false)
              setSqlScript("")
            }}
          >
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
            variant={copied ? "ghost" : "default"}
          >
            {loading ? (
              <>
                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                Generating SQL script...
              </>
            ) : copied ? (
              <>
                <ClipboardCheck className="mr-2 h-4 w-4" />
                Copied
              </>
            ) : sqlScript ? (
              "Copy generated SQL"
            ) : (
              "Generate SQL script"
            )}
          </Button>
        </div>
      </CardContent>
    </Card>
  )
}
