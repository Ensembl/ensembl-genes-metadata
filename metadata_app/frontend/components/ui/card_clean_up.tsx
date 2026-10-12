"use client"

import * as React from "react"
import {
  Code2,
  Copy,
  DatabaseZap,
  Loader2,
  ShieldCheck,
} from "lucide-react"
import { Button } from "@/components/ui/button"
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card"
import { Toaster } from "@/components/ui/sonner"
import { toast } from "sonner"



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

type DatabaseCleanupProps = {
  genebuilder?: string | null
}

export function DatabaseCleanup({ genebuilder: selectedGenebuilder }: DatabaseCleanupProps) {
  const [loading, setLoading] = React.useState(false)
  const [genebuilder, setGenebuilder] = React.useState<string>("")
  const [sqlScript, setSqlScript] = React.useState("")


  React.useEffect(() => {
    if (selectedGenebuilder === undefined) return

    setGenebuilder(selectedGenebuilder ?? "")
    setSqlScript("")
  }, [selectedGenebuilder])

  const copyGeneratedSql = async (script = sqlScript, showAlert = true) => {
    if (!script) return false

    try {
      await copyTextToClipboard(script)
      toast.success("SQL script copied", {
        description: "The cleanup script has been copied to your clipboard.",
      })
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
    <>
    <Card className="min-w-0 overflow-hidden bg-card">
      <CardHeader className="space-y-3 pb-3">
        <div className="flex min-w-0 items-start gap-3">
          <div className="flex size-10 shrink-0 items-center justify-center rounded-full bg-primary/10 text-primary">
            <DatabaseZap className="size-5" />
          </div>
          <div className="min-w-0">
            <CardTitle className="text-base">Database cleanup</CardTitle>
            <CardDescription className="mt-1">
              Generate cleanup SQL for the selected genebuilder.
            </CardDescription>
          </div>
        </div>
      </CardHeader>

      <CardContent className="space-y-3">
        <div className="rounded-md bg-background">
          <div className="flex items-center gap-3">
            <div className="min-w-0 flex-1">
              <p className="text-sm font-medium">SQL cleanup script will check for live databases on all genebuild servers and copy the statements to the clipboard</p>
              <p className="truncate text-xs text-muted-foreground">
                {genebuilder
                  ? `Selected genebuilder: ${genebuilder}`
                  : "Select a genebuilder at the top of the page."}
              </p>
            </div>
          </div>
        </div>

        <Button
          onClick={fetchCleanupData}
          disabled={!genebuilder || loading}
          className="w-full whitespace-normal"
        >
          {loading ? (
            <>
              <Loader2 className="mr-2 h-4 w-4 animate-spin" />
              Generating SQL script
            </>
          ) : sqlScript ? (
            <>
              <Copy className="mr-2 h-4 w-4" />
              Copy again
            </>
          ) : (
            <>
              <Code2 className="mr-2 h-4 w-4" />
              Generate SQL script
            </>
          )}
        </Button>

        <div className="flex items-center justify-center gap-2 text-xs text-muted-foreground">
          <ShieldCheck className="size-3.5" />
          <span>Always review generated scripts before execution. Anno pipe DBs should be reviewed manually as they are not listed here.</span>
        </div>
      </CardContent>
    </Card>
    <Toaster position="bottom-right" richColors />
    </>
  )
}
