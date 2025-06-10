"use client"
import * as React from "react"
import Link from "next/link"
import {
  Card,
  CardHeader,
  CardTitle,
  CardContent,
  CardFooter,
} from "@/components/ui/card"
import { Button } from "@/components/ui/button"

export function ReportCard() {
  return (
    <Card className="bg-primary text-primary-foreground">
      <CardHeader>
        <CardTitle className="text-primary-foreground">
          Generate reports
        </CardTitle>
      </CardHeader>

      <CardContent className="text-muted-foreground text-sm">
        Generate reports on available annotations by Genebuild or identify assemblies ready for annotation. Create and download publication-ready tables and visualizations.
      </CardContent>

      <CardFooter className="justify-end">
        <div className="flex gap-4">
          <Button variant="secondary">
          <Link href="/report/asm">
            Assembly report <span className="ml-1">→</span>
          </Link>
        </Button>
          <Button variant="secondary">
          <Link href="/report/anno">
            Annotation report <span className="ml-1">→</span>
          </Link>
        </Button>
            </div>
      </CardFooter>
    </Card>
  )
}