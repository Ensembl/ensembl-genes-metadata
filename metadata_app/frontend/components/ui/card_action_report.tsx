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
import { BackgroundGradient } from "@/components/ui/backround-gradient"



export function ReportCard() {
  return (
      <BackgroundGradient>
    <Card className="bg-primary border-none">
      <CardHeader>
        <CardTitle className="text-primary-foreground text-lg">
          Generate reports
        </CardTitle>
      </CardHeader>

      <CardContent className="text-sidebar-ring text-sm">
        Generate reports on available annotations by Genebuild or identify assemblies ready for annotation. Create and download publication-ready tables and visualizations.
      </CardContent>

      <CardFooter className="justify-end">
        <div className="flex w-full flex-col gap-3 sm:w-auto sm:flex-row sm:gap-4">
          <Button variant="secondary" className="w-full sm:w-auto">
          <Link href="/report/asm">
            Assembly report <span className="ml-1">→</span>
          </Link>
        </Button>
          <Button variant="secondary" className="w-full sm:w-auto">
          <Link href="/report/anno">
            Annotation report <span className="ml-1">→</span>
          </Link>
        </Button>
            </div>
      </CardFooter>
    </Card>
    </BackgroundGradient>
  )
}
