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

export function WelcomeCard() {
  return (
    <Card className="bg-primary text-primary-foreground">
      <CardHeader>
        <CardTitle className="text-primary-foreground text-lg">
          Welcome to Genebuild Metadata
        </CardTitle>
      </CardHeader>

      <CardContent className="text-sidebar-ring text-sm">
        Quickly access and explore genome assembly and annotation metadata. Filter assemblies and annotations by project, release date, taxonomy, and more. Download ready-to-use tables and figures for your reports. Always up-to-date.
      </CardContent>

      <CardFooter className="justify-end">
        <Button variant="secondary">
          <Link href="/assemblies">
            Search Assemblies <span className="ml-1">→</span>
          </Link>
        </Button>
      </CardFooter>
    </Card>
  )
}