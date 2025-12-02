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


export function WelcomeCard() {
  return (
      <BackgroundGradient>
    <Card className="bg-primary border-none ">
      <CardHeader>
        <CardTitle className="text-background text-3xl my-4">
          Welcome to Genebuild Metadata
        </CardTitle>
      </CardHeader>

      <CardContent className="text-background text-sm">
        Quickly access and explore genome assembly and annotation metadata. Filter assemblies and annotations by project, release date, taxonomy, and more. Download ready-to-use tables and figures for your reports. Always up-to-date.
      </CardContent>

      <CardFooter className="justify-end mt-6 mb-4">
        <Button variant="secondary">
          <Link href="/assemblies">
            Search Assemblies <span className="ml-1">→</span>
          </Link>
        </Button>
      </CardFooter>
    </Card>
          </BackgroundGradient>
  )
}