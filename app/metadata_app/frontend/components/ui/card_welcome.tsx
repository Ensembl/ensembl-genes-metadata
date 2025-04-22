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
        <CardTitle className="text-primary-foreground">
          Welcome to Genebuild Metadata
        </CardTitle>
      </CardHeader>

      <CardContent className="text-muted-foreground text-sm">
        Enim magna commodo minim et ut minim nisi aliquip ex ex. Ullamco sunt officia amet veniam.
        Labore pariatur commodo ut anim consectetur velit sunt irure esse pariatur deserunt.
      </CardContent>

      <CardFooter className="justify-end">
        <Button asChild className="bg-primary-foreground text-primary shadow hover:bg-primary-foreground/90 font-bold">
          <Link href="/assemblies">
            Search Assemblies <span className="ml-1">→</span>
          </Link>
        </Button>
      </CardFooter>
    </Card>
  )
}