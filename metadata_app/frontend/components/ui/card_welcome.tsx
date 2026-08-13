"use client"
import * as React from "react"
import Link from "next/link"
import { Button } from "@/components/ui/button"


export function WelcomeCard() {
  return (
    <div className="w-full mx-auto max-w-2xl px-4 py-14 sm:px-6 sm:py-20 lg:max-w-7xl lg:px-8">
      <div className="container">
        <div className="mx-auto max-w-2xl text-center">
          <div className="hidden sm:mb-4 sm:flex sm:justify-center">
            <div className="relative rounded-full border px-3 py-1 text-sm">
              Track annotations by bioprojects{" "}
              <Link href="/projects" className="font-semibold text-chart-5">
                <span aria-hidden="true" className="absolute inset-0" />
                Read more
              </Link>
            </div>
          </div>
          <div className="text-center">
            <h1 className="text-foreground text-4xl leading-tight font-bold md:text-5xl">
              Welcome to Genebuild Metadata
            </h1>
            <p className="text-muted-foreground mt-4 text-lg leading-relaxed">
              Quickly access and explore genome assembly and annotation metadata.
                Filter assemblies and annotations by project, release date, taxonomy, and more.
                Download ready-to-use tables and figures for your reports. Always up-to-date.
            </p>
            <div className="mt-10 flex items-center justify-center gap-x-3">
              <Button size="lg" asChild>
                <Link href="/find-gca">Check genome annotation status</Link>
              </Button>
              <Button
                size="lg"
                variant="ghost"
                className="font-semibold"
                asChild
              >
              </Button>
            </div>
          </div>
        </div>
      </div>
    </div>
  )
}
