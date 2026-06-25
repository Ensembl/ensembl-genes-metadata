"use client";

import React from "react";
import Link from "next/link";
import {
  Card,
  CardHeader,
  CardTitle,
  CardDescription,
  CardContent,
} from "@/components/ui/card";
import { ArrowRight } from "lucide-react";
import { BackgroundGradient } from "@/components/ui/backround-gradient"


import {CardsReportTable} from "@/components/ui/report_card_table"
import {ProjectsBar} from "@/components/ui/report_per_year"

export default function ReportSelectorPage() {
  const title = "Generate reports";
  const description =
    "This page lets you generate reports on biodiversity projects and more. Select your filtering criteria to create visual summaries and download tables and figures. Select from the functions bellow.";



  return (
    <div className="mt-10 flex w-full justify-center overflow-x-hidden px-4 sm:px-6 lg:px-8">
      <div className="grid w-full min-w-0 max-w-6xl gap-8">
        <h1 className="scroll-m-20 break-words text-3xl font-extrabold tracking-tight text-balance sm:text-4xl">{title}</h1>
        <p className="max-w-full break-words leading-7 [&:not(:first-child)]:mt-6">{description}</p>
        <div className="grid w-full min-w-0 auto-rows-fr grid-cols-1 items-stretch gap-4 md:grid-cols-3">
         <BackgroundGradient containerClassName="h-full min-w-0" className="flex h-full min-w-0">
            <Link href="/report/asm" className="group flex h-full w-full min-w-0">
            <Card className="relative h-full w-full min-w-0 cursor-pointer hover:shadow-lg transition-shadow dark:bg-secondary">
                <CardHeader>
                <CardTitle>Assemblies</CardTitle>
                <CardDescription>
                  Create a report on available assemblies
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>
             </BackgroundGradient>
        <BackgroundGradient containerClassName="h-full min-w-0" className="flex h-full min-w-0">
          <Link href="/report/anno" className="group flex h-full w-full min-w-0">
            <Card className="relative h-full w-full min-w-0 cursor-pointer hover:shadow-lg transition-shadow dark:bg-secondary">
              <CardHeader>
                <CardTitle>Annotations</CardTitle>
                <CardDescription>
                  Create a report on available annotations by Genebuild
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>
        </BackgroundGradient>
        <BackgroundGradient containerClassName="h-full min-w-0" className="flex h-full min-w-0">
          <Link href="/report/annotation-qc" className="group flex h-full w-full min-w-0">
            <Card className="relative h-full w-full min-w-0 cursor-pointer hover:shadow-lg transition-shadow dark:bg-secondary">
              <CardHeader>
                <CardTitle>Annotation QC</CardTitle>
                <CardDescription>
                  Review assembly and annotation quality metrics
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>
        </BackgroundGradient>
        </div>

        <h2 className="scroll-m-20 text-2xl font-semibold tracking-tight">
          Quick overview
        </h2>

        <div className="min-w-0">
          <ProjectsBar></ProjectsBar>
        </div>

        <div className="mb-8 min-w-0">
          <CardsReportTable></CardsReportTable>
        </div>


      </div>
    </div>
  );
}
