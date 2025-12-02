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
    <div className="flex items-center justify-center mt-15">
      <div className="grid max-w-6xl gap-8">
        <h1 className="scroll-m-20 text-4xl font-extrabold tracking-tight text-balance">{title}</h1>
        <p className="leading-7 [&:not(:first-child)]:mt-6">{description}</p>
        <div className="grid grid-cols-2 gap-4 justify-center">
         <BackgroundGradient>
            <Link href="/report/asm" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
                <CardHeader>
                <CardTitle>Non-annotated assemblies</CardTitle>
                <CardDescription>
                  Create a report on available non-annotated assemblies
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>
             </BackgroundGradient>
        <BackgroundGradient>
          <Link href="/report/anno" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
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
        </div>

        <h2 className="scroll-m-20 text-2xl font-semibold tracking-tight">
          Quick overview
        </h2>

        <div>
          <ProjectsBar></ProjectsBar>
        </div>

        <div className="mb-8">
          <CardsReportTable></CardsReportTable>
        </div>


      </div>
    </div>
  );
}