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

export default function ReportSelectorPage() {
  const title = "Generate reports";
  const description =
    "This page lets you generate reports on biodiversity projects and more. Select your filtering criteria to create visual summaries and download tables and figures. Select from the functions bellow.";



  return (
    <div className="flex items-center justify-center mt-15">
      <div className="container m-16 max-w-6xl">
        <h1 className="scroll-m-20 text-4xl font-extrabold tracking-tight text-balance">{title}</h1>
        <p className="leading-7 [&:not(:first-child)]:mt-6">{description}</p>
        <div className="grid grid-cols-2 gap-4 justify-center mt-8">
          <Link href="/report/asm" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>Assemblies</CardTitle>
                <CardDescription>
                  Create a report on available assemblies and find out what needs to be annotated
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

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
        </div>
      </div>
    </div>
  );
}