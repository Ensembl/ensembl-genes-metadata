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
  const title = "Biodiversity Projects";
  const description =
    "Select Biodiversity Project to track annotation status. High priority projects are marked with *. Shows primary, chromosome or complete genome level assemblies per project.";



  return (
    <div className="flex items-center justify-center mt-15">
      <div className="container m-16 max-w-6xl">
        <h1 className="scroll-m-20 text-4xl font-extrabold tracking-tight text-balance">{title}</h1>
        <p className="leading-7 [&:not(:first-child)]:mt-6">{description}</p>
        <div className="grid grid-cols-2 gap-4 justify-center mt-8">
          <Link href="/projects/erga-bge" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>ERGA-BGE*</CardTitle>
                <CardDescription>
                  European Reference Genome Atlas Biodiversity Genomics Europe Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

          <Link href="/projects/erga" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>ERGA*</CardTitle>
                <CardDescription>
                  European Reference Genome Atlas Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

          <Link href="/projects/erga-pilot" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>ERGA-pilot*</CardTitle>
                <CardDescription>
                  European Reference Genome Atlas Pilot Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

          <Link href="/projects/cbp" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>CBP*</CardTitle>
                <CardDescription>
                  Canadian BioGenome Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

          <Link href="/projects/ebp" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>EBP*</CardTitle>
                <CardDescription>
                  Earth BioGenome Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

          <Link href="/projects/vgp" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>VGP*</CardTitle>
                <CardDescription>
                  Vertebrate Genomes Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>


          <Link href="/projects/dtol" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>DTol</CardTitle>
                <CardDescription>
                  Darwin Tree of Life Project
                </CardDescription>
              </CardHeader>
              <CardContent className="absolute bottom-4 right-4">
                <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
              </CardContent>
            </Card>
          </Link>

          <Link href="/projects/asg" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
              <CardHeader>
                <CardTitle>ASG</CardTitle>
                <CardDescription>
                  Aquatic Symbiosis Genomics Project
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