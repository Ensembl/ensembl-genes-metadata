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
import { PROJECTS } from "@/features/projects/project-config";

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
          {PROJECTS.map((project) => (
            <Link
              key={project.slug}
              href={`/projects/${project.slug}`}
              className="group"
            >
              <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full dark:bg-secondary">
                <CardHeader>
                  <CardTitle>{project.title}</CardTitle>
                  <CardDescription>{project.description}</CardDescription>
                </CardHeader>
                <CardContent className="absolute bottom-4 right-4">
                  <ArrowRight className="w-5 h-5 text-muted-foreground group-hover:text-primary transition-colors" />
                </CardContent>
              </Card>
            </Link>
          ))}
        </div>
      </div>
    </div>
  );
}
