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
  return (
    <div className="flex items-center justify-center mt-60">
      <div className="container m-16 mt-10 max-w-6xl">
        <div className="grid grid-cols-2 gap-4 justify-center">
          <Link href="/report/asm" className="group">
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full">
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
            <Card className="relative hover:shadow-lg transition-shadow cursor-pointer h-full">
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