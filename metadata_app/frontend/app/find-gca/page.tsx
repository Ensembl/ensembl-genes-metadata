"use client";

import { useState } from "react";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Field } from "@/components/ui/field";
import {
  Timeline,
  TimelineContent,
  TimelineDate,
  TimelineHeader,
  TimelineIndicator,
  TimelineItem,
  TimelineSeparator,
  TimelineTitle,
} from "@/components/ui/blocks/timeline";
import { cleanPayload } from "@/features/shared/filter-utils";

type TimelineEvent = {
  date: string;
  event: string;
};

type Reason = {
  title: string;
  description: string;
  action: string;
};

type Result = {
  summary: {
    gca: string;
    gb_status?: string | null;
    genebuilder?: string | null;
    annotation_source?: string | null;
    annotation_method?: string | null;
    assembly_busco?: number | null;
    protein_busco?: number | null;
    ftp?: string | null;
    asm_level?: string | null;
    contig_n50?: number | null;
    short_read_paired_end_illumina_lowest?: number | null;
    short_read_paired_end_illumina?: number | null;
  };
  status: {
    title: string;
    description: string;
    action: string;
  };
  reasons: Reason[];
  timeline: TimelineEvent[];
};

const TIMELINE_STEP_COUNT = 4;

export default function ReportSelectorPage() {
  const [gca, setGca] = useState("");
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState<Result[]>([]);
  const [error, setError] = useState("");

  const handleSearch = async () => {
    setLoading(true);
    setError("");

    try {
      const gcaArray = gca
        .split(",")
        .map((x) => x.trim())
        .filter(Boolean);

      const payload = cleanPayload({
        gca: gcaArray,
      });

      const response = await fetch("/api/gca_lookup/gca_lookup", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Accept: "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!response.ok) {
        throw new Error();
      }

      const data = await response.json();

      setResults(Array.isArray(data) ? data : [data]);
    } catch {
      setError("Unable to retrieve annotation status.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="mx-auto mt-10 max-w-6xl space-y-8 px-4 pb-16 sm:pb-12">
      <div>
        <h1 className="text-4xl font-bold">
          Find genome annotation status
        </h1>

        <p className="text-muted-foreground">
          Search one or more assembly accessions separated by comma and space.
        </p>
      </div>

      <Field orientation="horizontal" className="max-w-3xl">
        <Input
          value={gca}
          onChange={(e) => setGca(e.target.value)}
          placeholder="GCA_123456789.1, GCA_987654321.1"
          className="flex-1"
        />

        <Button disabled={loading} onClick={handleSearch}>
          {loading ? "Searching..." : "Search"}
        </Button>
      </Field>

      {error && (
        <p className="text-destructive">{error}</p>
      )}

      {results.map((res) => (
        <Card
          key={res.summary.gca}
          className="gap-6"
        >
          <CardHeader className="border-b">
            <CardTitle className="text-2xl font-bold">
              {res.summary.gca}
            </CardTitle>

            <div className="space-y-1">
              <p className="text-lg font-medium">
                {res.status.title}
              </p>

              <CardDescription className="text-base text-foreground">
                {res.status.description}
              </CardDescription>

              <CardDescription>
                {res.status.action}
              </CardDescription>
            </div>
          </CardHeader>

          <CardContent className="space-y-8">

          {/* Summary */}

          <table className="w-full border-collapse text-sm">
            <tbody>
              {res.summary.genebuilder && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Genebuilder</td>
                  <td>{res.summary.genebuilder}</td>
                </tr>
              )}

              {res.summary.annotation_source && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Annotation source</td>
                  <td>{res.summary.annotation_source}</td>
                </tr>
              )}

              {res.summary.annotation_method && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Annotation method</td>
                  <td>{res.summary.annotation_method}</td>
                </tr>
              )}

              {res.summary.asm_level && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Assembly level</td>
                  <td>{res.summary.asm_level}</td>
                </tr>
              )}

              {res.summary.contig_n50 !== undefined && res.summary.contig_n50 !== null && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Contig N50</td>
                  <td>{res.summary.contig_n50}</td>
                </tr>
              )}

              {res.summary.short_read_paired_end_illumina_lowest !== undefined &&
                res.summary.short_read_paired_end_illumina_lowest !== null && (
                  <tr className="border-b">
                    <td className="py-2 font-medium">RNA-seq lowest taxon ID</td>
                    <td>{res.summary.short_read_paired_end_illumina_lowest}</td>
                  </tr>
                )}

              {res.summary.short_read_paired_end_illumina !== undefined &&
                res.summary.short_read_paired_end_illumina !== null && (
                  <tr className="border-b">
                    <td className="py-2 font-medium">RNA-seq genus taxon ID</td>
                    <td>{res.summary.short_read_paired_end_illumina}</td>
                  </tr>
                )}

              {res.summary.assembly_busco !== undefined && res.summary.assembly_busco !== null && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Genome BUSCO</td>
                  <td>{res.summary.assembly_busco}</td>
                </tr>
              )}

              {res.summary.protein_busco !== undefined && res.summary.protein_busco !== null && (
                <tr className="border-b">
                  <td className="py-2 font-medium">Protein BUSCO</td>
                  <td>{res.summary.protein_busco}</td>
                </tr>
              )}

              {res.summary.ftp?.trim() && (
                <tr>
                  <td className="py-2 font-medium">FTP</td>
                  <td>
                    <a
                      href={res.summary.ftp}
                      target="_blank"
                      rel="noopener noreferrer"
                      className="text-color-chart-5 underline"
                    >
                      Download annotation
                    </a>
                  </td>
                </tr>
              )}
            </tbody>
          </table>

          {/* Reasons */}

          {res.reasons.length > 0 && (
            <div className="space-y-3">
              <h3 className="text-lg font-semibold">
                Possible reasons
              </h3>

              {res.reasons.map((reason) => (
                <div
                  key={reason.title}
                  className="rounded-lg border bg-muted p-4"
                >
                  <p className="font-medium">{reason.title}</p>
                  <p>{reason.description}</p>
                  <p className="text-muted-foreground">
                    {reason.action}
                  </p>
                </div>
              ))}
            </div>
          )}

          {/* Timeline */}

          {res.timeline.length > 0 &&
            (() => {
              const completedSteps = res.timeline.length;

              return (
                <Timeline
                  orientation="horizontal"
                  defaultValue={completedSteps}
                  className="w-full"
                >
                  {Array.from({ length: TIMELINE_STEP_COUNT }).map((_, index) => {
                    const item = res.timeline[index];
                    const isCompleted = Boolean(item);

                    return (
                      <TimelineItem
                        key={`timeline-step-${index + 1}`}
                        step={index + 1}
                      >
                        <TimelineHeader>
                          <TimelineSeparator />
                          <TimelineDate>{item?.date ?? ""}</TimelineDate>
                          <TimelineTitle>
                            {item?.event ?? ""}
                          </TimelineTitle>
                          <TimelineIndicator
                            className={
                              isCompleted ? "bg-primary" : "bg-background"
                            }
                          />
                        </TimelineHeader>

                        <TimelineContent />
                      </TimelineItem>
                    );
                  })}
                </Timeline>
              );
            })()}
          </CardContent>
        </Card>
      ))}
    </div>
  );
}
