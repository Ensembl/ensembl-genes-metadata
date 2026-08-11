"use client";

import { useState } from "react";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
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
    gb_status: string;
    genebuilder: string;
    annotation_source: string;
    assembly_busco: number;
    protein_busco: number;
    ftp: string;
  };
  status: {
    title: string;
    description: string;
    action: string;
  };
  reasons: Reason[];
  timeline: TimelineEvent[];
};

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
    <div className="mx-auto mt-10 max-w-6xl space-y-8 px-4">
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
        <div
          key={res.summary.gca}
          className="space-y-8 rounded-xl border p-6"
        >
          <div>
            <h2 className="text-2xl font-bold">{res.summary.gca}</h2>

            <p className="text-lg font-medium">
              {res.status.title}
            </p>

            <p>{res.status.description}</p>

            <p className="text-muted-foreground">
              {res.status.action}
            </p>
          </div>

          {/* Summary */}

          <table className="w-full border-collapse text-sm">
            <tbody>
              <tr className="border-b">
                <td className="py-2 font-medium">Genebuilder</td>
                <td>{res.summary.genebuilder}</td>
              </tr>

              <tr className="border-b">
                <td className="py-2 font-medium">Annotation source</td>
                <td>{res.summary.annotation_source}</td>
              </tr>

              <tr className="border-b">
                <td className="py-2 font-medium">Genome BUSCO</td>
                <td>{res.summary.assembly_busco}%</td>
              </tr>

              <tr className="border-b">
                <td className="py-2 font-medium">Protein BUSCO</td>
                <td>{res.summary.protein_busco}%</td>
              </tr>

              <tr>
                <td className="py-2 font-medium">FTP</td>
                <td>
                  <a
                    href={res.summary.ftp}
                    className="text-blue-600 underline"
                  >
                    Download annotation
                  </a>
                </td>
              </tr>
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

          <Timeline
            orientation="horizontal"
            defaultValue={res.timeline.length}
            className="w-full"
          >
            {res.timeline.map((item, index) => (
              <TimelineItem
                key={`${item.date}-${index}`}
                step={index + 1}
              >
                <TimelineHeader>
                  <TimelineSeparator />
                  <TimelineDate>{item.date}</TimelineDate>
                  <TimelineTitle>{item.event}</TimelineTitle>
                  <TimelineIndicator />
                </TimelineHeader>

                <TimelineContent />
              </TimelineItem>
            ))}
          </Timeline>
        </div>
      ))}
    </div>
  );
}