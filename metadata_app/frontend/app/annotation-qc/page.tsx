"use client";

import React, { useMemo, useState } from "react";
import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  Line,
  LineChart,
  Scatter,
  ScatterChart,
  XAxis,
  YAxis,
} from "recharts";
import { Activity, Search, ShieldAlert } from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import {
  ChartConfig,
  ChartContainer,
  ChartTooltip,
  ChartTooltipContent,
} from "@/components/ui/chart";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";

type GenebuildStatus = "live" | "in progress" | "handover" | "blocked";

type QcRecord = {
  gca: string;
  species: string;
  taxonomy: string;
  taxonId: string;
  status: GenebuildStatus;
  contigN50Mb: number;
  scaffoldN50Mb: number;
  genomeSizeGb: number;
  gcPercent: number;
  geneCount: number;
  proteinCoding: number;
  transcriptCount: number;
  genomeBuscoComplete: number;
  proteinBuscoComplete: number;
  buscoComplete: number;
  buscoDuplicated: number;
  buscoFragmented: number;
  buscoMissing: number;
  repeatPercent: number;
  annotationScore: number;
  assemblyScore: number;
  lncRnaCount: number;
  pseudogeneCount: number;
  smallRnaCount: number;
  pcoa1: number;
  pcoa2: number;
  outlierScore: number;
};

type BuscoComparison = {
  name: string;
  genomeComplete: number;
  proteinComplete: number;
  complete: number;
  duplicated: number;
  fragmented: number;
  missing: number;
  difference: number;
  differenceFlagged: boolean;
  duplicationFlagged: boolean;
};

const QC_RECORDS: QcRecord[] = [
  {
    gca: "GCA_000001405.29",
    species: "Homo sapiens",
    taxonomy: "Primates",
    taxonId: "9606",
    status: "live",
    contigN50Mb: 57.9,
    scaffoldN50Mb: 67.8,
    genomeSizeGb: 3.1,
    gcPercent: 41.0,
    geneCount: 21421,
    proteinCoding: 19957,
    transcriptCount: 252316,
    genomeBuscoComplete: 98.5,
    proteinBuscoComplete: 98.1,
    buscoComplete: 98.1,
    buscoDuplicated: 1.2,
    buscoFragmented: 0.8,
    buscoMissing: 1.1,
    repeatPercent: 50.6,
    annotationScore: 96,
    assemblyScore: 94,
    lncRnaCount: 18422,
    pseudogeneCount: 14391,
    smallRnaCount: 7568,
    pcoa1: -0.18,
    pcoa2: 0.09,
    outlierScore: 0.6,
  },
  {
    gca: "GCA_000001635.9",
    species: "Mus musculus",
    taxonomy: "Rodentia",
    taxonId: "10090",
    status: "live",
    contigN50Mb: 32.4,
    scaffoldN50Mb: 61.5,
    genomeSizeGb: 2.73,
    gcPercent: 42.0,
    geneCount: 22612,
    proteinCoding: 21704,
    transcriptCount: 149882,
    genomeBuscoComplete: 97.9,
    proteinBuscoComplete: 97.4,
    buscoComplete: 97.4,
    buscoDuplicated: 1.6,
    buscoFragmented: 1.2,
    buscoMissing: 1.4,
    repeatPercent: 45.3,
    annotationScore: 93,
    assemblyScore: 91,
    lncRnaCount: 14180,
    pseudogeneCount: 10862,
    smallRnaCount: 6120,
    pcoa1: -0.11,
    pcoa2: 0.05,
    outlierScore: 0.8,
  },
  {
    gca: "GCA_000002035.6",
    species: "Danio rerio",
    taxonomy: "Actinopteri",
    taxonId: "7955",
    status: "handover",
    contigN50Mb: 1.7,
    scaffoldN50Mb: 55.8,
    genomeSizeGb: 1.37,
    gcPercent: 36.5,
    geneCount: 25143,
    proteinCoding: 24638,
    transcriptCount: 60321,
    genomeBuscoComplete: 94.0,
    proteinBuscoComplete: 96.2,
    buscoComplete: 96.2,
    buscoDuplicated: 3.8,
    buscoFragmented: 2.1,
    buscoMissing: 1.7,
    repeatPercent: 52.8,
    annotationScore: 89,
    assemblyScore: 84,
    lncRnaCount: 9480,
    pseudogeneCount: 5042,
    smallRnaCount: 2934,
    pcoa1: 0.2,
    pcoa2: -0.08,
    outlierScore: 1.3,
  },
  {
    gca: "GCA_000146045.2",
    species: "Arabidopsis thaliana",
    taxonomy: "Brassicaceae",
    taxonId: "3702",
    status: "live",
    contigN50Mb: 12.1,
    scaffoldN50Mb: 23.5,
    genomeSizeGb: 0.12,
    gcPercent: 36.1,
    geneCount: 27655,
    proteinCoding: 27445,
    transcriptCount: 48219,
    genomeBuscoComplete: 98.9,
    proteinBuscoComplete: 98.7,
    buscoComplete: 98.7,
    buscoDuplicated: 2.1,
    buscoFragmented: 0.7,
    buscoMissing: 0.6,
    repeatPercent: 15.4,
    annotationScore: 95,
    assemblyScore: 92,
    lncRnaCount: 4120,
    pseudogeneCount: 924,
    smallRnaCount: 831,
    pcoa1: -0.35,
    pcoa2: -0.24,
    outlierScore: 0.9,
  },
  {
    gca: "GCA_905319855.2",
    species: "Salmo salar",
    taxonomy: "Actinopteri",
    taxonId: "8030",
    status: "in progress",
    contigN50Mb: 4.2,
    scaffoldN50Mb: 35.6,
    genomeSizeGb: 2.97,
    gcPercent: 43.4,
    geneCount: 46722,
    proteinCoding: 44081,
    transcriptCount: 81792,
    genomeBuscoComplete: 90.1,
    proteinBuscoComplete: 94.5,
    buscoComplete: 94.5,
    buscoDuplicated: 28.4,
    buscoFragmented: 2.8,
    buscoMissing: 2.7,
    repeatPercent: 59.1,
    annotationScore: 82,
    assemblyScore: 80,
    lncRnaCount: 20410,
    pseudogeneCount: 12254,
    smallRnaCount: 6881,
    pcoa1: 0.44,
    pcoa2: -0.18,
    outlierScore: 2.7,
  },
  {
    gca: "GCA_964030725.1",
    species: "Bombus terrestris",
    taxonomy: "Arthropoda",
    taxonId: "30195",
    status: "blocked",
    contigN50Mb: 0.38,
    scaffoldN50Mb: 2.4,
    genomeSizeGb: 0.28,
    gcPercent: 35.2,
    geneCount: 14102,
    proteinCoding: 13218,
    transcriptCount: 18502,
    genomeBuscoComplete: 78.2,
    proteinBuscoComplete: 83.9,
    buscoComplete: 83.9,
    buscoDuplicated: 0.9,
    buscoFragmented: 7.4,
    buscoMissing: 8.7,
    repeatPercent: 18.2,
    annotationScore: 61,
    assemblyScore: 58,
    lncRnaCount: 3288,
    pseudogeneCount: 1844,
    smallRnaCount: 912,
    pcoa1: 0.72,
    pcoa2: 0.43,
    outlierScore: 3.8,
  },
  {
    gca: "GCA_963924645.1",
    species: "Canis lupus familiaris",
    taxonomy: "Carnivora",
    taxonId: "9615",
    status: "handover",
    contigN50Mb: 19.6,
    scaffoldN50Mb: 41.1,
    genomeSizeGb: 2.42,
    gcPercent: 41.6,
    geneCount: 20954,
    proteinCoding: 19885,
    transcriptCount: 73210,
    genomeBuscoComplete: 95.9,
    proteinBuscoComplete: 96.8,
    buscoComplete: 96.8,
    buscoDuplicated: 1.1,
    buscoFragmented: 1.5,
    buscoMissing: 1.7,
    repeatPercent: 41.7,
    annotationScore: 88,
    assemblyScore: 87,
    lncRnaCount: 8262,
    pseudogeneCount: 5940,
    smallRnaCount: 2446,
    pcoa1: -0.08,
    pcoa2: 0.18,
    outlierScore: 1.1,
  },
  {
    gca: "GCA_949987765.1",
    species: "Lynx lynx",
    taxonomy: "Carnivora",
    taxonId: "13125",
    status: "in progress",
    contigN50Mb: 8.3,
    scaffoldN50Mb: 31.4,
    genomeSizeGb: 2.51,
    gcPercent: 41.2,
    geneCount: 20318,
    proteinCoding: 19304,
    transcriptCount: 45622,
    genomeBuscoComplete: 89.4,
    proteinBuscoComplete: 93.8,
    buscoComplete: 93.8,
    buscoDuplicated: 1.0,
    buscoFragmented: 3.2,
    buscoMissing: 3.0,
    repeatPercent: 39.8,
    annotationScore: 81,
    assemblyScore: 79,
    lncRnaCount: 7140,
    pseudogeneCount: 4812,
    smallRnaCount: 2180,
    pcoa1: 0.08,
    pcoa2: 0.26,
    outlierScore: 1.9,
  },
];

const statusOptions: Array<"all" | GenebuildStatus> = [
  "all",
  "live",
  "in progress",
  "handover",
  "blocked",
];

const qcChartConfig = {
  genomeComplete: {
    label: "Genome complete",
    color: "var(--chart-2)",
  },
  proteinComplete: {
    label: "Protein complete",
    color: "var(--chart-4)",
  },
  mean: {
    label: "Mean",
    color: "var(--chart-2)",
  },
  median: {
    label: "Median",
    color: "var(--chart-4)",
  },
  assembly: {
    label: "Assembly",
    color: "var(--chart-2)",
  },
  annotation: {
    label: "Annotation",
    color: "var(--chart-4)",
  },
  pcoa1: {
    label: "PCoA 1",
  },
  pcoa2: {
    label: "PCoA 2",
  },
} satisfies ChartConfig;

const mean = (values: number[]) =>
  values.length ? values.reduce((total, value) => total + value, 0) / values.length : 0;

const median = (values: number[]) => {
  if (!values.length) return 0;
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 0
    ? (sorted[middle - 1] + sorted[middle]) / 2
    : sorted[middle];
};

const formatNumber = (value: number, digits = 1) =>
  new Intl.NumberFormat("en-GB", {
    maximumFractionDigits: digits,
    minimumFractionDigits: digits,
  }).format(value);

const getStatusBadge = (status: GenebuildStatus) => {
  if (status === "blocked") return "destructive";
  if (status === "live") return "default";
  return "secondary";
};

const getBuscoDifference = (record: QcRecord) =>
  Math.abs(record.genomeBuscoComplete - record.proteinBuscoComplete);

const toBuscoComparison = (
  record: QcRecord,
  averageDifference: number,
  averageDuplicated: number,
  name = record.species,
): BuscoComparison => {
  const difference = getBuscoDifference(record);

  return {
    name,
    genomeComplete: record.genomeBuscoComplete,
    proteinComplete: record.proteinBuscoComplete,
    complete: record.buscoComplete,
    duplicated: record.buscoDuplicated,
    fragmented: record.buscoFragmented,
    missing: record.buscoMissing,
    difference,
    differenceFlagged: difference > averageDifference,
    duplicationFlagged: record.buscoDuplicated > averageDuplicated,
  };
};

const getBuscoMedian = (
  records: QcRecord[],
  averageDifference: number,
  averageDuplicated: number,
  name: string,
): BuscoComparison => {
  const genomeComplete = median(records.map((record) => record.genomeBuscoComplete));
  const proteinComplete = median(records.map((record) => record.proteinBuscoComplete));
  const difference = Math.abs(genomeComplete - proteinComplete);

  return {
    name,
    genomeComplete,
    proteinComplete,
    complete: median(records.map((record) => record.buscoComplete)),
    duplicated: median(records.map((record) => record.buscoDuplicated)),
    fragmented: median(records.map((record) => record.buscoFragmented)),
    missing: median(records.map((record) => record.buscoMissing)),
    difference,
    differenceFlagged: difference > averageDifference,
    duplicationFlagged: median(records.map((record) => record.buscoDuplicated)) > averageDuplicated,
  };
};

export default function AnnotationQcPage() {
  const [gca, setGca] = useState("");
  const [taxonomy, setTaxonomy] = useState("Carnivora");
  const [status, setStatus] = useState<(typeof statusOptions)[number]>("all");
  const [submitted, setSubmitted] = useState({ gca: "", taxonomy: "Carnivora", status: "all" });

  const filteredRecords = useMemo(() => {
    const statusMatches = (record: QcRecord) =>
      submitted.status === "all" || record.status === submitted.status;

    if (submitted.gca.trim()) {
      const query = submitted.gca.trim().toLowerCase();
      return QC_RECORDS.filter(
        (record) => record.gca.toLowerCase() === query && statusMatches(record),
      );
    }

    const taxonQuery = submitted.taxonomy.trim().toLowerCase();
    return QC_RECORDS.filter(
      (record) =>
        statusMatches(record) &&
        (!taxonQuery ||
          record.taxonomy.toLowerCase().includes(taxonQuery) ||
          record.species.toLowerCase().includes(taxonQuery) ||
          record.taxonId === taxonQuery),
    );
  }, [submitted]);

  const singleRecord = submitted.gca.trim() ? filteredRecords[0] : null;

  const cohortStats = useMemo(() => {
    const metrics = [
      { key: "proteinCoding", label: "Coding genes" },
      { key: "lncRnaCount", label: "lncRNA genes" },
      { key: "pseudogeneCount", label: "Pseudogenes" },
      { key: "smallRnaCount", label: "Small RNA genes" },
      { key: "transcriptCount", label: "Transcripts" },
    ] as const;

    return metrics.map((metric) => {
      const values = filteredRecords.map((record) => record[metric.key]);
      const metricMean = mean(values);
      const metricMedian = median(values);

      return {
        metric: metric.label,
        mean: metricMean,
        median: metricMedian,
      };
    });
  }, [filteredRecords]);

  const outliers = useMemo(
    () => filteredRecords.filter((record) => record.outlierScore >= 2).sort((a, b) => b.outlierScore - a.outlierScore),
    [filteredRecords],
  );

  const scoreTrend = filteredRecords.map((record) => ({
    species: record.species.split(" ")[0],
    assembly: record.assemblyScore,
    annotation: record.annotationScore,
  }));

  const pcoaSummary = useMemo(() => {
    const pcoa1Values = filteredRecords.map((record) => record.pcoa1);
    const pcoa2Values = filteredRecords.map((record) => record.pcoa2);
    const formatRange = (values: number[]) =>
      values.length
        ? `${formatNumber(Math.min(...values), 2)} to ${formatNumber(Math.max(...values), 2)}`
        : "n/a";

    return {
      assemblies: filteredRecords.length,
      outliers: outliers.length,
      pcoa1Range: formatRange(pcoa1Values),
      pcoa2Range: formatRange(pcoa2Values),
    };
  }, [filteredRecords, outliers.length]);

  const buscoAverageDifference = useMemo(() => {
    const records = singleRecord
      ? QC_RECORDS.filter((record) => record.taxonomy === singleRecord.taxonomy)
      : filteredRecords;

    return mean(records.map(getBuscoDifference));
  }, [filteredRecords, singleRecord]);

  const buscoAverageDuplicated = useMemo(() => {
    const records = singleRecord
      ? QC_RECORDS.filter((record) => record.taxonomy === singleRecord.taxonomy)
      : filteredRecords;

    return mean(records.map((record) => record.buscoDuplicated));
  }, [filteredRecords, singleRecord]);

  const buscoComparison = useMemo(() => {
    if (singleRecord) {
      const taxonomyPeers = QC_RECORDS.filter(
        (record) => record.taxonomy === singleRecord.taxonomy,
      );

      return [
        toBuscoComparison(
          singleRecord,
          buscoAverageDifference,
          buscoAverageDuplicated,
          singleRecord.species,
        ),
        getBuscoMedian(
          taxonomyPeers,
          buscoAverageDifference,
          buscoAverageDuplicated,
          `${singleRecord.taxonomy} median`,
        ),
      ];
    }

    return filteredRecords.map((record) =>
      toBuscoComparison(
        record,
        buscoAverageDifference,
        buscoAverageDuplicated,
        record.species.split(" ")[0],
      ),
    );
  }, [buscoAverageDifference, buscoAverageDuplicated, filteredRecords, singleRecord]);

  const submitFilters = () => {
    setSubmitted({ gca, taxonomy, status });
  };

  return (
    <div className="min-h-screen flex justify-center px-4 py-10 sm:px-6 lg:px-8">
      <div className="w-full max-w-7xl">
        <div className="rounded-2xl bg-secondary p-4 gap-10 shadow-lg sm:p-8">
          <div className="mb-8">
            <h1 className="text-2xl font-bold mb-2">Annotation QC</h1>
            <p className="text-muted-foreground">
              Review assembly and annotation quality metrics by GCA accession or by taxonomy cohort.
              This view is wired to mock data for frontend validation.
            </p>
          </div>

          <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
            <div>
              <Label htmlFor="gca-accession">GCA accession</Label>
              <Input
                id="gca-accession"
                className="mt-3 bg-filter-input-bg dark:bg-transparent"
                placeholder="GCA_000001405.29"
                value={gca}
                onChange={(event) => setGca(event.target.value)}
              />
            </div>
            <div>
              <Label htmlFor="taxonomy-filter">Taxonomy</Label>
              <Input
                id="taxonomy-filter"
                className="mt-3 bg-filter-input-bg dark:bg-transparent"
                placeholder="Carnivora, Primates, 9606"
                value={taxonomy}
                onChange={(event) => setTaxonomy(event.target.value)}
                disabled={gca.trim().length > 0}
              />
            </div>
            <div>
              <Label>Genebuild status</Label>
              <Select value={status} onValueChange={(value) => setStatus(value as typeof status)}>
                <SelectTrigger className="mt-3 w-full bg-filter-input-bg dark:bg-transparent">
                  <SelectValue placeholder="Select status" />
                </SelectTrigger>
                <SelectContent>
                  {statusOptions.map((option) => (
                    <SelectItem key={option} value={option}>
                      {option === "all" ? "All statuses" : option}
                    </SelectItem>
                  ))}
                </SelectContent>
              </Select>
            </div>
          </div>

          <div className="mt-6 flex justify-end">
            <Button className="w-full sm:w-auto" onClick={submitFilters}>
              <Search className="mr-2 h-4 w-4" />
              Run QC report
            </Button>
          </div>
        </div>

        {submitted.gca.trim() && !singleRecord && (
          <Card className="mt-8">
            <CardHeader>
              <CardTitle>No GCA match</CardTitle>
              <CardDescription>
                The mock data does not include that accession with the selected status filter.
              </CardDescription>
            </CardHeader>
          </Card>
        )}

        {singleRecord && (
          <div className="mt-8 grid grid-cols-1 gap-6">
            <div className="flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
              <div>
                <h2 className="text-xl font-semibold">{singleRecord.species}</h2>
                <p className="text-sm text-muted-foreground">{singleRecord.gca}</p>
              </div>
              <Badge variant={getStatusBadge(singleRecord.status)}>{singleRecord.status}</Badge>
            </div>

            <div className="grid grid-cols-1 gap-4 md:grid-cols-4">
              <MetricCard label="Assembly QC" value={singleRecord.assemblyScore} suffix="/100" />
              <MetricCard label="Annotation QC" value={singleRecord.annotationScore} suffix="/100" />
              <MetricCard label="Protein BUSCO" value={singleRecord.proteinBuscoComplete} suffix="%" />
              <MetricCard label="Outlier score" value={singleRecord.outlierScore} />
            </div>

            <BuscoComparisonBlock
              data={buscoComparison}
              averageDifference={buscoAverageDifference}
              averageDuplicated={buscoAverageDuplicated}
              title="BUSCO comparison"
              description={`Selected assembly compared with the ${singleRecord.taxonomy} median.`}
            />

            <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
              <Card>
                <CardHeader>
                  <CardTitle>Assembly metrics</CardTitle>
                  <CardDescription>N50, genome size, GC, and repeat content.</CardDescription>
                </CardHeader>
                <CardContent>
                  <MetricRows
                    rows={[
                      ["Contig N50", `${formatNumber(singleRecord.contigN50Mb)} Mb`],
                      ["Scaffold N50", `${formatNumber(singleRecord.scaffoldN50Mb)} Mb`],
                      ["Genome size", `${formatNumber(singleRecord.genomeSizeGb, 2)} Gb`],
                      ["GC", `${formatNumber(singleRecord.gcPercent)}%`],
                      ["Repeats", `${formatNumber(singleRecord.repeatPercent)}%`],
                    ]}
                  />
                </CardContent>
              </Card>

              <Card>
                <CardHeader>
                  <CardTitle>Annotation metrics</CardTitle>
                  <CardDescription>Gene model and transcript-level indicators.</CardDescription>
                </CardHeader>
                <CardContent>
                  <MetricRows
                    rows={[
                      ["Genes", singleRecord.geneCount.toLocaleString("en-GB")],
                      ["Protein coding", singleRecord.proteinCoding.toLocaleString("en-GB")],
                      ["Transcripts", singleRecord.transcriptCount.toLocaleString("en-GB")],
                      ["Genome BUSCO complete", `${formatNumber(singleRecord.genomeBuscoComplete)}%`],
                      ["Protein BUSCO complete", `${formatNumber(singleRecord.proteinBuscoComplete)}%`],
                      ["Genome/protein delta", `${formatNumber(getBuscoDifference(singleRecord))}%`],
                      ["BUSCO duplicated", `${formatNumber(singleRecord.buscoDuplicated)}%`],
                      ["Taxonomy", `${singleRecord.taxonomy} (${singleRecord.taxonId})`],
                    ]}
                  />
                </CardContent>
              </Card>
            </div>
          </div>
        )}

        {!submitted.gca.trim() && (
          <div className="mt-8 grid grid-cols-1 gap-6">
            <div className="flex flex-col gap-3 sm:flex-row sm:items-center sm:justify-between">
              <div>
                <h2 className="text-xl font-semibold">Taxonomy cohort analysis</h2>
                <p className="text-sm text-muted-foreground">
                  {filteredRecords.length} assemblies matched {submitted.taxonomy || "all taxonomy"}.
                </p>
              </div>
              <Badge variant="outline">Mock cohort</Badge>
            </div>

            <div className="grid grid-cols-1 gap-4 md:grid-cols-4">
              <MetricCard label="Assemblies" value={filteredRecords.length} />
              <MetricCard label="Mean assembly QC" value={mean(filteredRecords.map((item) => item.assemblyScore))} />
              <MetricCard label="Median annotation QC" value={median(filteredRecords.map((item) => item.annotationScore))} />
              <MetricCard label="Outliers" value={outliers.length} />
            </div>

            <BuscoComparisonBlock
              data={buscoComparison}
              averageDifference={buscoAverageDifference}
              averageDuplicated={buscoAverageDuplicated}
              title="BUSCO comparison"
              description="Genome and protein BUSCO complete scores, plus C/D/F/M protein BUSCO breakdown across the matched cohort."
            />

            <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
              <Card>
                <CardHeader>
                  <CardTitle>Gene model means and medians</CardTitle>
                  <CardDescription>
                    Coding, non-coding, and biotype count summaries across the matched cohort.
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <ChartContainer config={qcChartConfig} className="h-72 w-full">
                      <BarChart data={cohortStats}>
                        <CartesianGrid strokeDasharray="3 3" />
                        <XAxis dataKey="metric" tick={{ fontSize: 11 }} interval={0} angle={-18} textAnchor="end" height={70} />
                        <YAxis tickFormatter={(value) => Number(value).toLocaleString("en-GB")} />
                        <ChartTooltip content={<ChartTooltipContent />} />
                        <Bar dataKey="mean" fill="var(--chart-2)" name="Mean" radius={[4, 4, 0, 0]} />
                        <Bar dataKey="median" fill="var(--chart-4)" name="Median" radius={[4, 4, 0, 0]} />
                      </BarChart>
                  </ChartContainer>
                  <div className="mt-4 overflow-x-auto rounded-md border">
                    <Table>
                      <TableHeader>
                        <TableRow>
                          <TableHead>Metric</TableHead>
                          <TableHead className="text-right">Mean</TableHead>
                          <TableHead className="text-right">Median</TableHead>
                        </TableRow>
                      </TableHeader>
                      <TableBody>
                        {cohortStats.map((item) => (
                          <TableRow key={item.metric}>
                            <TableCell>{item.metric}</TableCell>
                            <TableCell className="text-right">
                              {Math.round(item.mean).toLocaleString("en-GB")}
                            </TableCell>
                            <TableCell className="text-right">
                              {Math.round(item.median).toLocaleString("en-GB")}
                            </TableCell>
                          </TableRow>
                        ))}
                      </TableBody>
                    </Table>
                  </div>
                </CardContent>
              </Card>

              <Card>
                <CardHeader>
                  <CardTitle>PCoA metric space</CardTitle>
                  <CardDescription>Mock ordination using assembly and annotation QC features.</CardDescription>
                </CardHeader>
                <CardContent>
                  <ChartContainer config={qcChartConfig} className="h-80 w-full">
                    <ScatterChart>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis dataKey="pcoa1" name="PCoA 1" type="number" />
                      <YAxis dataKey="pcoa2" name="PCoA 2" type="number" />
                      <ChartTooltip
                        cursor={{ strokeDasharray: "3 3" }}
                        content={<ChartTooltipContent />}
                      />
                      <Scatter data={filteredRecords} name="Assembly">
                        {filteredRecords.map((entry) => (
                          <Cell
                            key={entry.gca}
                            fill={entry.outlierScore >= 2 ? "var(--destructive)" : "var(--chart-2)"}
                          />
                        ))}
                      </Scatter>
                    </ScatterChart>
                  </ChartContainer>
                  <div className="mt-4 grid grid-cols-1 gap-3 text-sm sm:grid-cols-2">
                    <div className="rounded-md border p-3">
                      <div className="text-muted-foreground">Assemblies plotted</div>
                      <div className="mt-1 text-lg font-semibold">{pcoaSummary.assemblies}</div>
                    </div>
                    <div className="rounded-md border p-3">
                      <div className="text-muted-foreground">Outliers</div>
                      <div className="mt-1 text-lg font-semibold">{pcoaSummary.outliers}</div>
                    </div>
                    <div className="rounded-md border p-3">
                      <div className="text-muted-foreground">PCoA 1 range</div>
                      <div className="mt-1 font-medium">{pcoaSummary.pcoa1Range}</div>
                    </div>
                    <div className="rounded-md border p-3">
                      <div className="text-muted-foreground">PCoA 2 range</div>
                      <div className="mt-1 font-medium">{pcoaSummary.pcoa2Range}</div>
                    </div>
                  </div>
                </CardContent>
              </Card>

              <Card>
                <CardHeader>
                  <CardTitle>QC score profile</CardTitle>
                  <CardDescription>Assembly and annotation score by matched assembly.</CardDescription>
                </CardHeader>
                <CardContent>
                  <ChartContainer config={qcChartConfig} className="h-80 w-full">
                    <LineChart data={scoreTrend}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis dataKey="species" />
                      <YAxis domain={[50, 100]} />
                      <ChartTooltip content={<ChartTooltipContent />} />
                      <Line type="monotone" dataKey="assembly" stroke="var(--chart-2)" strokeWidth={2} />
                      <Line type="monotone" dataKey="annotation" stroke="var(--chart-4)" strokeWidth={2} />
                    </LineChart>
                  </ChartContainer>
                </CardContent>
              </Card>

              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center gap-2">
                    <ShieldAlert className="h-5 w-5" />
                    Outlier analysis
                  </CardTitle>
                  <CardDescription>Records with elevated multivariate distance.</CardDescription>
                </CardHeader>
                <CardContent>
                  {outliers.length === 0 ? (
                    <div className="flex h-52 items-center justify-center text-sm text-muted-foreground">
                      No outliers in the current mock cohort.
                    </div>
                  ) : (
                    <div className="overflow-x-auto">
                      <Table>
                        <TableHeader>
                          <TableRow>
                            <TableHead>GCA</TableHead>
                            <TableHead>Species</TableHead>
                            <TableHead>Status</TableHead>
                            <TableHead className="text-right">Score</TableHead>
                          </TableRow>
                        </TableHeader>
                        <TableBody>
                          {outliers.map((record) => (
                            <TableRow key={record.gca}>
                              <TableCell className="font-mono text-xs">{record.gca}</TableCell>
                              <TableCell>{record.species}</TableCell>
                              <TableCell>
                                <Badge variant={getStatusBadge(record.status)}>{record.status}</Badge>
                              </TableCell>
                              <TableCell className="text-right">{formatNumber(record.outlierScore)}</TableCell>
                            </TableRow>
                          ))}
                        </TableBody>
                      </Table>
                    </div>
                  )}
                </CardContent>
              </Card>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}

function MetricCard({
  label,
  value,
  suffix = "",
}: {
  label: string;
  value: number;
  suffix?: string;
}) {
  return (
    <Card>
      <CardHeader className="pb-2">
        <CardDescription>{label}</CardDescription>
        <CardTitle className="flex items-baseline gap-1 text-2xl">
          <Activity className="h-4 w-4 text-muted-foreground" />
          {Number.isInteger(value) ? value.toLocaleString("en-GB") : formatNumber(value)}
          <span className="text-sm font-normal text-muted-foreground">{suffix}</span>
        </CardTitle>
      </CardHeader>
    </Card>
  );
}

function BuscoComparisonBlock({
  data,
  averageDifference,
  averageDuplicated,
  title,
  description,
}: {
  data: BuscoComparison[];
  averageDifference: number;
  averageDuplicated: number;
  title: string;
  description: string;
}) {
  return (
    <Card>
      <CardHeader>
        <CardTitle>{title}</CardTitle>
        <CardDescription>
          {description} Flag threshold: above {formatNumber(averageDifference)}% average
          genome/protein difference or above {formatNumber(averageDuplicated)}% average BUSCO duplication.
        </CardDescription>
      </CardHeader>
      <CardContent>
        <div className="grid grid-cols-1 gap-6">
          <div className="min-w-0">
            <ChartContainer config={qcChartConfig} className="h-96 w-full">
              <BarChart data={data}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis dataKey="name" tick={{ fontSize: 11 }} />
                <YAxis domain={[0, 100]} />
                <ChartTooltip cursor={false} content={<ChartTooltipContent />} />
                <Bar dataKey="genomeComplete" fill="var(--chart-2)" name="Genome complete %" radius={[4, 4, 0, 0]} />
                <Bar dataKey="proteinComplete" fill="var(--chart-4)" name="Protein complete %" radius={[4, 4, 0, 0]} />
              </BarChart>
            </ChartContainer>
          </div>

          <div className="max-h-96 overflow-auto rounded-md border">
            <Table>
              <TableHeader className="sticky top-0 bg-card">
                <TableRow>
                  <TableHead>Sample</TableHead>
                  <TableHead className="text-right">Genome C</TableHead>
                  <TableHead className="text-right">Protein C</TableHead>
                  <TableHead className="text-right">D</TableHead>
                  <TableHead className="text-right">Delta</TableHead>
                  <TableHead>Flags</TableHead>
                </TableRow>
              </TableHeader>
              <TableBody>
                {data.map((item) => (
                  <TableRow key={item.name}>
                    <TableCell className="max-w-64 truncate">{item.name}</TableCell>
                    <TableCell className="text-right">{formatNumber(item.genomeComplete)}%</TableCell>
                    <TableCell className="text-right">{formatNumber(item.proteinComplete)}%</TableCell>
                    <TableCell className="text-right">{formatNumber(item.duplicated)}%</TableCell>
                    <TableCell className="text-right">{formatNumber(item.difference)}%</TableCell>
                    <TableCell>
                      <div className="flex flex-wrap gap-1">
                        {item.differenceFlagged && (
                          <Badge variant="destructive">Delta above avg</Badge>
                        )}
                        {item.duplicationFlagged && (
                          <Badge variant="destructive">Duplication above avg</Badge>
                        )}
                        {!item.differenceFlagged && !item.duplicationFlagged && (
                          <Badge variant="outline">OK</Badge>
                        )}
                      </div>
                    </TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
          </div>

          <div className="max-h-72 overflow-auto rounded-md border p-3">
            <div className="grid grid-cols-1 gap-2 text-xs text-muted-foreground sm:grid-cols-2 lg:grid-cols-4">
                {data.map((item) => (
                  <div key={`${item.name}-breakdown`} className="rounded-md border p-2">
                    <div className="mb-1 truncate font-medium text-foreground">{item.name}</div>
                    <div>C {formatNumber(item.complete)}%</div>
                    <div>D {formatNumber(item.duplicated)}%</div>
                    <div>F {formatNumber(item.fragmented)}%</div>
                    <div>M {formatNumber(item.missing)}%</div>
                  </div>
                ))}
              </div>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

function MetricRows({ rows }: { rows: Array<[string, string]> }) {
  return (
    <div className="divide-y">
      {rows.map(([label, value]) => (
        <div key={label} className="flex items-center justify-between gap-6 py-3">
          <span className="text-sm text-muted-foreground">{label}</span>
          <span className="text-right text-sm font-medium">{value}</span>
        </div>
      ))}
    </div>
  );
}
