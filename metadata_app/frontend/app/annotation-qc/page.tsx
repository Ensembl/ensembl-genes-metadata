"use client";

import React, {useEffect, useMemo, useState} from "react";
import {Bar, BarChart, CartesianGrid, Scatter, ScatterChart, XAxis, YAxis,} from "recharts";
import {
  Activity,
  ArrowRight,
  ChevronLeft,
  ChevronRight,
  Dna,
  Download,
  Info,
  Layers3,
  Search,
  ShieldAlert,
  ShieldCheck
} from "lucide-react";
import {Badge} from "@/components/ui/badge";
import {Button} from "@/components/ui/button";
import {Card, CardContent, CardDescription, CardHeader, CardTitle} from "@/components/ui/card";
import {Input} from "@/components/ui/input";
import {Label} from "@/components/ui/label";
import MultipleSelector, {Option} from "@/components/ui/multi_select";
import {Select, SelectContent, SelectItem, SelectTrigger, SelectValue} from "@/components/ui/select";
import {ChartConfig, ChartContainer, ChartTooltip, ChartTooltipContent} from "@/components/ui/chart";
import {Table, TableBody, TableCell, TableHead, TableHeader, TableRow} from "@/components/ui/table";
import {Tooltip, TooltipContent, TooltipTrigger} from "@/components/ui/tooltip";

type GenebuildStatus = "live" | "in progress" | "handover" | "completed" | "blocked";
type ScreeningFilter = "all" | "review";
type SingleReportTab = "Overview" | "Gene structure" | "Quality & validation" | "All metrics";

type QcRecord = {
    gca: string;
    species: string;
    taxonomy: string;
    taxonId: string;
    status: GenebuildStatus;
    associatedProject: string | null;
    groupName: string | null;
    infraName: string | null;
    genebuilder: string | null;
    annotationSource: string | null;
    annotationMethod: string | null;
    dateStarted: string | null;
    releaseDate: string | null;
    dateStatusUpdate: string | null;
    lastGenebuildUpdate: string | null;
    assemblyBuscoLineage: string | null;
    assemblyBuscoVersion: string | null;
    proteinBuscoLineage: string | null;
    proteinBuscoVersion: string | null;
    contigN50Mb: number | null;
    scaffoldN50Mb: number | null;
    genomeSizeGb: number | null;
    gcPercent: number | null;
    geneCount: number | null;
    proteinCoding: number | null;
    transcriptCount: number | null;
    genomeBuscoComplete: number | null;
    genomeBuscoDuplicated: number | null;
    genomeBuscoFragmented: number | null;
    genomeBuscoMissing: number | null;
    proteinBuscoComplete: number | null;
    buscoComplete: number | null;
    buscoDuplicated: number | null;
    buscoFragmented: number | null;
    buscoMissing: number | null;
    repeatPercent: number | null;
    annotationScore: number | null;
    assemblyScore: number | null;
    lncRnaCount: number | null;
    pseudogeneCount: number | null;
    smallRnaCount: number | null;
    newMetrics: Record<string, number | string | null>;
    pcoa1: number;
    pcoa2: number;
    outlierScore: number;
};

type CohortMeans = Record<string, number | null>;

const statusOptions: Array<"all" | GenebuildStatus> = [
    "all",
    "live",
    "in progress",
    "handover",
    "completed",
    "blocked",
];

const annotationMethodOptions = ["ensembl", "projection_build", "full_genebuild", "anno", "helixer", "braker"] as const;

const chartConfig = {
    count: {label: "Annotations", color: "var(--chart-2)"},
    flagged: {label: "Needs review", color: "var(--destructive)"},
    genome: {label: "Genome BUSCO", color: "var(--chart-2)"},
    protein: {label: "Protein BUSCO", color: "var(--chart-4)"},
} satisfies ChartConfig;

const PAGE_SIZE = 25;

const numbers = (values: Array<number | null | undefined>) =>
    values.filter((value): value is number => typeof value === "number" && Number.isFinite(value));

const quantile = (values: Array<number | null | undefined>, fraction: number) => {
    const valid = numbers(values).sort((a, b) => a - b);
    if (!valid.length) return 0;
    const position = (valid.length - 1) * fraction;
    const lower = Math.floor(position);
    const upper = Math.ceil(position);
    return valid[lower] + (valid[upper] - valid[lower]) * (position - lower);
};

const median = (values: Array<number | null | undefined>) => quantile(values, 0.5);

const formatNumber = (value: number | null | undefined, digits = 1) =>
    value === null || value === undefined || !Number.isFinite(value)
        ? "n/a"
        : new Intl.NumberFormat("en-GB", {
            maximumFractionDigits: digits,
            minimumFractionDigits: digits,
        }).format(value);

const formatCompact = (value: number) =>
    new Intl.NumberFormat("en-GB", {notation: "compact", maximumFractionDigits: 1}).format(value);

const formatMetricLabel = (name: string) =>
    name.replace(/^Genebuild\.Stats\./i, "").replace(/^Stats\./i, "").replace(/[_-]+/g, " ").replace(/\b\w/g, (letter) => letter.toUpperCase());

const formatMetricValue = (value: number | string | null) => {
    if (typeof value === "number") return formatNumber(value, 2);
    if (typeof value === "string" && value.trim() !== "" && Number.isFinite(Number(value))) {
        return formatNumber(Number(value), 2);
    }
    return value ?? "n/a";
};

const getStatusBadge = (status: GenebuildStatus) => {
    if (status === "blocked") return "destructive";
    if (status === "live") return "default";
    return "secondary";
};

const buscoDelta = (record: QcRecord) =>
    record.genomeBuscoComplete === null || record.proteinBuscoComplete === null
        ? null
        : Math.abs(record.genomeBuscoComplete - record.proteinBuscoComplete);

const needsReview = (record: QcRecord) => record.outlierScore >= 2;

const buscoDeltaThreshold = (records: QcRecord[]) => {
    const deltas = records.map(buscoDelta).filter((value): value is number => value !== null && Number.isFinite(value));
    if (deltas.length < 2) return Number.POSITIVE_INFINITY;
    const mean = deltas.reduce((total, value) => total + value, 0) / deltas.length;
    const standardDeviation = Math.sqrt(deltas.reduce((total, value) => total + (value - mean) ** 2, 0) / deltas.length);
    return mean + 2 * standardDeviation;
};

const needsBuscoReview = (record: QcRecord, threshold: number) => {
    const delta = buscoDelta(record);
    return delta !== null && delta > threshold;
};

const OUTLIER_SCORE_HELP =
    "This score compares the annotation with live annotations in the same clade. It combines standardized differences in protein BUSCO completeness, genome BUSCO completeness, gene count, and transcript count. A score near 0 means no measurable difference or insufficient live peers; a score of 2 or more flags the annotation for review. It is a screening signal, not a pass/fail judgement.";

const reviewReasons = (record: QcRecord) => {
    const reasons: string[] = [];
    if (record.outlierScore >= 2) reasons.push("multivariate outlier");
    if ((record.buscoMissing ?? 0) >= 5) reasons.push("BUSCO missing");
    if ((record.buscoDuplicated ?? 0) >= 5) reasons.push("BUSCO duplication");
    if ((buscoDelta(record) ?? 0) >= 5) reasons.push("BUSCO delta");
    return reasons.length ? reasons : ["within cohort range"];
};

const exportHeader = (key: string) =>
    key.replace(/([a-z])([A-Z])/g, "$1_$2").replace(/\./g, "_").toLowerCase();

const downloadQcCsv = (records: QcRecord[]) => {
    if (!records.length) return;
    const preferredColumns = [
        "gca", "species", "taxonId", "taxonomy", "status", "proteinBuscoComplete",
        "buscoMissing", "outlierScore", "reviewSignals",
    ];
    const metricColumns = [...new Set(records.flatMap((record) => Object.keys(record.newMetrics)))].sort();
    const baseColumns = [...new Set(records.flatMap((record) => Object.keys(record)))].filter(
        (key) => key !== "newMetrics" && !preferredColumns.includes(key),
    );
    const columns = [...preferredColumns, ...baseColumns, ...metricColumns.map((key) => `newMetrics.${key}`)];
    const valueFor = (record: QcRecord, column: string) => {
        if (column === "reviewSignals") return reviewReasons(record).join("; ");
        if (column.startsWith("newMetrics.")) return record.newMetrics[column.slice("newMetrics.".length)] ?? "";
        return (record as unknown as Record<string, unknown>)[column] ?? "";
    };
    const escape = (value: unknown) => `"${String(value).replace(/"/g, '""')}"`;
    const csv = [
        columns.map((column) => escape(exportHeader(column))).join(","),
        ...records.map((record) => columns.map((column) => escape(valueFor(record, column))).join(",")),
    ].join("\n");
    const url = URL.createObjectURL(new Blob([csv], {type: "text/csv;charset=utf-8"}));
    const link = document.createElement("a");
    link.href = url;
    link.download = `annotation-qc-${new Date().toISOString().slice(0, 10)}.csv`;
    link.click();
    URL.revokeObjectURL(url);
};

const getRecordText = (record: QcRecord) =>
    `${record.gca} ${record.species} ${record.taxonomy} ${record.taxonId}`.toLowerCase();

export default function AnnotationQcPage() {
    const [gca, setGca] = useState("");
    const [taxonId, setTaxonId] = useState("");
    const [clade, setClade] = useState<Option[]>([]);
    const [cladeOptions, setCladeOptions] = useState<Option[]>([]);
    const [status, setStatus] = useState<(typeof statusOptions)[number]>("all");
    const [annotationMethod, setAnnotationMethod] = useState<(typeof annotationMethodOptions)[number]>("ensembl");
    const [submitted, setSubmitted] = useState({
        gca: "",
        taxonId: "",
        clade: "",
        status: "all",
        annotationMethod: "ensembl"
    });
    const [records, setRecords] = useState<QcRecord[]>([]);
    const [comparisonRecords, setComparisonRecords] = useState<QcRecord[]>([]);
    const [cohortMeans, setCohortMeans] = useState<CohortMeans>({});
    const [loading, setLoading] = useState(false);
    const [errorMessage, setErrorMessage] = useState<string | null>(null);
    const [hasRun, setHasRun] = useState(false);
    const [search, setSearch] = useState("");
    const [screeningFilter, setScreeningFilter] = useState<ScreeningFilter>("all");
    const [sortBy, setSortBy] = useState("outlierScore");
    const [page, setPage] = useState(1);

    const loadRecords = async (
        requestedGca: string,
        requestedTaxonIds: number[] | null,
        requestedClade: string,
        requestedStatus: (typeof statusOptions)[number],
        requestedAnnotationMethod: string,
    ) => {
        setLoading(true);
        setErrorMessage(null);
        setRecords([]);
        setComparisonRecords([]);
        setCohortMeans({});
        try {
            const response = await fetch("/api/qc/annotation/filter", {
                method: "POST",
                headers: {"Content-Type": "application/json", accept: "application/json"},
                body: JSON.stringify({
                    gca: requestedGca.trim() ? [requestedGca.trim()] : null,
                    taxon_id: requestedTaxonIds,
                    clade: requestedClade || null,
                    status: requestedStatus === "all" ? null : requestedStatus,
                    annotation_method: requestedAnnotationMethod ? [requestedAnnotationMethod] : null,
                }),
            });
            if (!response.ok) throw new Error(await response.text());
            const result = await response.json();
            setRecords(result.records ?? []);
            setComparisonRecords(result.cohort_records ?? result.records ?? []);
            setCohortMeans(result.cohort_means ?? {});
            setPage(1);
        } catch (error) {
            setErrorMessage(error instanceof Error ? error.message : String(error));
            setRecords([]);
        } finally {
            setLoading(false);
        }
    };

    useEffect(() => {
        void fetch("/api/qc/clades")
            .then(async (response) => {
                if (!response.ok) throw new Error(await response.text());
                return response.json();
            })
            .then((result) => setCladeOptions((result.clades ?? []).map((value: string) => ({value, label: value}))))
            .catch((error) => setErrorMessage(error instanceof Error ? error.message : String(error)));
    }, []);

    const cohortRecords = useMemo(() => {
        const cladeQuery = submitted.clade.toLowerCase();
        const gcaQuery = submitted.gca.trim().toLowerCase();

        return records.filter(
            (record) =>
                (submitted.status === "all" || record.status === submitted.status) &&
                (!gcaQuery || record.gca.toLowerCase() === gcaQuery) &&
                (!cladeQuery || record.taxonomy.toLowerCase() === cladeQuery),
        );
    }, [records, submitted]);

    const reviewCount = useMemo(() => cohortRecords.filter(needsReview).length, [cohortRecords]);
    const singleRecord = submitted.gca.trim() ? cohortRecords[0] : null;

    const metricSummary = useMemo(() => {
        const metrics = [
            ["Protein BUSCO complete", cohortRecords.map((record) => record.proteinBuscoComplete), "%"],
            ["Genome BUSCO complete", cohortRecords.map((record) => record.genomeBuscoComplete), "%"],
            ["BUSCO missing", cohortRecords.map((record) => record.buscoMissing), "%"],
            ["Genes", cohortRecords.map((record) => record.geneCount), ""],
            ["Transcripts", cohortRecords.map((record) => record.transcriptCount), ""],
        ] as const;

        return metrics.map(([label, values, suffix]) => ({
            label,
            median: median(values),
            p10: quantile(values, 0.1),
            p90: quantile(values, 0.9),
            suffix,
        }));
    }, [cohortRecords]);

    const reviewQueue = useMemo(() => {
        const query = search.trim().toLowerCase();
        const filtered = cohortRecords.filter(
            (record) => (!query || getRecordText(record).includes(query)) && (screeningFilter === "all" || needsReview(record)),
        );
        return [...filtered].sort((a, b) => {
            if (sortBy === "buscoMissing") return (b.buscoMissing ?? -1) - (a.buscoMissing ?? -1);
            if (sortBy === "species") return a.species.localeCompare(b.species);
            return b.outlierScore - a.outlierScore;
        });
    }, [cohortRecords, screeningFilter, search, sortBy]);

    const pageCount = Math.max(1, Math.ceil(reviewQueue.length / PAGE_SIZE));
    const visibleQueue = reviewQueue.slice((page - 1) * PAGE_SIZE, page * PAGE_SIZE);

    const submitFilters = async () => {
        const requestedTaxonIds = taxonId.split(",").map((value) => Number.parseInt(value.trim(), 10)).filter((value) => Number.isFinite(value));
        const requestedClade = clade[0]?.value ?? "";
        const isCohortQuery = !gca.trim() && Boolean(taxonId.trim() || requestedClade);
        const requestedStatus = isCohortQuery ? status : "all";
        const requestedAnnotationMethod = isCohortQuery ? annotationMethod : "";
        setSubmitted({
            gca,
            taxonId,
            clade: requestedClade,
            status: requestedStatus,
            annotationMethod: requestedAnnotationMethod
        });
        setHasRun(true);
        await loadRecords(gca, requestedTaxonIds.length ? requestedTaxonIds : null, requestedClade, requestedStatus, requestedAnnotationMethod);
    };

    const primaryFilterCount = [gca.trim(), taxonId.trim(), clade[0]?.value ?? ""].filter(Boolean).length;
    const showAdditionalFilters = !gca.trim() && Boolean(taxonId.trim() || clade[0]?.value);
    const primaryFilterMessage = primaryFilterCount === 0
        ? "Choose one primary filter to begin."
        : primaryFilterCount > 1
            ? "Please fill out only one primary filter."
            : "";

    return (
        <div className="min-h-screen flex justify-center px-4 py-10 sm:px-6 lg:px-8">
            <div className="w-full max-w-7xl">
                <Card className="bg-secondary shadow-lg">
                    <CardHeader className="mb-4">
                        <CardTitle>Annotation QC screening</CardTitle>
                        <CardDescription>Screen a clade or taxonomy cohort for annotations that differ from their peers
                            or get a single annotation report of all metrics.</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <div className="space-y-6">
                            <section aria-labelledby="primary-query-heading">
                                <div className="mb-4">
                                    <div className="flex items-center gap-2"><h3 id="primary-query-heading"
                                                                                 className="text-sm font-semibold">Primary
                                        query</h3><Tooltip><TooltipTrigger asChild><Button type="button" variant="ghost"
                                                                                           size="icon"
                                                                                           className="h-5 w-5"
                                                                                           aria-label="Explain QC cohort selection"><Info
                                        className="h-3.5 w-3.5"/></Button></TooltipTrigger><TooltipContent
                                        className="max-w-xs leading-relaxed">GCA searches use live annotations in the
                                        GCA&apos;s configured clade as the comparison cohort. Clade searches use the
                                        selected clade cohort. Taxon ID searches use all annotations returned under that
                                        taxon ID, including descendants, as one taxon cohort.</TooltipContent></Tooltip>
                                    </div>
                                    <p className="mt-1 text-sm text-muted-foreground">Choose exactly one: enter a GCA
                                        accession, a taxon ID, or select a clade. Leave the other two empty.</p></div>
                                <div className="grid grid-cols-1 gap-4 md:grid-cols-3">
                                    <div><Label htmlFor="gca-accession">GCA accession</Label><Input id="gca-accession"
                                                                                                    className="mt-3 bg-filter-input-bg dark:bg-transparent"
                                                                                                    placeholder="GCA_000001405.29"
                                                                                                    value={gca}
                                                                                                    onChange={(event) => setGca(event.target.value)}/>
                                    </div>
                                    <div><Label htmlFor="taxon-id-filter">Taxon ID</Label><Input id="taxon-id-filter"
                                                                                                 className="mt-3 bg-filter-input-bg dark:bg-transparent"
                                                                                                 placeholder="9606"
                                                                                                 value={taxonId}
                                                                                                 onChange={(event) => setTaxonId(event.target.value)}/>
                                    </div>
                                    <div><Label>Clade</Label><MultipleSelector
                                        className="mt-3 min-h-8 bg-filter-input-bg dark:bg-transparent"
                                        commandProps={{className: "bg-transparent"}} value={clade}
                                        options={cladeOptions} defaultOptions={cladeOptions} maxSelected={1}
                                        hidePlaceholderWhenSelected placeholder="Type to find a clade..."
                                        triggerSearchOnFocus onSearchSync={(value) => {
                                        const query = value.trim().toLowerCase();
                                        return cladeOptions.filter((option) => option.label.toLowerCase().includes(query));
                                    }} onChange={setClade}/></div>
                                </div>
                                {primaryFilterMessage && <p className="mt-3 text-sm text-destructive"
                                                            role="alert">{primaryFilterMessage}</p>}
                            </section>
                            {showAdditionalFilters &&
                                <section aria-labelledby="additional-filters-heading" className="border-t pt-5">
                                    <div className="mb-4"><h3 id="additional-filters-heading"
                                                              className="text-sm font-semibold">Additional filters</h3>
                                        <p className="mt-1 text-sm text-muted-foreground">Optionally narrow the taxon or
                                            clade cohort by genebuild status and annotation method.</p></div>
                                    <div className="grid grid-cols-1 gap-4 md:grid-cols-2">
                                        <div><Label>Genebuild status</Label><Select value={status}
                                                                                    onValueChange={(value) => setStatus(value as typeof status)}><SelectTrigger
                                            className="mt-3 w-full bg-filter-input-bg dark:bg-transparent"><SelectValue
                                            placeholder="Select status"/></SelectTrigger><SelectContent>{statusOptions.map((option) =>
                                            <SelectItem key={option}
                                                        value={option}>{option === "all" ? "All statuses" : option}</SelectItem>)}</SelectContent></Select>
                                        </div>
                                        <div><Label>Annotation method</Label><Select value={annotationMethod}
                                                                                     onValueChange={(value) => setAnnotationMethod(value as (typeof annotationMethodOptions)[number])}><SelectTrigger
                                            className="mt-3 w-full bg-filter-input-bg dark:bg-transparent"><SelectValue/></SelectTrigger><SelectContent>{annotationMethodOptions.map((option) =>
                                            <SelectItem key={option}
                                                        value={option}>{option}</SelectItem>)}</SelectContent></Select>
                                        </div>
                                    </div>
                                </section>}
                        </div>
                        <div className="mt-6 flex justify-end"><Button className="w-full sm:w-auto"
                                                                       onClick={submitFilters}
                                                                       disabled={loading || primaryFilterCount !== 1}><Search
                            className="mr-1 h-4 w-4"/>Run QC</Button></div>
                        {errorMessage && <p className="mt-4 text-sm text-destructive" role="alert">{errorMessage}</p>}
                    </CardContent>
                </Card>

                {hasRun && loading &&
                    <Card className="mt-8"><CardHeader><CardTitle>Running QC</CardTitle><CardDescription>Retrieving
                        annotations and calculating quality metrics...</CardDescription></CardHeader></Card>}
                {hasRun && !loading && !errorMessage && submitted.gca.trim() && !singleRecord &&
                    <Card className="mt-8"><CardHeader><CardTitle>No GCA match</CardTitle><CardDescription>No annotation
                        matched that accession and filter.</CardDescription></CardHeader></Card>}
                {hasRun && !loading && !errorMessage && !submitted.gca.trim() &&
                    <CohortDashboard records={cohortRecords} reviewCount={reviewCount} metricSummary={metricSummary}
                                     cohortMeans={cohortMeans} multivariateRecords={cohortRecords}
                                     selectedStatus={submitted.status} search={search} setSearch={(value) => {
                        setSearch(value);
                        setPage(1);
                    }} screeningFilter={screeningFilter} setScreeningFilter={(value) => {
                        setScreeningFilter(value);
                        setPage(1);
                    }} sortBy={sortBy} setSortBy={(value) => {
                        setSortBy(value);
                        setPage(1);
                    }} page={page} pageCount={pageCount} reviewQueueLength={reviewQueue.length}
                                     visibleQueue={visibleQueue} downloadRecords={reviewQueue} setPage={setPage}/>}
                {hasRun && !loading && !errorMessage && singleRecord &&
                    <SingleRecord record={singleRecord} cohortMeans={cohortMeans} cohortRecords={comparisonRecords}/>}
            </div>
        </div>
    );
}

function PointLegend({flaggedLabel = "Needs review (outlier score ≥ 2)", flaggedColor = "var(--destructive)"}: { flaggedLabel?: string; flaggedColor?: string }) {
    return <div className="mt-3 flex flex-wrap gap-x-5 gap-y-2 text-xs text-muted-foreground"><span
        className="flex items-center gap-2"><span className="h-2.5 w-2.5 rounded-full"
                                                  style={{backgroundColor: "var(--chart-2)"}}/>Within cohort range</span><span
        className="flex items-center gap-2"><span className="h-2.5 w-2.5 rounded-full"
                                                  style={{backgroundColor: flaggedColor}}/>{flaggedLabel}</span>
    </div>;
}

function MultivariateTooltip({active, payload}: { active?: boolean; payload?: Array<{ payload?: QcRecord }> }) {
    const record = payload?.[0]?.payload;
    if (!active || !record) return null;
    return <div className="min-w-52 rounded-lg border bg-background px-3 py-2 text-xs shadow-xl">
        <div className="font-mono font-medium">{record.gca}</div>
        <div className="mt-1 text-muted-foreground">{record.species}</div>
        <div className="mt-2 grid gap-1">
            <span>Profile axis 1: <strong>{formatNumber(record.pcoa1, 2)}</strong></span><span>Profile axis 2: <strong>{formatNumber(record.pcoa2, 2)}</strong></span><span>Outlier score: <strong>{formatNumber(record.outlierScore, 2)}</strong></span>
        </div>
    </div>;
}

function BuscoTooltip({active, payload, threshold}: {
    active?: boolean;
    payload?: Array<{ payload?: QcRecord }>;
    threshold: number
}) {
    const record = payload?.[0]?.payload;
    if (!active || !record) return null;
    return <div className="min-w-52 rounded-lg border bg-background px-3 py-2 text-xs shadow-xl">
        <div className="font-mono font-medium">{record.gca}</div>
        <div className="mt-2 grid gap-1">
            <span>Genome BUSCO: <strong>{formatNumber(record.genomeBuscoComplete)}%</strong></span><span>Protein BUSCO: <strong>{formatNumber(record.proteinBuscoComplete)}%</strong></span><span>Completeness delta: <strong>{formatNumber(buscoDelta(record), 2)}%</strong></span><span>Cohort flag threshold: <strong>{formatNumber(threshold, 2)}%</strong></span>
        </div>
    </div>;
}

function CohortDashboard({
                             records,
                             reviewCount,
                             metricSummary,
                             cohortMeans,
                             multivariateRecords,
                             selectedStatus,
                             search,
                             setSearch,
                             screeningFilter,
                             setScreeningFilter,
                             sortBy,
                             setSortBy,
                             page,
                             pageCount,
                             reviewQueueLength,
                             visibleQueue,
                             downloadRecords,
                             setPage,
                         }: {
    records: QcRecord[];
    reviewCount: number;
    metricSummary: Array<{ label: string; median: number; p10: number; p90: number; suffix: string }>;
    cohortMeans: CohortMeans;
    multivariateRecords: QcRecord[];
    selectedStatus: string;
    search: string;
    setSearch: (value: string) => void;
    screeningFilter: ScreeningFilter;
    setScreeningFilter: (value: ScreeningFilter) => void;
    sortBy: string;
    setSortBy: (value: string) => void;
    page: number;
    pageCount: number;
    reviewQueueLength: number;
    visibleQueue: QcRecord[];
    downloadRecords: QcRecord[];
    setPage: (value: number) => void;
}) {
    const reviewPercent = records.length ? (reviewCount / records.length) * 100 : 0;
    const buscoThreshold = buscoDeltaThreshold(records);
    return (
        <div className="mt-8 space-y-6" aria-live="polite">
            <div className="flex flex-col gap-2 sm:flex-row sm:items-end sm:justify-between">
                <div><h2 className="text-xl font-semibold">Cohort screening</h2><p
                    className="text-sm text-muted-foreground">{records.length.toLocaleString("en-GB")} annotations
                    loaded. Review threshold: multivariate score ≥ 2.</p></div>
                <Badge variant={reviewCount ? "destructive" : "outline"}>{reviewCount.toLocaleString("en-GB")} need
                    review</Badge></div>
            <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4"><MetricCard label="Annotations"
                                                                                              value={records.length}/><MetricCard
                label="Needs review" value={reviewCount}
                suffix={records.length ? `${formatNumber(reviewPercent)}%` : ""}/><MetricCard
                label="Median protein BUSCO" value={median(records.map((record) => record.proteinBuscoComplete))}
                suffix="%"/><MetricCard label="Genebuild status"
                                        value={selectedStatus === "all" ? "All statuses" : selectedStatus}/></div>

            <div className="grid grid-cols-1 gap-6 lg:grid-cols-2">
                <Card>
                    <CardHeader>
                        <PlotTitle title="Multivariate annotation profile"
                                   help="Each dot is an annotation positioned by the two coordinates calculated from its clade-standardised metric profile. Dots that are far from the main cloud are more unusual in the multivariate analysis."/>
                        <CardDescription>All annotations are shown. Red points meet the outlier
                            threshold.</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <ChartContainer config={chartConfig} className="h-72 w-full mb-2">
                            <ScatterChart accessibilityLayer>
                                <CartesianGrid strokeDasharray="3 3"/>
                                <XAxis type="number" dataKey="pcoa1" name="Profile axis 1"/>
                                <YAxis type="number" dataKey="pcoa2" name="Profile axis 2"/>
                                <ChartTooltip cursor={{strokeDasharray: "3 3"}} content={<MultivariateTooltip/>}/>
                                <Scatter data={multivariateRecords.filter((record) => !needsReview(record))}
                                         fill="var(--chart-2)" name="Within range"/>
                                <Scatter data={multivariateRecords.filter(needsReview)} fill="var(--destructive)"
                                         name="Needs review"/>
                            </ScatterChart></ChartContainer><PointLegend/>
                    </CardContent>
                </Card>
                <Card>
                    <CardHeader>
                        <PlotTitle title="BUSCO concordance"
                                   help="Each point is one annotation. The x-axis is genome BUSCO completeness and the y-axis is protein BUSCO completeness. Red points have an unusually large absolute difference between annotation and assembly completeness for this cohort; high completeness alone is not flagged. Hover over a point to see its GCA."/>
                        <CardDescription>Genome versus protein completeness. Red points exceed the cohort mean BUSCO
                            delta by more than two standard deviations.</CardDescription>
                    </CardHeader>
                    <CardContent>
                        <ChartContainer config={chartConfig} className="h-72 w-full">
                            <ScatterChart accessibilityLayer>
                                <CartesianGrid strokeDasharray="3 3"/>
                                <XAxis type="number" dataKey="genomeBuscoComplete" name="Genome BUSCO"
                                       domain={["auto", "auto"]} unit="%"/>
                                <YAxis type="number" dataKey="proteinBuscoComplete" name="Protein BUSCO"
                                       domain={["auto", "auto"]} unit="%"/><ChartTooltip
                                cursor={{strokeDasharray: "3 3"}} content={<BuscoTooltip threshold={buscoThreshold}/>}/>
                                <Scatter
                                    data={multivariateRecords.filter((record) => !needsBuscoReview(record, buscoThreshold))}
                                    fill="var(--chart-2)" name="BUSCO delta within range"/>
                                <Scatter
                                    data={multivariateRecords.filter((record) => needsBuscoReview(record, buscoThreshold))}
                                    fill="var(--destructive)" name="Unusual BUSCO delta"/>
                            </ScatterChart></ChartContainer><PointLegend
                        flaggedLabel="Unusual BUSCO delta (mean + 2 SD)"/>
                    </CardContent>
                </Card>
                <Card><CardHeader><PlotTitle title="Metric reference"
                                             help="P10 and P90 show the central 80% range of the current cohort, with the median in the middle. These values describe the cohort and are not universal pass/fail limits."/><CardDescription>Cohort
                    median and central 80% range for screening
                    context.</CardDescription></CardHeader><CardContent><Table><TableHeader><TableRow><TableHead>Metric</TableHead><TableHead
                    className="text-right">P10</TableHead><TableHead className="text-right">Median</TableHead><TableHead
                    className="text-right">P90</TableHead></TableRow></TableHeader><TableBody>{metricSummary.map((item) =>
                    <TableRow key={item.label}><TableCell>{item.label}</TableCell><TableCell
                        className="text-right">{formatNumber(item.p10)}{item.suffix}</TableCell><TableCell
                        className="text-right font-medium">{formatNumber(item.median)}{item.suffix}</TableCell><TableCell
                        className="text-right">{formatNumber(item.p90)}{item.suffix}</TableCell></TableRow>)}</TableBody></Table></CardContent></Card>
                <Card><CardHeader><CardTitle>Mean feature counts</CardTitle><CardDescription>Mean gene-structure counts
                    in the current
                    cohort.</CardDescription></CardHeader><CardContent><Table><TableHeader><TableRow><TableHead>Feature</TableHead><TableHead
                    className="text-right">Mean</TableHead></TableRow></TableHeader><TableBody>{[["Genes", "geneCount"], ["Protein coding", "proteinCoding"], ["Transcripts", "transcriptCount"], ["lncRNA", "lncRnaCount"], ["Pseudogenes", "pseudogeneCount"], ["Small RNA", "smallRnaCount"]].map(([label, metric]) =>
                    <TableRow key={metric}><TableCell>{label}</TableCell><TableCell
                        className="text-right font-medium">{formatNumber(cohortMeans[metric])}</TableCell></TableRow>)}</TableBody></Table></CardContent></Card>
            </div>

            <Card>
                <CardHeader>
                    <div className="flex flex-col gap-3 lg:flex-row lg:items-end lg:justify-between">
                        <div>
                            <CardTitle className="flex items-center gap-2"><ShieldAlert className="h-5 w-5"/>Review
                                queue</CardTitle>
                            <CardDescription>Search and rank records without rendering the full
                                cohort.</CardDescription>
                        </div>
                        <div className="flex flex-col gap-3 sm:flex-row">
                            <Input aria-label="Search annotations" className="w-full sm:w-64"
                                   placeholder="Search GCA, species, clade..." value={search}
                                   onChange={(event) => setSearch(event.target.value)}/>
                            <Select value={screeningFilter}
                                    onValueChange={(value) => setScreeningFilter(value as ScreeningFilter)}><SelectTrigger
                                className="w-full sm:w-40"><SelectValue/></SelectTrigger><SelectContent><SelectItem
                                value="all">All records</SelectItem><SelectItem value="review">Needs review</SelectItem></SelectContent></Select>
                            <Select value={sortBy} onValueChange={setSortBy}><SelectTrigger
                                className="w-full sm:w-44"><SelectValue/></SelectTrigger><SelectContent><SelectItem
                                value="outlierScore">Highest score</SelectItem><SelectItem value="buscoMissing">BUSCO
                                missing</SelectItem><SelectItem value="species">Species A–Z</SelectItem></SelectContent></Select>
                            <Button variant="outline" onClick={() => downloadQcCsv(downloadRecords)}
                                    disabled={!downloadRecords.length}><Download className="mr-1 h-4 w-4"/>Download CSV</Button>
                        </div>
                    </div>
                </CardHeader>
                <CardContent>
                    <div className="overflow-x-auto">
                        <Table><TableHeader><TableRow><TableHead>GCA</TableHead><TableHead>Species</TableHead><TableHead>Taxon
                            ID</TableHead><TableHead>Clade</TableHead><TableHead>Status</TableHead><TableHead
                            className="text-right">Protein BUSCO</TableHead><TableHead
                            className="text-right">Missing</TableHead><TableHead
                            className="text-right">Score</TableHead><TableHead>Signals</TableHead></TableRow></TableHeader><TableBody>{visibleQueue.map((record) =>
                            <TableRow key={`${record.gca}-${record.annotationMethod ?? ""}`}><TableCell
                                className="font-mono text-xs whitespace-nowrap">{record.gca}</TableCell><TableCell
                                className="max-w-56 truncate">{record.species}</TableCell><TableCell>{record.taxonId}</TableCell><TableCell>{record.taxonomy}</TableCell><TableCell><Badge
                                variant={getStatusBadge(record.status)}>{record.status}</Badge></TableCell><TableCell
                                className="text-right">{formatNumber(record.proteinBuscoComplete)}%</TableCell><TableCell
                                className="text-right">{formatNumber(record.buscoMissing)}%</TableCell><TableCell
                                className="text-right font-medium">{formatNumber(record.outlierScore)}</TableCell><TableCell
                                className="max-w-52 text-xs text-muted-foreground">{reviewReasons(record).join(", ")}</TableCell></TableRow>)}{!visibleQueue.length &&
                            <TableRow><TableCell colSpan={9} className="h-24 text-center text-muted-foreground">No
                                records match the current queue filters.</TableCell></TableRow>}</TableBody></Table>
                    </div>
                    <div className="mt-4 flex items-center justify-between text-sm text-muted-foreground">
                        <span>{reviewQueueLength.toLocaleString("en-GB")} records in queue · page {page} of {pageCount}</span>
                        <div className="flex gap-2"><Button variant="outline" size="sm"
                                                            onClick={() => setPage(Math.max(1, page - 1))}
                                                            disabled={page <= 1}><ChevronLeft className="h-4 w-4"/>Previous</Button><Button
                            variant="outline" size="sm" onClick={() => setPage(Math.min(pageCount, page + 1))}
                            disabled={page >= pageCount}>Next<ChevronRight className="h-4 w-4"/></Button></div>
                    </div>
                </CardContent>
            </Card>
        </div>
    );
}

function SingleRecord({record, cohortMeans, cohortRecords}: {
    record: QcRecord;
    cohortMeans: CohortMeans;
    cohortRecords: QcRecord[]
}) {
    const [tab, setTab] = useState<SingleReportTab>("Overview");
    const buscoProfile = [
        {
            metric: "Assembly BUSCO",
            singleCopy: Math.max(0, (record.genomeBuscoComplete ?? 0) - (record.genomeBuscoDuplicated ?? 0)),
            duplicated: record.genomeBuscoDuplicated ?? 0,
            fragmented: record.genomeBuscoFragmented ?? 0,
            missing: record.genomeBuscoMissing ?? 0
        },
        {
            metric: "Annotation BUSCO",
            singleCopy: Math.max(0, (record.buscoComplete ?? 0) - (record.buscoDuplicated ?? 0)),
            duplicated: record.buscoDuplicated ?? 0,
            fragmented: record.buscoFragmented ?? 0,
            missing: record.buscoMissing ?? 0
        },
    ];
    const countProfile = [
        {metric: "Genes", value: record.geneCount ?? 0},
        {metric: "Protein coding", value: record.proteinCoding ?? 0},
        {metric: "Transcripts", value: record.transcriptCount ?? 0},
        {metric: "lncRNA", value: record.lncRnaCount ?? 0},
        {metric: "Pseudogenes", value: record.pseudogeneCount ?? 0},
        {metric: "Small RNA", value: record.smallRnaCount ?? 0},
    ];
    const featureComparison = [
        ["Genes", "geneCount", record.geneCount],
        ["Protein coding", "proteinCoding", record.proteinCoding],
        ["Transcripts", "transcriptCount", record.transcriptCount],
        ["lncRNA", "lncRnaCount", record.lncRnaCount],
        ["Pseudogenes", "pseudogeneCount", record.pseudogeneCount],
        ["Small RNA", "smallRnaCount", record.smallRnaCount],
    ] as const;
    const targetGca = record.gca.toLowerCase();
    const cohortBuscoRecords = cohortRecords.filter((item) => item.genomeBuscoComplete !== null && item.proteinBuscoComplete !== null);
    const cohortProfileRecords = cohortRecords.filter((item) => Number.isFinite(item.pcoa1) && Number.isFinite(item.pcoa2));
    const isTarget = (item: QcRecord) => item.gca.toLowerCase() === targetGca;

    return (
        <div className="mt-8 space-y-6" aria-live="polite">
            <Card>
                <CardHeader>
                    <div className="flex flex-col gap-5 sm:flex-row sm:items-center sm:justify-between">
                        <div className="flex items-center gap-4">
                            <div
                                className="flex h-14 w-14 shrink-0 items-center justify-center rounded-full bg-secondary text-primary">
                                <Dna className="h-7 w-7"/></div>
                            <div>
                                <div className="mb-1 flex flex-wrap items-center gap-2 text-xs text-muted-foreground">
                                    <Badge variant="outline">{record.taxonomy}</Badge><Badge
                                    variant="outline">Taxon {record.taxonId}</Badge><Badge variant="outline"
                                                                                           className="font-mono">{record.gca}</Badge><Badge
                                    variant={getStatusBadge(record.status)}>{record.status}</Badge></div>
                                <CardTitle className="text-2xl">{record.species}</CardTitle><CardDescription>Annotation
                                quality report for one assembly</CardDescription></div>
                        </div>
                    </div>
                </CardHeader>
                <CardContent>
                    <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 lg:grid-cols-4">
                        <MetricCard label="Genome BUSCO complete" value={record.genomeBuscoComplete} suffix="%"/>
                        <MetricCard label="Protein BUSCO complete" value={record.proteinBuscoComplete} suffix="%"/>
                        <MetricCard label="BUSCO delta" value={buscoDelta(record)} suffix="%"/>
                        <MetricCard label="Clade outlier score" value={record.outlierScore} help={OUTLIER_SCORE_HELP}/>
                    </div>
                    <div className="mt-4 flex flex-wrap gap-2 text-sm text-muted-foreground">
                        {reviewReasons(record).map((reason) => <Badge key={reason}
                                                                      variant={needsReview(record) ? "destructive" : "outline"}>{reason}</Badge>)}
                    </div>
                </CardContent>
            </Card>

            <div className="flex gap-1 overflow-x-auto border-b" role="tablist"
                 aria-label="Single annotation report sections">{(["Overview", "Gene structure", "Quality & validation", "All metrics"] as SingleReportTab[]).map((item) =>
                <Button key={item} role="tab" aria-selected={tab === item} variant="ghost"
                        className={`shrink-0 rounded-none border-b-2 px-3 ${tab === item ? "border-primary text-primary" : "border-transparent text-muted-foreground"}`}
                        onClick={() => setTab(item)}>{item}{item === "All metrics" && <span
                    className="ml-2 text-xs text-muted-foreground">{Object.keys(record.newMetrics).length}</span>}</Button>)}</div>

            {tab === "Overview" && <>
                <div className="grid grid-cols-1 gap-6 lg:grid-cols-[1.4fr_1fr]">
                    <Card><CardHeader><CardTitle className="flex items-center gap-2"><Layers3
                        className="h-5 w-5 text-primary"/>Annotation composition</CardTitle><CardDescription>Gene and
                        transcript content from the new metrics
                        table.</CardDescription></CardHeader><CardContent><MetricRows
                        rows={[["Total genes", record.geneCount], ["Protein-coding genes", record.proteinCoding], ["Transcripts", record.transcriptCount], ["lncRNA", record.lncRnaCount], ["Pseudogenes", record.pseudogeneCount], ["Small RNA", record.smallRnaCount]]}/></CardContent></Card>
                    <Card><CardHeader><CardTitle className="flex items-center gap-2"><ShieldCheck
                        className="h-5 w-5 text-primary"/>Quality & completeness</CardTitle><CardDescription>Independent
                        quality signals for this annotation.</CardDescription></CardHeader><CardContent><MetricRows
                        rows={[["Assembly BUSCO", record.genomeBuscoComplete, "%"], ["Annotation BUSCO", record.proteinBuscoComplete, "%"], ["BUSCO delta", buscoDelta(record), "%"], ["Clade outlier score", record.outlierScore]]}/><Button
                        variant="link" className="mt-3 px-0" onClick={() => setTab("Quality & validation")}>View quality
                        details <ArrowRight className="ml-1 h-4 w-4"/></Button></CardContent></Card>
                </div>
                <Card><CardHeader><CardTitle>Assembly context</CardTitle><CardDescription>Assembly measurements used
                    when interpreting annotation quality.</CardDescription></CardHeader><CardContent><MetricRows
                    rows={[["Genome size", record.genomeSizeGb, " Gb"], ["Contig N50", record.contigN50Mb, " Mb"], ["Scaffold N50", record.scaffoldN50Mb, " Mb"], ["GC", record.gcPercent, "%"], ["Repeats", record.repeatPercent, "%"]]}/></CardContent></Card>
            </>}

            {tab === "Gene structure" &&
                <div className="grid grid-cols-1 gap-6 lg:grid-cols-[1fr_1.35fr]"><Card><CardHeader><CardTitle>Annotation
                    detail</CardTitle><CardDescription>Counts and feature composition returned for this
                    annotation.</CardDescription></CardHeader><CardContent><MetricRows
                    rows={[["Genes", record.geneCount], ["Protein coding", record.proteinCoding], ["Transcripts", record.transcriptCount], ["lncRNA", record.lncRnaCount], ["Pseudogenes", record.pseudogeneCount], ["Small RNA", record.smallRnaCount]]}/></CardContent></Card><Card><CardHeader><PlotTitle
                    title="Annotation scale"
                    help="Annotation counts use a normal linear axis. Smaller categories may appear shorter when counts differ greatly."/><CardDescription>Counts
                    from the new metrics table on a linear
                    scale.</CardDescription></CardHeader><CardContent><ChartContainer config={chartConfig}
                                                                                      className="h-80 w-full"><BarChart
                    data={countProfile} layout="vertical" accessibilityLayer><CartesianGrid
                    strokeDasharray="3 3"/><XAxis type="number" tickFormatter={formatCompact}/><YAxis dataKey="metric"
                                                                                                      type="category"
                                                                                                      width={100}
                                                                                                      tick={{fontSize: 11}}/><ChartTooltip
                    content={<ChartTooltipContent/>}/><Bar dataKey="value" fill="var(--chart-2)" name="Count"
                                                           radius={[0, 4, 4, 0]}/></BarChart></ChartContainer></CardContent></Card>
                </div>}

            {tab === "Gene structure" &&
                <Card><CardHeader><CardTitle>Comparison with cohort mean</CardTitle><CardDescription>Feature counts
                    compared with the mean of the relevant QC cohort. The selected annotation is excluded from the mean
                    when peer records are available.</CardDescription></CardHeader><CardContent>
                    <div className="overflow-x-auto">
                        <Table><TableHeader><TableRow><TableHead>Feature</TableHead><TableHead className="text-right">This
                            annotation</TableHead><TableHead className="text-right">Cohort mean</TableHead><TableHead
                            className="text-right">Difference</TableHead></TableRow></TableHeader><TableBody>{featureComparison.map(([label, metric, value]) => {
                            const mean = cohortMeans[metric];
                            const difference = typeof value === "number" && typeof mean === "number" ? value - mean : null;
                            return <TableRow key={metric}><TableCell>{label}</TableCell><TableCell
                                className="text-right">{formatNumber(value)}</TableCell><TableCell
                                className="text-right">{formatNumber(mean)}</TableCell><TableCell
                                className={`text-right ${difference !== null && difference < 0 ? "text-destructive" : ""}`}>{formatNumber(difference)}</TableCell></TableRow>;
                        })}</TableBody></Table></div>
                </CardContent></Card>}

            {tab === "Quality & validation" && <>
                <div className="grid grid-cols-1 gap-6 lg:grid-cols-[1.35fr_1fr]"><Card><CardHeader><PlotTitle
                    title="BUSCO profile"
                    help="Each row is a stacked BUSCO run. Complete BUSCOs are split into single-copy and duplicated; duplicated is not added on top of complete because it is a subset of complete."/><CardDescription>Assembly
                    and annotation BUSCO composition on the same percentage
                    scale.</CardDescription></CardHeader><CardContent><ChartContainer config={chartConfig}
                                                                                      className="h-72 w-full"><BarChart
                    data={buscoProfile} layout="vertical" accessibilityLayer><CartesianGrid
                    strokeDasharray="3 3"/><XAxis type="number" domain={[0, 100]}
                                                  tickFormatter={(value) => `${Math.round(value)}%`}/><YAxis
                    dataKey="metric" type="category" width={110}/><ChartTooltip content={<ChartTooltipContent/>}/><Bar
                    dataKey="singleCopy" stackId="busco" fill="var(--chart-2)" name="Single-copy complete"/><Bar
                    dataKey="duplicated" stackId="busco" fill="var(--chart-4)" name="Duplicated complete"/><Bar
                    dataKey="fragmented" stackId="busco" fill="var(--chart-3)" name="Fragmented"/><Bar dataKey="missing"
                                                                                                       stackId="busco"
                                                                                                       fill="var(--destructive)"
                                                                                                       name="Missing"/></BarChart></ChartContainer></CardContent></Card><Card><CardHeader><CardTitle>Screening
                    interpretation</CardTitle><CardDescription>Compared with live annotations in the same clade. A score
                    of 2 or more prioritises manual review; it is not an automatic
                    failure.</CardDescription></CardHeader><CardContent><MetricRows
                    rows={[["Clade outlier score", record.outlierScore], ["BUSCO delta", buscoDelta(record), "%"], ["BUSCO missing", record.buscoMissing, "%"], ["BUSCO duplicated", record.buscoDuplicated, "%"], ["Review status", needsReview(record) ? "Needs review" : "Within cohort range"]]}/>
                    <div
                        className="mt-4 rounded-md bg-secondary/50 p-3 text-sm text-muted-foreground">{OUTLIER_SCORE_HELP}</div>
                    <div className="mt-4 flex flex-wrap gap-2">{reviewReasons(record).map((reason) => <Badge
                        key={reason} variant={needsReview(record) ? "destructive" : "outline"}>{reason}</Badge>)}</div>
                </CardContent></Card></div>
                <div className="grid grid-cols-1 gap-6 lg:grid-cols-2"><Card>
                    <CardHeader>
                        <PlotTitle title="BUSCO cohort comparison"
                                   help="Each point is one annotation in the comparison cohort. The selected GCA is highlighted. The tooltip reports the completeness delta between assembly and annotation BUSCO."/>
                        <CardDescription>Assembly BUSCO versus annotation BUSCO for the cohort.</CardDescription>
                    </CardHeader><CardContent><ChartContainer config={chartConfig} className="h-80 w-full">
                    <ScatterChart accessibilityLayer>
                        <CartesianGrid strokeDasharray="3 3"/>
                        <XAxis type="number" dataKey="genomeBuscoComplete" name="Assembly BUSCO"
                               domain={["auto", "auto"]} unit="%"/>
                        <YAxis type="number" dataKey="proteinBuscoComplete" name="Annotation BUSCO"
                               domain={["auto", "auto"]} unit="%"/>
                        <ChartTooltip cursor={{strokeDasharray: "3 3"}}
                                      content={<BuscoTooltip threshold={buscoDeltaThreshold(cohortBuscoRecords)}/>}/>
                        <Scatter data={cohortBuscoRecords.filter((item) => !isTarget(item))} fill="var(--chart-2)"
                                 name="Cohort annotations"/>
                        <Scatter data={cohortBuscoRecords.filter(isTarget)} fill="var(--foreground)" name="Selected GCA"/></ScatterChart></ChartContainer><PointLegend
                    flaggedLabel="Selected GCA" flaggedColor="var(--foreground)"/></CardContent></Card><Card><CardHeader><PlotTitle
                    title="Multivariate cohort position"
                    help="This shows the selected GCA against the same multivariate cohort used for screening. The selected GCA is highlighted; the other annotations are shown for context."/><CardDescription>Multivariate
                    profile with the selected GCA
                    highlighted.</CardDescription></CardHeader><CardContent><ChartContainer config={chartConfig}
                                                                                            className="h-80 w-full"><ScatterChart
                    accessibilityLayer><CartesianGrid strokeDasharray="3 3"/><XAxis type="number" dataKey="pcoa1"
                                                                                    name="Profile axis 1"/><YAxis
                    type="number" dataKey="pcoa2" name="Profile axis 2"/><ChartTooltip cursor={{strokeDasharray: "3 3"}}
                                                                                       content={
                                                                                           <MultivariateTooltip/>}/><Scatter
                    data={cohortProfileRecords.filter((item) => !isTarget(item))} fill="var(--chart-2)"
                    name="Cohort annotations"/><Scatter data={cohortProfileRecords.filter(isTarget)}
                                                        fill="var(--foreground)"
                                                        name="Selected GCA"/></ScatterChart></ChartContainer><PointLegend
                    flaggedLabel="Selected GCA" flaggedColor="var(--foreground)"/></CardContent></Card></div>
            </>}

            {tab === "All metrics" && <>
                <Card><CardHeader><CardTitle>Build and provenance</CardTitle><CardDescription>Registry metadata for this
                    annotation.</CardDescription></CardHeader><CardContent>
                    <div className="grid grid-cols-1 gap-x-8 gap-y-1 md:grid-cols-2"><MetricRows
                        rows={[["BioProject", record.associatedProject], ["Group", record.groupName], ["Infrastructure", record.infraName], ["Genebuilder", record.genebuilder], ["Annotation source", record.annotationSource], ["Annotation method", record.annotationMethod]]}/><MetricRows
                        rows={[["Started", record.dateStarted], ["Release date", record.releaseDate], ["Status updated", record.dateStatusUpdate], ["Last genebuild update", record.lastGenebuildUpdate], ["Assembly BUSCO lineage", record.assemblyBuscoLineage], ["Assembly BUSCO version", record.assemblyBuscoVersion], ["Annotation BUSCO lineage", record.proteinBuscoLineage], ["Annotation BUSCO version", record.proteinBuscoVersion]]}/>
                    </div>
                </CardContent></Card>
                <Card><CardHeader><CardTitle>Complete metric detail</CardTitle><CardDescription>Assembly, annotation
                    BUSCO, and every metric returned from the new metrics
                    table.</CardDescription></CardHeader><CardContent>
                    <div className="grid grid-cols-1 gap-x-8 gap-y-1 md:grid-cols-2"><MetricRows
                        rows={[["Assembly BUSCO complete", record.genomeBuscoComplete, "%"], ["Assembly BUSCO duplicated", record.genomeBuscoDuplicated, "%"], ["Assembly BUSCO fragmented", record.genomeBuscoFragmented, "%"], ["Assembly BUSCO missing", record.genomeBuscoMissing, "%"], ["Annotation BUSCO complete", record.proteinBuscoComplete, "%"], ["Annotation BUSCO duplicated", record.buscoDuplicated, "%"], ["Annotation BUSCO fragmented", record.buscoFragmented, "%"], ["Annotation BUSCO missing", record.buscoMissing, "%"], ["Review signals", reviewReasons(record).join(", ")]]}/>
                    </div>
                    <div className="mt-6"><h3 className="mb-2 text-sm font-medium">New metrics</h3><MetricRows
                        rows={Object.entries(record.newMetrics).sort(([left], [right]) => left.localeCompare(right)).map(([name, value]) => [formatMetricLabel(name), formatMetricValue(value)])}/>
                    </div>
                </CardContent></Card>
            </>}
        </div>
    );
}

function MetricCard({label, value, suffix = "", help}: {
    label: string;
    value: number | string | null | undefined;
    suffix?: string;
    help?: string
}) {
    const displayValue = typeof value === "number" ? formatNumber(value) : value ?? "n/a";
    return <Card><CardHeader className="pb-2"><CardDescription className="flex items-center gap-1">{label}{help &&
        <Tooltip><TooltipTrigger asChild><Button type="button" variant="ghost" size="icon" className="h-5 w-5"
                                                 aria-label={`Explain ${label}`}><Info
            className="h-3.5 w-3.5"/></Button></TooltipTrigger><TooltipContent
            className="max-w-xs leading-relaxed">{help}</TooltipContent></Tooltip>}</CardDescription><CardTitle
        className="flex items-baseline gap-1 text-2xl"><Activity
        className="h-4 w-4 text-muted-foreground"/>{displayValue}<span
        className="text-sm font-normal text-muted-foreground">{suffix}</span></CardTitle></CardHeader></Card>;
}

function PlotTitle({title, help}: { title: string; help: string }) {
    return <div className="flex items-center gap-2"><CardTitle>{title}</CardTitle><Tooltip><TooltipTrigger
        asChild><Button type="button" variant="ghost" size="icon" className="h-6 w-6"
                        aria-label={`Explain ${title}`}><Info
        className="h-4 w-4"/></Button></TooltipTrigger><TooltipContent
        className="max-w-xs leading-relaxed">{help}</TooltipContent></Tooltip></div>;
}

function MetricRows({rows}: { rows: Array<[string, number | string | null | undefined, string?]> }) {
    return <div className="divide-y">{rows.map(([label, value, suffix = ""]) => <div key={label}
                                                                                     className="flex items-center justify-between gap-6 py-3">
        <span className="text-sm text-muted-foreground">{label}</span><span
        className="text-right text-sm font-medium">{typeof value === "number" ? `${formatNumber(value)}${suffix}` : value ?? "n/a"}</span>
    </div>)}</div>;
}
