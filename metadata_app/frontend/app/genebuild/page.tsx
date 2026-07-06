"use client";

import React, {useEffect, useMemo, useState} from "react";
import { BackgroundGradient } from "@/components/ui/backround-gradient"

import GridList02 from "@/components/ui/blocks/select_user"
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Checkbox } from "@/components/ui/checkbox";
import { Input } from "@/components/ui/input";
import {
  Card,
  CardAction,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import {
  Dialog,
  DialogClose,
  DialogContent,
  DialogDescription,
  DialogFooter,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
} from "@/components/ui/dialog";
import {
  Select,
  SelectContent,
  SelectGroup,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import {
  Activity,
  AlertTriangle,
  LucideHeartCrack,
  Smile,
  CircleDashed,
  ListChecks,
  LockKeyhole,
  Download,
  RefreshCw,
  ArrowUpDown,
  UserRound,
} from "lucide-react";

type DashboardItem = {
  id: string;
  gca: string;
  scientific_name: string;
  status: string;
  queue: string;
  action: string;
  priority: "high" | "medium" | "normal";
  project?: string;
  method?: string;
  updated?: string;
  daysSinceUpdate?: number | null;
};

type AnnotationOverview = {
  genebuild_status_id: number;
  gca: string | null;
  scientific_name: string | null;
  gb_status: string | null;
  dashboard_status: string | null;
  queue: string | null;
  next_action: string | null;
  priority: "high" | "medium" | "normal";
  bioproject_name?: string | null;
  annotation_method?: string | null;
  date_status_update?: string | null;
  days_since_update?: number | null;
};

type NeedActionApiItem = {
  gca: string | null;
  scientific_name: string | null;
  gb_status: string | null;
  bioproject_name?: string | null;
  date_status_update?: string | null;
  days_since_update?: number | null;
};

type StatusSummaryItem = {
  gb_status: string;
  main: number;
  anno: number;
  hprc: number;
  bioprojects?: string;
};

type OverviewSortKey = "status" | "priority" | "updated";
type SortDirection = "asc" | "desc";

const formatStatus = (status: string) =>
  status
    .replace(/_/g, " ")
    .replace(/\b\w/g, (letter) => letter.toUpperCase());

const stringValue = (value: string | null | undefined, fallback = "Unknown") =>
  value && value.trim() ? value : fallback;

const getPriorityBadge = (priority: DashboardItem["priority"]) => {
  if (priority === "high") return "destructive";
  if (priority === "medium") return "secondary";
  return "outline";
};

const GENEBUILD_PASSWORD = process.env.NEXT_PUBLIC_GENEBUILD_PASSWORD ?? "genebuild";
const GENEBUILD_AUTH_KEY = "genebuildDashboardUnlocked";


export default function Page() {
  const title = "Genebuild dashboard";
  const [selectedGenebuilder, setSelectedGenebuilder] = useState<string | null>(null);
  const [authChecked, setAuthChecked] = useState(false);
  const [isUnlocked, setIsUnlocked] = useState(false);
  const [password, setPassword] = useState("");
  const [passwordError, setPasswordError] = useState<string | null>(null);

  const handleGenebuilderChange = (user: {
    name: string;
    role: string;
    imageUrl: string;
  }) => {
    setSelectedGenebuilder(user.role);
    // Save to browser
  localStorage.setItem("selectedGenebuilder", user.role);
  };
    const [horeadyCount, setCountHOR] = useState<number>(0);
    const [dataCount, setCountData] = useState<number>(0);
    const [pendingCount, setPending] = useState<number>(0);
    const [annotationOverview, setAnnotationOverview] = useState<AnnotationOverview[]>([]);
    const [backendStatusSummary, setBackendStatusSummary] = useState<StatusSummaryItem[]>([]);
    const [staleDataItems, setStaleDataItems] = useState<NeedActionApiItem[]>([]);
    const [stalePendingItems, setStalePendingItems] = useState<NeedActionApiItem[]>([]);

  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [selectedAnnotationIds, setSelectedAnnotationIds] = useState<string[]>([]);
  const [selectedStatus, setSelectedStatus] = useState<string>("");
  const [statusDialogOpen, setStatusDialogOpen] = useState(false);
  const [overviewSort, setOverviewSort] = useState<{
    key: OverviewSortKey;
    direction: SortDirection;
  }>({ key: "priority", direction: "asc" });


  const handleGetHO = async () => {
    if (!selectedGenebuilder) return;
    setLoading(true);
    setError(null);
    try {
      const payload = { genebuilder: selectedGenebuilder }; // string

      const res = await fetch("/api/handover/handover/genebuilder", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(payload),
      });

      if (!res.ok) {
        throw new Error("Failed to load dashboard data");
      }

      const result = await res.json();
      console.log("API response:", result);

      if (result.df_ready) {
        setCountHOR(result.count_ho_ready);
        setCountData(result.count_data);
        setPending(result.count_pending);
        setAnnotationOverview(result.annotation_overview ?? []);
        setBackendStatusSummary(result.status_summary ?? []);
        setStaleDataItems(result.list_data ?? []);
        setStalePendingItems(result.list_pending ?? []);
        setSelectedAnnotationIds([]);


      } else {
        setCountHOR(0);
        setCountData(0);
        setPending(0);
        setAnnotationOverview([]);
        setBackendStatusSummary([]);
        setStaleDataItems([]);
        setStalePendingItems([]);
        setSelectedAnnotationIds([]);
        setError("No dashboard data found for this genebuilder.");
      }

    } catch (error) {
      console.error("Error fetching data:", error);
      setError("Unable to load dashboard data. Please try again.");
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
  const saved = localStorage.getItem("selectedGenebuilder");
  if (saved) {
    setSelectedGenebuilder(saved);
  }

  setIsUnlocked(sessionStorage.getItem(GENEBUILD_AUTH_KEY) === "true");
  setAuthChecked(true);
}, []);

  useEffect(() => {
    if (isUnlocked && selectedGenebuilder) handleGetHO();
  }, [isUnlocked, selectedGenebuilder]);

  const handlePasswordSubmit = (event: React.FormEvent<HTMLFormElement>) => {
    event.preventDefault();

    if (password === GENEBUILD_PASSWORD) {
      sessionStorage.setItem(GENEBUILD_AUTH_KEY, "true");
      setIsUnlocked(true);
      setPassword("");
      setPasswordError(null);
      return;
    }

    setPasswordError("Incorrect password.");
  };

  const dashboardItems = useMemo<DashboardItem[]>(() => {
    if (annotationOverview.length) {
      return annotationOverview.map((item) => ({
        id: `overview-${item.genebuild_status_id}`,
        gca: stringValue(item.gca, ""),
        scientific_name: stringValue(item.scientific_name),
        status: stringValue(item.dashboard_status),
        queue: stringValue(item.queue),
        action: stringValue(item.next_action, "Review current status and next action"),
        priority: item.priority,
        project: stringValue(item.bioproject_name),
        method: item.annotation_method ?? undefined,
        updated: item.date_status_update ?? undefined,
        daysSinceUpdate: item.days_since_update,
      }));
    }

    return [];
  }, [annotationOverview]);

  const needActionItems = useMemo<DashboardItem[]>(() => {
    const dataItems = staleDataItems.map((item, index) => ({
      id: `data-${item.gca}-${item.gb_status}-${index}`,
      gca: stringValue(item.gca, ""),
      scientific_name: stringValue(item.scientific_name),
      status: stringValue(item.gb_status),
      queue: "Data quality",
      action: "Check BUSCO, transcript/protein evidence, or mark abandoned",
      priority: "high" as const,
      project: stringValue(item.bioproject_name),
      updated: item.date_status_update ?? undefined,
      daysSinceUpdate: item.days_since_update,
    }));

    const pendingItems = stalePendingItems.map((item, index) => ({
      id: `pending-${item.gca}-${item.gb_status}-${index}`,
      gca: stringValue(item.gca, ""),
      scientific_name: stringValue(item.scientific_name),
      status: stringValue(item.gb_status),
      queue: "In progress",
      action: "Follow up stale in-progress annotation",
      priority: "high" as const,
      project: stringValue(item.bioproject_name),
      updated: item.date_status_update ?? undefined,
      daysSinceUpdate: item.days_since_update,
    }));

    return [...dataItems, ...pendingItems];
  }, [staleDataItems, stalePendingItems]);
  const selectableDashboardItems = useMemo(
    () => dashboardItems.filter((item) => item.gca),
    [dashboardItems],
  );
  const selectedOverviewItems = useMemo(
    () => dashboardItems.filter((item) => selectedAnnotationIds.includes(item.id) && item.gca),
    [dashboardItems, selectedAnnotationIds],
  );
  const statusChangeItems = useMemo(
    () => selectedOverviewItems.filter((item) => item.method),
    [selectedOverviewItems],
  );
  const sortedDashboardItems = useMemo(() => {
    const priorityRank: Record<DashboardItem["priority"], number> = {
      high: 0,
      medium: 1,
      normal: 2,
    };

    return [...dashboardItems].sort((a, b) => {
      let comparison = 0;

      if (overviewSort.key === "priority") {
        comparison = priorityRank[a.priority] - priorityRank[b.priority];
      } else if (overviewSort.key === "updated") {
        const aTime = a.updated ? new Date(a.updated).getTime() : Number.NEGATIVE_INFINITY;
        const bTime = b.updated ? new Date(b.updated).getTime() : Number.NEGATIVE_INFINITY;
        comparison = aTime - bTime;
      } else {
        comparison = stringValue(a.status).localeCompare(stringValue(b.status));
      }

      if (comparison === 0) {
        comparison = stringValue(a.scientific_name).localeCompare(stringValue(b.scientific_name));
      }

      return overviewSort.direction === "asc" ? comparison : -comparison;
    });
  }, [dashboardItems, overviewSort]);
  const allOverviewRowsSelected =
    selectableDashboardItems.length > 0 &&
    selectableDashboardItems.every((item) => selectedAnnotationIds.includes(item.id));

  const toggleOverviewSort = (key: OverviewSortKey) => {
    setOverviewSort((current) => ({
      key,
      direction: current.key === key && current.direction === "asc" ? "desc" : "asc",
    }));
  };

  const sortLabel = (key: OverviewSortKey) =>
    overviewSort.key === key ? (overviewSort.direction === "asc" ? "ascending" : "descending") : "none";

  const SortableHead = ({
    sortKey,
    children,
  }: {
    sortKey: OverviewSortKey;
    children: React.ReactNode;
  }) => (
    <Button
      type="button"
      variant="ghost"
      size="sm"
      className="-ml-3 h-8 px-3"
      onClick={() => toggleOverviewSort(sortKey)}
      aria-label={`Sort by ${children}`}
      aria-sort={sortLabel(sortKey) as "ascending" | "descending" | "none" | "other"}
    >
      {children}
      <ArrowUpDown className="size-3.5" />
    </Button>
  );

  const toggleAnnotationSelection = (id: string, checked: boolean) => {
    setSelectedAnnotationIds((current) =>
      checked ? Array.from(new Set([...current, id])) : current.filter((itemId) => itemId !== id),
    );
  };

  const toggleAllAnnotationRows = (checked: boolean) => {
    setSelectedAnnotationIds(checked ? selectableDashboardItems.map((item) => item.id) : []);
  };

  const handleChangeStatus = async () => {
    if (!selectedGenebuilder) {
      alert("No genebuilder specified");
      return;
    }

    if (!statusChangeItems.length) {
      alert("Select at least one annotation with an annotation method");
      return;
    }

    if (!selectedStatus) {
      alert("Please select a new status");
      return;
    }

    const items = statusChangeItems.map((item) => ({
      gca: item.gca,
      annotation_method: item.method as string,
    }));

    try {
      const response = await fetch("/api/handover/handover/change_status", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
        },
        body: JSON.stringify({ genebuilder: selectedGenebuilder, items, new_status: selectedStatus }),
      });

      if (!response.ok) {
        throw new Error("Failed to update status");
      }

      alert(`Updated ${items.length} records successfully`);
      setStatusDialogOpen(false);
      setSelectedAnnotationIds([]);
      setSelectedStatus("");
      await handleGetHO();
    } catch (error) {
      console.error(error);
      alert("Error updating records");
    }
  };

  const handleDownloadSelectedGcas = () => {
    const selectedGcas = Array.from(
      new Set(selectedOverviewItems.map((item) => item.gca).filter(Boolean)),
    );

    if (!selectedGcas.length) return;

    const blob = new Blob([`${selectedGcas.join("\n")}\n`], {
      type: "text/plain;charset=utf-8",
    });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    const genebuilder = selectedGenebuilder ?? "genebuilder";
    const safeGenebuilder = genebuilder.replace(/[^a-z0-9_-]+/gi, "_").toLowerCase();

    link.href = url;
    link.download = `${safeGenebuilder}_selected_gcas.txt`;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
  };

  const totalTracked = backendStatusSummary.reduce(
    (total, item) => total + item.main + item.anno + item.hprc,
    0,
  );
  const actionNeeded = dataCount + pendingCount;
  const suggestions = [
    "Add due dates or expected handover week so delayed annotations are visible before they become stale.",
    "Track a blocker reason and next action owner for every non-handover annotation.",
    "Show last status update age and highlight records with no movement in 30, 90, and 180 days.",
    "Add quick filters for project, clade, annotation method, and priority.",
  ];


  return (
    <div className="relative min-h-screen">
    <div
      className={`flex w-full justify-center overflow-x-hidden px-4 py-8 transition duration-200 sm:px-6 lg:px-8 ${
        isUnlocked ? "" : "pointer-events-none select-none blur-sm"
      }`}
      aria-hidden={!isUnlocked}
    >
      <div className="w-full max-w-6xl min-w-0">
        <div className="mb-8 flex min-w-0 flex-col gap-5 lg:flex-row lg:items-end lg:justify-between">
          <div className="max-w-3xl min-w-0">
            <div className="mb-3 flex items-center gap-2 text-sm font-medium text-muted-foreground">
              <Activity className="size-4" />
              Personal annotation tracking
            </div>
            <h1 className="scroll-m-20 text-3xl font-extrabold tracking-tight text-balance sm:text-4xl">{title}</h1>
            <p className="mt-3 text-muted-foreground">
              See assigned annotations, current status queues, and handover-ready cores in one workspace.
            </p>
          </div>
          <div className="flex min-w-0 flex-col gap-2 sm:flex-row sm:items-center">
            {selectedGenebuilder ? (
              <Button variant="outline" onClick={handleGetHO} disabled={loading}>
                <RefreshCw className="size-4" />
                {loading ? "Refreshing" : "Refresh"}
              </Button>
            ) : null}
            <GridList02 value={selectedGenebuilder} onSelect={handleGenebuilderChange} placeholder="Select a Genebuilder" />
          </div>
        </div>

      {selectedGenebuilder ? (
        <div className="grid min-w-0 gap-8">
          {error ? (
            <Card className="border-destructive/40 bg-destructive/5">
              <CardContent className="flex items-start gap-3 pt-6">
                <AlertTriangle className="mt-0.5 size-5 text-destructive" />
                <p className="text-sm">{error}</p>
              </CardContent>
            </Card>
          ) : null}

          <div className="grid min-w-0 gap-4 sm:grid-cols-2 xl:grid-cols-4">
            <BackgroundGradient>
            <Card className="min-w-0">
              <CardHeader>
                <CardDescription>Genebuilder</CardDescription>
                <CardTitle className="flex items-center gap-2 text-2xl">
                  <UserRound className="size-5 text-muted-foreground" />
                  {selectedGenebuilder}
                </CardTitle>
              </CardHeader>
            </Card>
              </BackgroundGradient>
            <BackgroundGradient>
            <Card className="min-w-0">
              <CardHeader>
                <CardDescription>Tracked annotations</CardDescription>
                <CardTitle className="text-2xl">{totalTracked.toLocaleString()}</CardTitle>
              </CardHeader>
            </Card>
              </BackgroundGradient>
            <BackgroundGradient>
            <Card className="min-w-0">
              <CardHeader>
                <CardDescription>Need action</CardDescription>
                <CardTitle className="flex items-center gap-2 text-2xl">
                  {actionNeeded.toLocaleString()}
                  {actionNeeded ? <LucideHeartCrack className="size-5 text-destructive" /> : <Smile className="size-5 text-primary" />}
                </CardTitle>
              </CardHeader>
            </Card>
              </BackgroundGradient>
            <BackgroundGradient>
            <Card className="min-w-0">
              <CardHeader>
                <CardDescription>Ready to hand over</CardDescription>
                <CardTitle className="text-2xl">{horeadyCount.toLocaleString()}</CardTitle>
              </CardHeader>
            </Card>
            </BackgroundGradient>
          </div>

          <Card className="min-w-0 overflow-hidden">
            <CardHeader>
              <CardTitle>Registry status breakdown</CardTitle>
              <CardDescription>Breakdown of annotations into dashboard groups, with unique BioProjects listed for each group.</CardDescription>
            </CardHeader>
            <CardContent className="min-w-0 overflow-x-auto">
              <Table className="min-w-[620px]">
                <TableHeader>
                  <TableRow>
                    <TableHead>Status</TableHead>
                    <TableHead className="text-right">Main</TableHead>
                    <TableHead className="text-right">Anno</TableHead>
                    <TableHead className="text-right">HPRC</TableHead>
                    <TableHead>BioProjects</TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {backendStatusSummary.length ? (
                    backendStatusSummary.map((item) => (
                      <TableRow key={item.gb_status}>
                        <TableCell>
                          <Badge variant="outline">{formatStatus(item.gb_status)}</Badge>
                        </TableCell>
                        <TableCell className="text-right font-semibold">{item.main.toLocaleString()}</TableCell>
                        <TableCell className="text-right font-semibold">{item.anno.toLocaleString()}</TableCell>
                        <TableCell className="text-right font-semibold">{item.hprc.toLocaleString()}</TableCell>
                        <TableCell className="max-w-[520px] whitespace-normal break-words text-muted-foreground">
                          {item.bioprojects ?? "Unknown"}
                        </TableCell>
                      </TableRow>
                    ))
                  ) : (
                    <TableRow>
                      <TableCell colSpan={5} className="h-20 text-center text-muted-foreground">
                        No status summary returned.
                      </TableCell>
                    </TableRow>
                  )}
                </TableBody>
              </Table>
            </CardContent>
          </Card>

          <Card className="min-w-0 overflow-hidden">
            <CardHeader>
              <CardTitle>Need action annotations</CardTitle>
              <CardDescription>Annotations older than 6 months that are in progress or blocked by data quality.</CardDescription>
            </CardHeader>
            <CardContent className="min-w-0 overflow-x-auto">
              <Table className="min-w-[760px]">
                <TableHeader>
                  <TableRow>
                    <TableHead>Annotation</TableHead>
                    <TableHead>Status</TableHead>
                    <TableHead>BioProject</TableHead>
                    <TableHead>Next action</TableHead>
                    <TableHead>Last update</TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {needActionItems.length ? (
                    needActionItems.map((item) => (
                      <TableRow key={`need-action-${item.id}`}>
                        <TableCell className="max-w-[280px] whitespace-normal">
                          <div className="font-medium">{item.scientific_name}</div>
                          <div className="break-all text-xs text-muted-foreground">{item.gca}</div>
                        </TableCell>
                        <TableCell>
                          <Badge variant={item.priority === "high" ? "destructive" : "secondary"}>
                            {formatStatus(item.status)}
                          </Badge>
                        </TableCell>
                        <TableCell className="max-w-[220px] whitespace-normal break-words text-muted-foreground">
                          {item.project ?? "Unknown"}
                        </TableCell>
                        <TableCell className="max-w-[320px] whitespace-normal text-muted-foreground">
                          {item.action}
                        </TableCell>
                        <TableCell className="whitespace-nowrap text-muted-foreground">
                          {item.updated ?? "Unknown"}
                          {typeof item.daysSinceUpdate === "number" ? (
                            <span className="block text-xs">{item.daysSinceUpdate} days ago</span>
                          ) : null}
                        </TableCell>
                      </TableRow>
                    ))
                  ) : (
                    <TableRow>
                      <TableCell colSpan={5} className="h-20 text-center text-muted-foreground">
                        No annotations need action.
                      </TableCell>
                    </TableRow>
                  )}
                </TableBody>
              </Table>
            </CardContent>
          </Card>

          <Card id="annotation-overview" className="min-w-0 overflow-hidden">
            <CardHeader>
              <CardTitle>Annotation overview</CardTitle>
              <CardDescription>All assigned annotations returned by the handover service, ordered by records that need attention first.</CardDescription>
              <CardAction>
                <div className="flex flex-col gap-2 sm:flex-row">
                <Button
                  type="button"
                  variant="outline"
                  disabled={!selectedOverviewItems.length}
                  onClick={handleDownloadSelectedGcas}
                >
                  <Download className="size-4" />
                  Download GCAs
                </Button>
                <Dialog open={statusDialogOpen} onOpenChange={setStatusDialogOpen}>
                  <DialogTrigger asChild>
                    <Button disabled={!statusChangeItems.length}>
                      Change status
                    </Button>
                  </DialogTrigger>
                  <DialogContent className="sm:max-w-[425px]">
                    <DialogHeader>
                      <DialogTitle>Change status in the registry</DialogTitle>
                      <DialogDescription>
                        Select the new registry status for the selected annotations.
                      </DialogDescription>
                    </DialogHeader>
                    <div className="grid gap-4">
                      <Select value={selectedStatus} onValueChange={setSelectedStatus}>
                        <SelectTrigger className="w-full">
                          <SelectValue placeholder="Select status" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectGroup>
                            <SelectItem value="abandoned">Abandoned</SelectItem>
                            <SelectItem value="poor_genome_busco">Low genome BUSCO</SelectItem>
                            <SelectItem value="insufficient_data">Insufficient data</SelectItem>
                            <SelectItem value="check_busco">Low protein BUSCO</SelectItem>
                          </SelectGroup>
                        </SelectContent>
                      </Select>
                      <p className="text-sm text-muted-foreground">
                        {statusChangeItems.length} annotation
                        {statusChangeItems.length !== 1 ? "s" : ""} can be updated.
                      </p>
                      {selectedOverviewItems.length > statusChangeItems.length ? (
                        <p className="text-sm text-muted-foreground">
                          {selectedOverviewItems.length - statusChangeItems.length} selected row
                          {selectedOverviewItems.length - statusChangeItems.length !== 1 ? "s" : ""} can only be used for GCA download.
                        </p>
                      ) : null}
                    </div>
                    <DialogFooter>
                      <DialogClose asChild>
                        <Button variant="outline">Cancel</Button>
                      </DialogClose>
                      <Button type="button" onClick={handleChangeStatus}>
                        Apply
                      </Button>
                    </DialogFooter>
                  </DialogContent>
                </Dialog>
                </div>
              </CardAction>
            </CardHeader>
            <CardContent className="min-w-0 overflow-x-auto">
              <Table className="min-w-[940px]">
                <TableHeader>
                  <TableRow>
                    <TableHead className="w-10">
                      <Checkbox
                        checked={allOverviewRowsSelected || (selectedAnnotationIds.length > 0 && "indeterminate")}
                        onCheckedChange={(checked) => toggleAllAnnotationRows(checked === true)}
                        aria-label="Select all annotations"
                      />
                    </TableHead>
                    <TableHead>Annotation</TableHead>
                    <TableHead>
                      <SortableHead sortKey="status">Status</SortableHead>
                    </TableHead>
                    <TableHead>Queue</TableHead>
                    <TableHead>Next action</TableHead>
                    <TableHead>
                      <SortableHead sortKey="priority">Priority</SortableHead>
                    </TableHead>
                    <TableHead>
                      <SortableHead sortKey="updated">Last update</SortableHead>
                    </TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {dashboardItems.length ? (
                    sortedDashboardItems.map((item) => (
                      <TableRow
                        key={item.id}
                        className={item.gca ? "cursor-pointer select-none" : "cursor-not-allowed opacity-70"}
                        onClick={() => {
                          if (!item.gca) return;
                          toggleAnnotationSelection(item.id, !selectedAnnotationIds.includes(item.id));
                        }}
                        aria-selected={selectedAnnotationIds.includes(item.id)}
                      >
                        <TableCell onClick={(event) => event.stopPropagation()}>
                          <Checkbox
                            checked={selectedAnnotationIds.includes(item.id)}
                            disabled={!item.gca}
                            onCheckedChange={(checked) => toggleAnnotationSelection(item.id, checked === true)}
                            aria-label={`Select ${item.gca}`}
                          />
                        </TableCell>
                        <TableCell className="max-w-[300px] whitespace-normal">
                          <div className="font-medium">{item.scientific_name}</div>
                          <div className="break-all text-xs text-muted-foreground">
                            {item.gca}
                            {item.project ? ` · ${item.project}` : ""}
                            {item.method ? ` · ${item.method}` : ""}
                          </div>
                        </TableCell>
                        <TableCell>
                          <Badge variant={item.status === "handover ready" ? "default" : "secondary"}>
                            {formatStatus(item.status)}
                          </Badge>
                        </TableCell>
                        <TableCell className="whitespace-nowrap">{item.queue}</TableCell>
                        <TableCell className="max-w-[320px] whitespace-normal text-muted-foreground">{item.action}</TableCell>
                        <TableCell>
                          <Badge variant={getPriorityBadge(item.priority)}>
                            {formatStatus(item.priority)}
                          </Badge>
                        </TableCell>
                        <TableCell className="whitespace-nowrap text-muted-foreground">
                          {item.updated ?? "Unknown"}
                          {typeof item.daysSinceUpdate === "number" ? (
                            <span className="block text-xs">{item.daysSinceUpdate} days ago</span>
                          ) : null}
                        </TableCell>
                      </TableRow>
                    ))
                  ) : (
                    <TableRow>
                      <TableCell colSpan={7} className="h-24 text-center text-muted-foreground">
                        {loading ? "Loading annotations..." : "No annotations found for this genebuilder."}
                      </TableCell>
                    </TableRow>
                  )}
                </TableBody>
              </Table>
            </CardContent>
          </Card>

          <Card className="min-w-0 overflow-hidden">
            <CardHeader>
              <CardTitle>Useful additions</CardTitle>
              <CardDescription>Fields that would make this a stronger personalised annotation tracker.</CardDescription>
            </CardHeader>
            <CardContent>
              <div className="grid gap-3 md:grid-cols-2">
                {suggestions.map((suggestion) => (
                  <div key={suggestion} className="flex gap-3 rounded-lg border p-4">
                    <ListChecks className="mt-0.5 size-5 shrink-0 text-primary" />
                    <p className="text-sm text-muted-foreground">{suggestion}</p>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
          </div>
      ) : (
        <Card className="mt-8 border-dashed">
          <CardContent className="flex min-h-64 flex-col items-center justify-center gap-3 text-center">
            <CircleDashed className="size-8 text-muted-foreground" />
            <div>
              <p className="font-medium">Select a genebuilder to open their dashboard</p>
              <p className="text-sm text-muted-foreground">The dashboard will show annotation queues, status summaries, and handover-ready records.</p>
            </div>
          </CardContent>
        </Card>
      )}
        </div>
    </div>
      {!isUnlocked && authChecked ? (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-background/55 px-4 backdrop-blur-[2px]">
          <Card className="w-full max-w-sm shadow-2xl">
            <CardHeader>
              <div className="mb-2 flex size-10 items-center justify-center rounded-md bg-primary/10 text-primary">
                <LockKeyhole className="size-5" />
              </div>
              <CardTitle>Genebuild dashboard</CardTitle>
              <CardDescription>Enter the shared password to continue.</CardDescription>
            </CardHeader>
            <CardContent>
              <form className="grid gap-4" onSubmit={handlePasswordSubmit}>
                <div className="grid gap-2">
                  <Input
                    autoFocus
                    type="password"
                    value={password}
                    onChange={(event) => {
                      setPassword(event.target.value);
                      setPasswordError(null);
                    }}
                    placeholder="Password"
                    aria-label="Password"
                    aria-invalid={Boolean(passwordError)}
                  />
                  {passwordError ? (
                    <p className="text-sm text-destructive">{passwordError}</p>
                  ) : null}
                </div>
                <Button type="submit" className="w-full">
                  Unlock dashboard
                </Button>
              </form>
            </CardContent>
          </Card>
        </div>
      ) : null}
    </div>
  );
}
