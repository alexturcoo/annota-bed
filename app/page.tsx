"use client";

import React, { FormEvent, useMemo, useState } from "react";
import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from "recharts";

type Row = {
  region_id: number;
  input_chr: string;
  input_start: number;
  input_end: number;
  region_size: number;
  gene: string | null;
  strand: string | null;
  feature_biotype: string | null;
  ensembl_id: string | null;
  hugo: string | null;
  tx_overlap_pct: number;
  exon_overlap_pct: number;
  cds_overlap_pct: number;
  priority_score: number;
  ccre_overlap_count: number;
  ccre_best_type: string | null;
  ccre_types: string | null;
  ccre_ids: string | null;
  ccre_max_overlap_bp: number;
};

type Summary = {
  total_regions: number;
  reported_rows: number;
  annotated_regions: number;
  unannotated_regions: number;
  annotated_region_pct: number;
  unique_genes: number;
  median_region_bp: number;
  mean_region_bp: number;
  mean_tx_overlap_pct: number;
  mean_exon_overlap_pct: number;
  mean_cds_overlap_pct: number;
  by_biotype: Record<string, number>;
  by_biotype_regions?: Record<string, number>;
  by_chrom: Record<string, number>;
  region_classes: Record<string, number>;
};

type CcreEnrichment = {
  type: string;
  observed_count: number;
  region_count: number;
  observed_fraction: number;
  region_fraction: number;
  catalog_fraction: number;
  fold_enrichment: number | null;
  region_fold: number | null;
};

type CcreSummary = {
  mode: "off" | "unannotated" | "all";
  enabled: boolean;
  available: boolean;
  reason: string | null;
  tested_regions: number;
  overlap_regions: number;
  overlap_region_pct: number;
  input_overlap_region_pct: number;
  total_overlaps: number;
  by_type: Record<string, number>;
  by_type_regions: Record<string, number>;
  permutation: {
    enabled: boolean;
    method: string;
    permutations: number;
    seed: number;
    tested_regions: number;
    results: Array<{
      type: string;
      observed: number;
      expected: number;
      fold_enrichment: number;
      empirical_p: number;
    }>;
  } | null;
  enrichment: CcreEnrichment[];
};

type ReferencePayload = {
  genome: string;
  label: string;
  source: string;
};

const REFERENCES = [
  {
    key: "hg38",
    label: "GRCh38 / hg38",
    detail: "genes + ENCODE cCREs",
  },
  {
    key: "hg37",
    label: "GRCh37 / hg19",
    detail: "gene annotations",
  },
  {
    key: "chm13",
    label: "T2T-CHM13v2.0",
    detail: "gene annotations",
  },
];

const CHROM_ORDER = [
  "1",
  "2",
  "3",
  "4",
  "5",
  "6",
  "7",
  "8",
  "9",
  "10",
  "11",
  "12",
  "13",
  "14",
  "15",
  "16",
  "17",
  "18",
  "19",
  "20",
  "21",
  "22",
  "X",
  "Y",
  "MT",
  "M",
];

const CLASS_ORDER = [
  "CDS overlap",
  "exonic non-CDS",
  "transcript/intronic",
  "cCRE only",
  "unannotated",
];

const PALETTE = [
  "#2563eb",
  "#0f766e",
  "#dc2626",
  "#ca8a04",
  "#7c3aed",
  "#0891b2",
  "#16a34a",
  "#ea580c",
  "#4f46e5",
  "#64748b",
];

const numberFormatter = new Intl.NumberFormat("en-US");

function fmtNumber(value: number | null | undefined) {
  return numberFormatter.format(Math.round(Number(value || 0)));
}

function fmtPct(value: number | null | undefined) {
  return `${Number(value || 0).toFixed(1)}%`;
}

function fmtDecimal(value: number | null | undefined) {
  return Number(value || 0).toFixed(2);
}

function chrSortKey(chr: string) {
  const normalized = chr.replace(/^chr/i, "");
  const idx = CHROM_ORDER.indexOf(normalized);
  return idx === -1 ? 999 : idx;
}

function functionalClass(row: Row) {
  if (Number(row.cds_overlap_pct || 0) > 0) return "CDS overlap";
  if (Number(row.exon_overlap_pct || 0) > 0) return "exonic non-CDS";
  if (Number(row.tx_overlap_pct || 0) > 0) return "transcript/intronic";
  if (Number(row.ccre_overlap_count || 0) > 0) return "cCRE only";
  return "unannotated";
}

function ChartPanel({
  title,
  metric,
  children,
}: {
  title: string;
  metric?: string;
  children: React.ReactNode;
}) {
  return (
    <div className="rounded-lg border border-slate-200 bg-white p-4 shadow-sm">
      <div className="mb-3 flex items-baseline justify-between gap-3">
        <h3 className="text-sm font-semibold text-slate-900">{title}</h3>
        {metric && <span className="text-xs font-medium text-slate-500">{metric}</span>}
      </div>
      <div className="h-72 min-w-0">{children}</div>
    </div>
  );
}

function EmptyChart({ label }: { label: string }) {
  return (
    <div className="flex h-full items-center justify-center rounded-md border border-dashed border-slate-200 text-sm text-slate-400">
      {label}
    </div>
  );
}

function MetricCard({
  label,
  value,
  accent,
}: {
  label: string;
  value: string;
  accent?: string;
}) {
  return (
    <div className="rounded-lg border border-slate-200 bg-white p-4 shadow-sm">
      <div className="text-xs font-medium uppercase text-slate-500">{label}</div>
      <div className="mt-2 flex items-end gap-2">
        <div className="text-2xl font-semibold text-slate-950">{value}</div>
        {accent && <div className="pb-1 text-xs font-medium text-slate-500">{accent}</div>}
      </div>
    </div>
  );
}

export default function Home() {
  const [bedFile, setBedFile] = useState<File | null>(null);
  const [gtfBgz, setGtfBgz] = useState<File | null>(null);
  const [gtfTbi, setGtfTbi] = useState<File | null>(null);
  const [genome, setGenome] = useState("hg38");
  const [customReference, setCustomReference] = useState(false);
  const [ccreMode, setCcreMode] = useState<"off" | "unannotated" | "all">("off");
  const [runCcrePermutation, setRunCcrePermutation] = useState(true);
  const [ccrePermutations, setCcrePermutations] = useState(1000);
  const [ambiguities, setAmbiguities] = useState("best_all");
  const [query, setQuery] = useState("");
  const [onlyCcre, setOnlyCcre] = useState(false);
  const [loading, setLoading] = useState(false);
  const [rows, setRows] = useState<Row[]>([]);
  const [summary, setSummary] = useState<Summary | null>(null);
  const [ccre, setCcre] = useState<CcreSummary | null>(null);
  const [reference, setReference] = useState<ReferencePayload | null>(null);
  const [csvPath, setCsvPath] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const ccreSelectable = genome === "hg38";

  const handleSubmit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();
    if (!bedFile) {
      setError("Choose a BED file first.");
      return;
    }
    if (customReference && (!gtfBgz || !gtfTbi)) {
      setError("Custom references need both the bgzipped GTF and matching .tbi index.");
      return;
    }

    setError(null);
    setLoading(true);
    setRows([]);
    setSummary(null);
    setCcre(null);
    setReference(null);
    setCsvPath(null);

    const form = new FormData();
    form.append("file", bedFile);
    form.append("genome", genome);
    form.append("reference_mode", customReference ? "custom" : "builtin");
    form.append("ccre_mode", ccreSelectable ? ccreMode : "off");
    form.append("ccre_permutation", ccreSelectable && ccreMode !== "off" && runCcrePermutation ? "true" : "false");
    form.append("ccre_permutations", String(ccrePermutations));
    form.append("ambiguities", ambiguities);
    if (customReference && gtfBgz && gtfTbi) {
      form.append("gtf_bgz", gtfBgz);
      form.append("gtf_tbi", gtfTbi);
    }

    try {
      const response = await fetch("/api/annotate", { method: "POST", body: form });
      const contentType = response.headers.get("content-type") || "";
      const data = contentType.includes("application/json")
        ? await response.json()
        : { error: await response.text() };
      if (!response.ok) {
        throw new Error(data?.error || "Annotation failed");
      }
      setRows(data.rows || []);
      setSummary(data.summary || null);
      setCcre(data.ccre || null);
      setReference(data.reference || null);
      setCsvPath(data.csv_download_path || null);
      setQuery("");
      setOnlyCcre(false);
    } catch (err: any) {
      setError(err.message || "Request failed");
    } finally {
      setLoading(false);
    }
  };

  const chromChartData = useMemo(() => {
    if (!summary) return [];
    return Object.entries(summary.by_chrom)
      .map(([chrom, count]) => ({ chrom, count }))
      .sort((a, b) => chrSortKey(a.chrom) - chrSortKey(b.chrom) || a.chrom.localeCompare(b.chrom));
  }, [summary]);

  const classChartData = useMemo(() => {
    if (!summary) return [];
    const known = CLASS_ORDER.map((name) => ({ name, value: summary.region_classes[name] || 0 }));
    const extra = Object.entries(summary.region_classes)
      .filter(([name]) => !CLASS_ORDER.includes(name))
      .map(([name, value]) => ({ name, value }));
    return [...known, ...extra].filter((item) => item.value > 0);
  }, [summary]);

  const ccreChartData = useMemo(() => {
    if (!ccre) return [];
    if (ccre.permutation?.results?.length) {
      return ccre.permutation.results
        .filter((item) => item.observed > 0 || item.expected > 0)
        .slice(0, 8)
        .map((item) => ({
          name: `${item.type}${item.empirical_p < 0.05 ? " *" : ""}`,
          fold: item.fold_enrichment,
          observed: item.observed,
          expected: item.expected,
          p: item.empirical_p,
        }));
    }
    return [...ccre.enrichment]
      .filter((item) => item.region_count > 0)
      .sort((a, b) => (b.region_fold || 0) - (a.region_fold || 0))
      .slice(0, 8)
      .map((item) => ({
        name: item.type,
        fold: item.region_fold || 0,
        recordFold: item.fold_enrichment || 0,
        count: item.observed_count,
        regions: item.region_count,
        p: null,
      }));
  }, [ccre]);

  const geneNames = useMemo(() => {
    return Array.from(new Set(rows.map((row) => row.hugo || row.gene || "").filter(Boolean))).slice(0, 28);
  }, [rows]);

  const sortedRows = useMemo(() => {
    return [...rows].sort((a, b) => {
      const chrA = a.input_chr.replace(/^chr/i, "");
      const chrB = b.input_chr.replace(/^chr/i, "");
      const chrDelta = chrSortKey(chrA) - chrSortKey(chrB);
      if (chrDelta !== 0) return chrDelta;
      if (chrA !== chrB) return chrA.localeCompare(chrB);
      if (a.input_start !== b.input_start) return a.input_start - b.input_start;
      return (b.priority_score || 0) - (a.priority_score || 0);
    });
  }, [rows]);

  const filteredRows = useMemo(() => {
    const q = query.trim().toLowerCase();
    return sortedRows.filter((row) => {
      if (onlyCcre && !row.ccre_overlap_count) return false;
      if (!q) return true;
      const haystack = [
        row.input_chr,
        row.input_start,
        row.input_end,
        row.hugo,
        row.gene,
        row.ensembl_id,
        row.feature_biotype,
        row.ccre_types,
        row.ccre_ids,
      ]
        .filter(Boolean)
        .join(" ")
        .toLowerCase();
      return haystack.includes(q);
    });
  }, [onlyCcre, query, sortedRows]);

  const displayedRows = filteredRows.slice(0, 500);

  return (
    <main className="min-h-screen bg-[#f6f8fb] text-slate-950">
      <div className="mx-auto max-w-7xl px-4 py-8 sm:px-6 lg:px-8">
        <header className="mb-6 flex flex-col gap-4 lg:flex-row lg:items-end lg:justify-between">
          <div>
            <h1 className="text-3xl font-semibold text-slate-950">annota-bed</h1>
            <p className="mt-2 max-w-2xl text-sm leading-6 text-slate-600">
              BED annotation against installed gene references, with optional GRCh38 cCRE context.
            </p>
          </div>
          {reference && (
            <div className="rounded-lg border border-slate-200 bg-white px-4 py-3 text-sm shadow-sm">
              <span className="font-medium text-slate-900">{reference.label}</span>
              <span className="ml-2 text-slate-500">{reference.source}</span>
            </div>
          )}
        </header>

        <form onSubmit={handleSubmit} className="rounded-lg border border-slate-200 bg-white p-5 shadow-sm">
          <div className="grid gap-5 lg:grid-cols-[1fr_1.2fr]">
            <div>
              <label className="block text-sm font-semibold text-slate-900">Input BED</label>
              <input
                type="file"
                accept=".bed,.tsv,.csv,.txt"
                onChange={(event) => setBedFile(event.target.files?.[0] || null)}
                className="mt-2 w-full cursor-pointer rounded-lg border border-slate-200 bg-slate-50 p-2 text-sm text-slate-700 file:mr-3 file:rounded-md file:border-0 file:bg-slate-900 file:px-3 file:py-2 file:text-sm file:font-medium file:text-white"
                required
              />
              <div className="mt-4 grid gap-3 sm:grid-cols-2">
                <label className="block">
                  <span className="text-sm font-semibold text-slate-900">Ambiguities</span>
                  <select
                    value={ambiguities}
                    onChange={(event) => setAmbiguities(event.target.value)}
                    className="mt-2 w-full rounded-lg border border-slate-200 bg-white px-3 py-2 text-sm text-slate-900 shadow-sm"
                  >
                    <option value="best_all">Best tied matches</option>
                    <option value="best_one">Single best match</option>
                    <option value="all">All overlaps</option>
                  </select>
                </label>
                <label className="block">
                  <span className={ccreSelectable ? "text-sm font-semibold text-slate-900" : "text-sm font-semibold text-slate-400"}>
                    cCRE intersect
                  </span>
                  <select
                    value={ccreSelectable ? ccreMode : "off"}
                    disabled={!ccreSelectable}
                    onChange={(event) => setCcreMode(event.target.value as "off" | "unannotated" | "all")}
                    className="mt-2 w-full rounded-lg border border-slate-200 bg-white px-3 py-2 text-sm text-slate-900 shadow-sm disabled:bg-slate-100 disabled:text-slate-400"
                  >
                    <option value="off">Off</option>
                    <option value="unannotated">Unannotated regions</option>
                    <option value="all">All regions</option>
                  </select>
                </label>
              </div>
              {ccreSelectable && ccreMode !== "off" && (
                <div className="mt-3 grid gap-3 rounded-lg border border-slate-200 bg-slate-50 p-3 sm:grid-cols-[1fr_9rem]">
                  <label className="inline-flex min-h-[38px] items-center gap-2 text-sm font-medium text-slate-700">
                    <input
                      type="checkbox"
                      checked={runCcrePermutation}
                      onChange={(event) => setRunCcrePermutation(event.target.checked)}
                      className="h-4 w-4 accent-slate-900"
                    />
                    Permutation enrichment
                  </label>
                  <label className="block">
                    <span className="sr-only">Permutation count</span>
                    <input
                      type="number"
                      min={100}
                      max={5000}
                      step={100}
                      value={ccrePermutations}
                      disabled={!runCcrePermutation}
                      onChange={(event) => setCcrePermutations(Number(event.target.value))}
                      className="h-10 w-full rounded-lg border border-slate-200 bg-white px-3 text-sm text-slate-900 shadow-sm disabled:bg-slate-100 disabled:text-slate-400"
                    />
                  </label>
                </div>
              )}
            </div>

            <div>
              <div className="flex items-center justify-between gap-3">
                <label className="block text-sm font-semibold text-slate-900">Installed reference</label>
                <label className="inline-flex items-center gap-2 text-sm font-medium text-slate-600">
                  <input
                    type="checkbox"
                    checked={customReference}
                    onChange={(event) => setCustomReference(event.target.checked)}
                    className="h-4 w-4 accent-slate-900"
                  />
                  Custom GTF
                </label>
              </div>
              <div className="mt-2 grid gap-3 md:grid-cols-3">
                {REFERENCES.map((ref) => {
                  const selected = genome === ref.key;
                  return (
                    <button
                      type="button"
                      key={ref.key}
                      onClick={() => setGenome(ref.key)}
                      className={`rounded-lg border p-3 text-left transition ${
                        selected
                          ? "border-slate-950 bg-slate-950 text-white"
                          : "border-slate-200 bg-slate-50 text-slate-800 hover:border-slate-400"
                      }`}
                    >
                      <div className="text-sm font-semibold">{ref.label}</div>
                      <div className={`mt-1 text-xs ${selected ? "text-slate-300" : "text-slate-500"}`}>
                        {ref.detail}
                      </div>
                    </button>
                  );
                })}
              </div>

              {customReference && (
                <div className="mt-4 grid gap-3 md:grid-cols-2">
                  <label className="block">
                    <span className="text-sm font-semibold text-slate-900">GTF .bgz</span>
                    <input
                      type="file"
                      accept=".bgz,.gz"
                      onChange={(event) => setGtfBgz(event.target.files?.[0] || null)}
                      className="mt-2 w-full cursor-pointer rounded-lg border border-slate-200 bg-slate-50 p-2 text-sm text-slate-700 file:mr-3 file:rounded-md file:border-0 file:bg-slate-200 file:px-3 file:py-2 file:text-sm file:font-medium file:text-slate-900"
                    />
                  </label>
                  <label className="block">
                    <span className="text-sm font-semibold text-slate-900">Tabix .tbi</span>
                    <input
                      type="file"
                      accept=".tbi"
                      onChange={(event) => setGtfTbi(event.target.files?.[0] || null)}
                      className="mt-2 w-full cursor-pointer rounded-lg border border-slate-200 bg-slate-50 p-2 text-sm text-slate-700 file:mr-3 file:rounded-md file:border-0 file:bg-slate-200 file:px-3 file:py-2 file:text-sm file:font-medium file:text-slate-900"
                    />
                  </label>
                </div>
              )}
            </div>
          </div>

          <div className="mt-5 flex flex-col gap-3 border-t border-slate-100 pt-5 sm:flex-row sm:items-center sm:justify-between">
            <button
              type="submit"
              disabled={loading || !bedFile || (customReference && (!gtfBgz || !gtfTbi))}
              className="inline-flex h-10 items-center justify-center rounded-lg bg-slate-950 px-5 text-sm font-semibold text-white shadow-sm transition hover:bg-slate-800 disabled:cursor-not-allowed disabled:opacity-45"
            >
              {loading ? "Annotating..." : "Annotate BED"}
            </button>
            {error && <p className="text-sm font-medium text-red-600">{error}</p>}
          </div>
        </form>

        {summary && (
          <section className="mt-6 space-y-6">
            <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-4">
              <MetricCard label="Input regions" value={fmtNumber(summary.total_regions)} />
              <MetricCard
                label="Gene annotated"
                value={fmtPct(summary.annotated_region_pct)}
                accent={`${fmtNumber(summary.annotated_regions)} regions`}
              />
              <MetricCard label="Unique genes" value={fmtNumber(summary.unique_genes)} />
              <MetricCard
                label="cCRE overlap"
                value={ccre?.enabled ? fmtPct(ccre.overlap_region_pct) : "Off"}
                accent={ccre?.enabled ? `${fmtNumber(ccre.overlap_regions)} of ${fmtNumber(ccre.tested_regions)} tested` : ccre?.reason || undefined}
              />
            </div>

            <div className="grid gap-4 lg:grid-cols-3">
              <MetricCard label="Median region bp" value={fmtNumber(summary.median_region_bp)} />
              <MetricCard label="Mean transcript overlap" value={fmtPct(summary.mean_tx_overlap_pct)} />
              <MetricCard label="Reported rows" value={fmtNumber(summary.reported_rows)} />
            </div>

            <div className="grid gap-4 lg:grid-cols-2">
              <ChartPanel title="Functional class">
                {classChartData.length ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart
                      layout="vertical"
                      data={classChartData}
                      margin={{ top: 4, right: 20, left: 8, bottom: 8 }}
                    >
                      <CartesianGrid stroke="#e2e8f0" strokeDasharray="3 3" horizontal={false} />
                      <XAxis type="number" allowDecimals={false} tick={{ fill: "#64748b", fontSize: 12 }} />
                      <YAxis
                        type="category"
                        dataKey="name"
                        width={132}
                        tick={{ fill: "#334155", fontSize: 12 }}
                      />
                      <Tooltip
                        cursor={{ fill: "#f1f5f9" }}
                        contentStyle={{ borderRadius: 8, borderColor: "#e2e8f0", color: "#0f172a" }}
                      />
                      <Bar dataKey="value" radius={[0, 4, 4, 0]}>
                        {classChartData.map((_, index) => (
                          <Cell key={index} fill={PALETTE[index % PALETTE.length]} />
                        ))}
                      </Bar>
                    </BarChart>
                  </ResponsiveContainer>
                ) : (
                  <EmptyChart label="No class data" />
                )}
              </ChartPanel>

              <ChartPanel title="Chromosomes">
                {chromChartData.length ? (
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={chromChartData} margin={{ top: 8, right: 12, left: 0, bottom: 8 }}>
                      <CartesianGrid stroke="#e2e8f0" strokeDasharray="3 3" vertical={false} />
                      <XAxis dataKey="chrom" interval={0} tick={{ fill: "#64748b", fontSize: 11 }} />
                      <YAxis allowDecimals={false} width={42} tick={{ fill: "#64748b", fontSize: 12 }} />
                      <Tooltip
                        cursor={{ fill: "#f1f5f9" }}
                        contentStyle={{ borderRadius: 8, borderColor: "#e2e8f0", color: "#0f172a" }}
                      />
                      <Bar dataKey="count" fill="#0f766e" radius={[4, 4, 0, 0]} />
                    </BarChart>
                  </ResponsiveContainer>
                ) : (
                  <EmptyChart label="No chromosome data" />
                )}
              </ChartPanel>

              {ccre?.enabled && (
                <ChartPanel
                  title={ccre.permutation ? "cCRE permutation enrichment" : "cCRE class signal"}
                  metric={
                    ccre.permutation
                      ? `${fmtNumber(ccre.permutation.permutations)} permutations`
                      : `${fmtNumber(ccre.overlap_regions)} regions / ${fmtNumber(ccre.total_overlaps)} records`
                  }
                >
                  {ccreChartData.length ? (
                    <ResponsiveContainer width="100%" height="100%">
                      <BarChart layout="vertical" data={ccreChartData} margin={{ top: 4, right: 20, left: 8, bottom: 8 }}>
                        <CartesianGrid stroke="#e2e8f0" strokeDasharray="3 3" horizontal={false} />
                        <XAxis type="number" tick={{ fill: "#64748b", fontSize: 12 }} />
                        <YAxis
                          type="category"
                          dataKey="name"
                          width={118}
                          tick={{ fill: "#334155", fontSize: 12 }}
                        />
                        <Tooltip
                          cursor={{ fill: "#f1f5f9" }}
                          formatter={(value, name) => {
                            if (name === "fold") {
                              return [
                                `${fmtDecimal(Number(value))}x`,
                                ccre.permutation ? "observed / expected" : "region fold vs catalog",
                              ];
                            }
                            if (name === "recordFold") return [`${fmtDecimal(Number(value))}x`, "record fold vs catalog"];
                            if (name === "expected") return [fmtDecimal(Number(value)), "mean expected"];
                            if (name === "p") return [fmtDecimal(Number(value)), "empirical p"];
                            return [value, name];
                          }}
                          contentStyle={{ borderRadius: 8, borderColor: "#e2e8f0", color: "#0f172a" }}
                        />
                        <Bar dataKey="fold" radius={[0, 4, 4, 0]}>
                          {ccreChartData.map((_, index) => (
                            <Cell key={index} fill={PALETTE[(index + 3) % PALETTE.length]} />
                          ))}
                        </Bar>
                      </BarChart>
                    </ResponsiveContainer>
                  ) : (
                    <EmptyChart label={ccre.reason || "No cCRE overlaps"} />
                  )}
                </ChartPanel>
              )}
            </div>

            {geneNames.length > 0 && (
              <div className="rounded-lg border border-slate-200 bg-white p-4 shadow-sm">
                <div className="mb-3 flex items-center justify-between gap-3">
                  <h3 className="text-sm font-semibold text-slate-900">Genes</h3>
                  {csvPath && (
                    <a
                      href={`/api/download?path=${encodeURIComponent(csvPath)}`}
                      className="rounded-lg bg-slate-950 px-3 py-2 text-sm font-semibold text-white hover:bg-slate-800"
                    >
                      Download CSV
                    </a>
                  )}
                </div>
                <div className="flex flex-wrap gap-2">
                  {geneNames.map((gene) => (
                    <span
                      key={gene}
                      className="rounded-md border border-slate-200 bg-slate-50 px-2.5 py-1 text-xs font-medium text-slate-700"
                    >
                      {gene}
                    </span>
                  ))}
                </div>
              </div>
            )}
          </section>
        )}

        {sortedRows.length > 0 && (
          <section className="mt-6 rounded-lg border border-slate-200 bg-white shadow-sm">
            <div className="flex flex-col gap-3 border-b border-slate-200 p-4 lg:flex-row lg:items-center lg:justify-between">
              <div>
                <h2 className="text-lg font-semibold text-slate-950">Results</h2>
                <p className="mt-1 text-sm text-slate-500">
                  {fmtNumber(filteredRows.length)} rows{filteredRows.length > displayedRows.length ? " / first 500 shown" : ""}
                </p>
              </div>
              <div className="flex flex-col gap-2 sm:flex-row sm:items-center">
                <input
                  value={query}
                  onChange={(event) => setQuery(event.target.value)}
                  placeholder="Search gene, transcript, cCRE..."
                  className="h-10 w-full rounded-lg border border-slate-200 px-3 text-sm text-slate-900 shadow-sm sm:w-72"
                />
                <label className="inline-flex h-10 items-center gap-2 rounded-lg border border-slate-200 px-3 text-sm font-medium text-slate-700">
                  <input
                    type="checkbox"
                    checked={onlyCcre}
                    onChange={(event) => setOnlyCcre(event.target.checked)}
                    className="h-4 w-4 accent-slate-900"
                  />
                  cCRE only
                </label>
              </div>
            </div>

            <div className="max-h-[560px] overflow-auto">
              <table className="min-w-[1180px] w-full text-left text-sm">
                <thead className="sticky top-0 z-10 bg-slate-100 text-xs uppercase text-slate-500">
                  <tr>
                    {[
                      "Region",
                      "Gene",
                      "Functional class",
                      "Strand",
                      "Tx %",
                      "Exon %",
                      "CDS %",
                      "Score",
                      "cCRE",
                      "cCRE IDs",
                      "Transcript",
                    ].map((header) => (
                      <th key={header} className="whitespace-nowrap px-3 py-3 font-semibold">
                        {header}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody className="divide-y divide-slate-100">
                  {displayedRows.map((row, index) => (
                    <tr key={`${row.region_id}-${row.ensembl_id || "none"}-${index}`} className="hover:bg-slate-50">
                      <td className="whitespace-nowrap px-3 py-3 font-mono text-xs text-slate-700">
                        chr{row.input_chr}:{fmtNumber(row.input_start)}-{fmtNumber(row.input_end)}
                      </td>
                      <td className="px-3 py-3 font-medium text-slate-900">{row.hugo || row.gene || "-"}</td>
                      <td className="px-3 py-3 text-slate-700">{functionalClass(row)}</td>
                      <td className="px-3 py-3 text-slate-700">{row.strand || "-"}</td>
                      <td className="px-3 py-3 text-slate-700">{fmtDecimal(row.tx_overlap_pct)}</td>
                      <td className="px-3 py-3 text-slate-700">{fmtDecimal(row.exon_overlap_pct)}</td>
                      <td className="px-3 py-3 text-slate-700">{fmtDecimal(row.cds_overlap_pct)}</td>
                      <td className="px-3 py-3 text-slate-700">{fmtDecimal(row.priority_score)}</td>
                      <td className="px-3 py-3">
                        {row.ccre_overlap_count ? (
                          <span className="rounded-md bg-emerald-50 px-2 py-1 text-xs font-semibold text-emerald-700">
                            {row.ccre_types || row.ccre_best_type} ({row.ccre_overlap_count})
                          </span>
                        ) : (
                          <span className="text-slate-400">-</span>
                        )}
                      </td>
                      <td className="max-w-[220px] truncate px-3 py-3 text-slate-700" title={row.ccre_ids || ""}>
                        {row.ccre_ids || "-"}
                      </td>
                      <td className="max-w-[220px] truncate px-3 py-3 font-mono text-xs text-slate-600" title={row.ensembl_id || ""}>
                        {row.ensembl_id || "-"}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </section>
        )}
      </div>
    </main>
  );
}
