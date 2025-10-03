"use client";

import React, { useState } from "react";
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, CartesianGrid, Cell, Label,
} from "recharts";

type Row = {
  input_chr: string;
  input_start: number;
  input_end: number;
  gene: string | null;
  strand: string | null;
  feature_biotype: string | null;
  ensembl_id: string | null;
  hugo: string | null;
  tx_overlap_pct: number;
  exon_overlap_pct: number;
  cds_overlap_pct: number;
  priority_score: number;
};

type Summary = {
  unique_genes: number;
  by_biotype: Record<string, number>;
};

export default function Home() {
  const [bedFile, setBedFile] = useState<File | null>(null);
  const [gtfBgz, setGtfBgz] = useState<File | null>(null);
  const [gtfTbi, setGtfTbi] = useState<File | null>(null);
  const [ambiguities, setAmbiguities] = useState("best_all");
  const [loading, setLoading] = useState(false);
  const [rows, setRows] = useState<Row[]>([]);
  const [summary, setSummary] = useState<Summary | null>(null);
  const [csvPath, setCsvPath] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!bedFile) return setError("Please choose a BED file first.");
    if (!gtfBgz || !gtfTbi) return setError("Please upload both GTF (.gtf.bgz) and its .tbi index.");
    setError(null);
    setLoading(true);
    setRows([]); setSummary(null); setCsvPath(null);

    const fd = new FormData();
    fd.append("file", bedFile);
    fd.append("gtf_bgz", gtfBgz);
    fd.append("gtf_tbi", gtfTbi);
    fd.append("ambiguities", ambiguities);

    try {
      const res = await fetch("/api/annotate", { method: "POST", body: fd });
      const ct = res.headers.get("content-type") || "";
      const data = ct.includes("application/json") ? await res.json() : { error: await res.text() };
      if (!res.ok) throw new Error(data?.error || "Annotation failed");
      setRows(data.rows || []);
      setSummary(data.summary || null);
      setCsvPath(data.csv_download_path || null);
    } catch (err: any) {
      setError(err.message || "Request failed");
    } finally {
      setLoading(false);
    }
  };

  // ===== Charts data =====
  const biotypeChartData = summary
    ? Object.entries(summary.by_biotype).map(([biotype, count]) => ({ biotype, count }))
    : [];

  const CHROM_ORDER = ["1","2","3","4","5","6","7","8","9","10","11","12",
                       "13","14","15","16","17","18","19","20","21","22","X","Y"];
  const chrCounts = rows.reduce<Record<string, number>>((acc, r) => {
    const c = (r.input_chr || "").replace(/^chr/i, "");
    acc[c] = (acc[c] || 0) + 1;
    return acc;
  }, {});
  const chromChartData = CHROM_ORDER
    .filter((c) => chrCounts[c] > 0)
    .map((c) => ({ chrom: c, count: chrCounts[c] }));

  const geneNames = Array.from(new Set(rows.map((r) => r.hugo || r.gene || "").filter(Boolean))).slice(0, 30);
  const BIOTYPE_COLORS = ["#4f46e5","#06b6d4","#10b981","#f59e0b","#ef4444","#a855f7","#22c55e","#eab308","#3b82f6","#f97316"];

  return (
    <main className="min-h-screen px-6 py-10 max-w-7xl mx-auto">
      <h1 className="text-3xl font-semibold mb-2">annota-bed</h1>
      <p className="italic text-gray-500 mb-6">
        Fast and easy annotation of BED files for contextual genomic analysis
      </p>

      {/* Uploads + options */}
      <form onSubmit={handleSubmit} className="rounded-2xl border p-5 space-y-5 bg-white/50">
        <div className="grid gap-5 md:grid-cols-3">
          <div className="md:col-span-3">
            <label className="block text-sm font-medium mb-1">Upload BED (chr, start, end)</label>
            <input
              type="file" accept=".bed,.tsv,.csv,.txt"
              onChange={(e) => setBedFile(e.target.files?.[0] || null)}
              className="w-full cursor-pointer file:mr-4 file:py-2 file:px-4 file:rounded-lg file:border-0 file:bg-gray-100 file:text-sm"
              required
            />
          </div>

          <div>
            <label className="block text-sm font-medium mb-1">GTF (bgzip, e.g. *.gtf.bgz)</label>
            <input
              type="file" accept=".bgz,.gz"
              onChange={(e) => setGtfBgz(e.target.files?.[0] || null)}
              className="w-full cursor-pointer file:mr-4 file:py-2 file:px-4 file:rounded-lg file:border-0 file:bg-gray-100 file:text-sm"
              required
            />
            <p className="text-xs text-gray-500 mt-1">
              Must be bgzip-compressed and coordinate-sorted.
            </p>
          </div>

          <div>
            <label className="block text-sm font-medium mb-1">Tabix index (.tbi)</label>
            <input
              type="file" accept=".tbi"
              onChange={(e) => setGtfTbi(e.target.files?.[0] || null)}
              className="w-full cursor-pointer file:mr-4 file:py-2 file:px-4 file:rounded-lg file:border-0 file:bg-gray-100 file:text-sm"
              required
            />
            <p className="text-xs text-gray-500 mt-1">
              Must match the uploaded GTF (same basename, created by <code>tabix -p gff</code>).
            </p>
          </div>

          <div>
            <label className="block text-sm font-medium mb-1">Ambiguities</label>
            <select
              value={ambiguities}
              onChange={(e) => setAmbiguities(e.target.value)}
              className="w-full border rounded-md p-2 text-black bg-white"
            >
              <option value="best_all">best_all (default)</option>
              <option value="best_one">best_one</option>
              <option value="all">all</option>
            </select>
          </div>
        </div>

        <button
          type="submit"
          disabled={loading || !bedFile || !gtfBgz || !gtfTbi}
          className="px-4 py-2 rounded-lg bg-black text-white disabled:opacity-50"
        >
          {loading ? "Annotating…" : "Annotate"}
        </button>

        {error && <p className="text-red-600 text-sm">{error}</p>}
      </form>

      {/* Summary + plots */}
      {summary && (
        <section className="mt-8 space-y-6">
          <div className="flex items-center justify-between">
            <h2 className="text-xl font-semibold">Summary</h2>
            <div className="text-sm">
              <span className="font-medium">Unique genes: </span>
              {summary.unique_genes}
            </div>
          </div>

          {/* Two plots side-by-side (bigger boxes + unclipped labels) */}
          <div className="grid gap-6 md:grid-cols-1 lg:grid-cols-2">
            {/* Biotype */}
            <div className="border rounded-2xl p-4 bg-white h-[32rem] min-h-[28rem]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={biotypeChartData} margin={{ top: 20, right: 30, left: 30, bottom: 100 }}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="biotype" interval={0} angle={-40} textAnchor="end" height={100} tick={{ fontSize: 14 }}>
                    <Label value="Biotype" position="bottom" offset={30} style={{ textAnchor: "middle" }} />
                  </XAxis>
                  <YAxis tick={{ fontSize: 14 }}>
                    <Label value="Count" angle={-90} position="insideLeft" style={{ textAnchor: "middle" }} />
                  </YAxis>
                  <Tooltip />
                  <Bar dataKey="count">
                    {biotypeChartData.map((_, idx) => (
                      <Cell key={idx} fill={BIOTYPE_COLORS[idx % BIOTYPE_COLORS.length]} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>

            {/* Chromosome distribution */}
            <div className="border rounded-2xl p-4 bg-white h-[32rem] min-h-[28rem]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={chromChartData} margin={{ top: 20, right: 30, left: 30, bottom: 100 }}>
                  <CartesianGrid strokeDasharray="3 3" />
                  <XAxis dataKey="chrom" interval={0} angle={-40} textAnchor="end" height={100} tick={{ fontSize: 14 }}>
                    <Label value="Chromosome" position="bottom" offset={30} style={{ textAnchor: "middle" }} />
                  </XAxis>
                  <YAxis tick={{ fontSize: 14 }}>
                    <Label value="Row count" angle={-90} position="insideLeft" style={{ textAnchor: "middle" }} />
                  </YAxis>
                  <Tooltip />
                  <Bar dataKey="count" fill="#475569" />
                </BarChart>
              </ResponsiveContainer>
            </div>
          </div>

          {/* Gene list */}
          {geneNames.length > 0 && (
            <div className="border rounded-2xl p-4 bg-white">
              <div className="text-sm text-gray-600 font-medium mb-1">Example genes:</div>
              <div className="text-sm text-gray-700 whitespace-nowrap overflow-x-auto">
                {geneNames.join(", ")}
              </div>
            </div>
          )}

          {/* Download */}
          {csvPath && (
            <a
              href={`/api/download?path=${encodeURIComponent(csvPath)}`}
              className="inline-block px-4 py-2 rounded-lg bg-gray-800 text-white"
            >
              Download CSV
            </a>
          )}
        </section>
      )}

      {/* Results table */}
      {rows.length > 0 && (
        <section className="mt-8">
          <h2 className="text-xl font-semibold mb-3">Results</h2>
          <div className="border rounded-2xl bg-white max-h=[520px] overflow-y-auto">
            <table className="min-w-full text-sm text-gray-800">
              <thead className="bg-gray-100 sticky top-0 z-10">
                <tr>
                  {[
                    "chr","start","end","hugo/gene","strand",
                    "biotype","tx_overlap%","exon_overlap%","cds_overlap%","transcript_id",
                  ].map((h) => (
                    <th key={h} className="text-left font-medium p-2">{h}</th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {rows.map((r, i) => (
                  <tr key={i} className="odd:bg-white even:bg-gray-50">
                    <td className="p-2">{r.input_chr}</td>
                    <td className="p-2">{r.input_start}</td>
                    <td className="p-2">{r.input_end}</td>
                    <td className="p-2">{r.hugo ?? r.gene ?? "-"}</td>
                    <td className="p-2">{r.strand ?? "-"}</td>
                    <td className="p-2">{r.feature_biotype ?? "-"}</td>
                    <td className="p-2">{r.tx_overlap_pct}</td>
                    <td className="p-2">{r.exon_overlap_pct}</td>
                    <td className="p-2">{r.cds_overlap_pct}</td>
                    <td className="p-2">{r.ensembl_id ?? "-"}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </section>
      )}
    </main>
  );
}