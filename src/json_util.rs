// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Zero-dependency reproducible structured JSON summary export helpers for batch
//! statistical, biodiversity (Darwin Core / eDNA), population-genetics, and CRISPR runs.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::crispr_core::{CrisprGuideCandidate, CrisprNuclease};
use crate::dwc_core::DwcDataset;
use crate::edna_diversity::CommunityDiversityReport;
use crate::popgen_core::PopGenFastaReport;

/// Escapes reserved JSON string characters (`"`, `\`, control bytes) per RFC 8259.
pub fn json_escape(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len() + 8);
    for ch in raw.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

/// Formats a floating-point value for valid RFC 8259 JSON (`null` for `NaN` or `±Inf`).
pub fn format_json_f64(value: f64) -> String {
    if value.is_finite() {
        let normalized = if value == 0.0 { 0.0 } else { value };
        format!("{normalized:.6}")
    } else {
        String::from("null")
    }
}

/// Typed scalar value inside a reproducible [`BatchRunSummary`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BatchMetricValue {
    Int(i64),
    Usize(usize),
    Float(f64),
    Bool(bool),
    Text(String),
}

impl BatchMetricValue {
    fn to_json_literal(&self) -> String {
        match self {
            Self::Int(v) => v.to_string(),
            Self::Usize(v) => v.to_string(),
            Self::Float(v) => format_json_f64(*v),
            Self::Bool(v) => v.to_string(),
            Self::Text(s) => format!("\"{}\"", json_escape(s)),
        }
    }
}

/// Deterministic structured summary record for batch OSABT runs (`BTreeMap`-ordered keys).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BatchRunSummary {
    pub schema_version: String,
    pub tool: String,
    pub dataset_label: String,
    pub record_count: usize,
    pub metrics: BTreeMap<String, BatchMetricValue>,
    pub notes: Vec<String>,
}

impl BatchRunSummary {
    /// Creates a new reproducible batch summary envelope.
    pub fn new(
        tool: impl Into<String>,
        dataset_label: impl Into<String>,
        record_count: usize,
    ) -> Self {
        Self {
            schema_version: String::from("osabt-batch-summary/v1"),
            tool: tool.into(),
            dataset_label: dataset_label.into(),
            record_count,
            metrics: BTreeMap::new(),
            notes: Vec::new(),
        }
    }

    /// Inserts or replaces a named scalar metric (sorted deterministically by key).
    pub fn with_metric(mut self, key: impl Into<String>, value: BatchMetricValue) -> Self {
        self.metrics.insert(key.into(), value);
        self
    }

    /// Appends a diagnostic or provenance note.
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }

    /// Serializes the summary to deterministic pretty-printed JSON.
    pub fn to_json(&self) -> String {
        let mut out = String::with_capacity(512);
        let _ = writeln!(out, "{{");
        let _ = writeln!(
            out,
            "  \"schema_version\": \"{}\",",
            json_escape(&self.schema_version)
        );
        let _ = writeln!(out, "  \"tool\": \"{}\",", json_escape(&self.tool));
        let _ = writeln!(
            out,
            "  \"dataset_label\": \"{}\",",
            json_escape(&self.dataset_label)
        );
        let _ = writeln!(out, "  \"record_count\": {},", self.record_count);
        let _ = writeln!(out, "  \"metrics\": {{");
        let metric_len = self.metrics.len();
        for (idx, (k, v)) in self.metrics.iter().enumerate() {
            let comma = if idx + 1 < metric_len { "," } else { "" };
            let _ = writeln!(
                out,
                "    \"{}\": {}{}",
                json_escape(k),
                v.to_json_literal(),
                comma
            );
        }
        let _ = writeln!(out, "  }},");
        let _ = writeln!(out, "  \"notes\": [");
        let note_len = self.notes.len();
        for (idx, note) in self.notes.iter().enumerate() {
            let comma = if idx + 1 < note_len { "," } else { "" };
            let _ = writeln!(out, "    \"{}\"{}", json_escape(note), comma);
        }
        let _ = writeln!(out, "  ]");
        let _ = writeln!(out, "}}");
        out
    }
}

/// Exports a reproducible structured JSON summary for a batch statistical run.
pub fn export_statistical_summary_json(
    tool: &str,
    dataset_label: &str,
    record_count: usize,
    metrics: &[(&str, f64)],
    notes: &[String],
) -> String {
    let mut summary = BatchRunSummary::new(tool, dataset_label, record_count);
    for &(k, v) in metrics {
        summary = summary.with_metric(k, BatchMetricValue::Float(v));
    }
    for n in notes {
        summary = summary.with_note(n.clone());
    }
    summary.to_json()
}

/// Exports a reproducible structured JSON audit summary for a Darwin Core biodiversity dataset.
pub fn export_biodiversity_summary_json(dataset: &DwcDataset) -> String {
    let s = dataset.summary();
    let mut issue_counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for iss in &dataset.issues {
        *issue_counts.entry(iss.code.as_str()).or_insert(0) += 1;
    }

    let mut out = String::with_capacity(1024);
    let _ = writeln!(out, "{{");
    let _ = writeln!(out, "  \"schema_version\": \"osabt-dwc-summary/v1\",");
    let _ = writeln!(out, "  \"tool\": \"dwc\",");
    let _ = writeln!(
        out,
        "  \"dataset_name\": \"{}\",",
        json_escape(&dataset.dataset_name)
    );
    let _ = writeln!(out, "  \"total_records\": {},", s.total_records);
    let _ = writeln!(out, "  \"clean_records\": {},", s.clean_records);
    let _ = writeln!(
        out,
        "  \"records_with_warnings\": {},",
        s.records_with_warnings
    );
    let _ = writeln!(out, "  \"records_with_errors\": {},", s.records_with_errors);
    let _ = writeln!(
        out,
        "  \"georeferenced_records\": {},",
        s.georeferenced_records
    );
    let _ = writeln!(out, "  \"auto_fixed_count\": {},", s.auto_fixed_count);
    let _ = writeln!(out, "  \"unique_taxa\": {},", s.unique_taxa);
    let _ = writeln!(out, "  \"synonyms_resolved\": {},", s.synonyms_resolved);
    let _ = writeln!(out, "  \"issue_breakdown\": {{");
    let n_codes = issue_counts.len();
    for (idx, (code, count)) in issue_counts.iter().enumerate() {
        let comma = if idx + 1 < n_codes { "," } else { "" };
        let _ = writeln!(out, "    \"{}\": {}{}", json_escape(code), count, comma);
    }
    let _ = writeln!(out, "  }}");
    let _ = writeln!(out, "}}");
    out
}

/// Exports a reproducible structured JSON summary for a batch CRISPR guide scan & `CFD` run.
pub fn export_crispr_summary_json(
    dataset_label: &str,
    nuclease: CrisprNuclease,
    target_len_bp: usize,
    guides: &[CrisprGuideCandidate],
    max_guides: usize,
) -> String {
    let mut ranked: Vec<&CrisprGuideCandidate> = guides.iter().collect();
    ranked.sort_by(|a, b| {
        b.composite_score
            .total_cmp(&a.composite_score)
            .then_with(|| a.guide_id.cmp(&b.guide_id))
    });

    let zero_offtarget = guides.iter().filter(|g| g.off_target_count == 0).count();
    let cbe_clean = guides
        .iter()
        .filter(|g| g.base_edit.cbe_clean_single_edit)
        .count();
    let abe_clean = guides
        .iter()
        .filter(|g| g.base_edit.abe_clean_single_edit)
        .count();
    let top_n: Vec<&CrisprGuideCandidate> = ranked.into_iter().take(max_guides.max(1)).collect();

    let mut out = String::with_capacity(1024);
    let _ = writeln!(out, "{{");
    let _ = writeln!(out, "  \"schema_version\": \"osabt-crispr-summary/v1\",");
    let _ = writeln!(out, "  \"tool\": \"crispr\",");
    let _ = writeln!(
        out,
        "  \"dataset_label\": \"{}\",",
        json_escape(dataset_label)
    );
    let _ = writeln!(
        out,
        "  \"nuclease\": \"{}\",",
        json_escape(nuclease.label())
    );
    let _ = writeln!(out, "  \"target_len_bp\": {},", target_len_bp);
    let _ = writeln!(out, "  \"candidate_count\": {},", guides.len());
    let _ = writeln!(out, "  \"zero_offtarget_count\": {},", zero_offtarget);
    let _ = writeln!(out, "  \"cbe_clean_count\": {},", cbe_clean);
    let _ = writeln!(out, "  \"abe_clean_count\": {},", abe_clean);
    let _ = writeln!(out, "  \"top_guides\": [");
    let count = top_n.len();
    for (idx, g) in top_n.iter().enumerate() {
        let comma = if idx + 1 < count { "," } else { "" };
        let _ = writeln!(
            out,
            "    {{\"guide_id\": \"{}\", \"strand\": \"{}\", \"cut_site_1idx\": {}, \"spacer_seq\": \"{}\", \"pam_seq\": \"{}\", \"gc_percent\": {}, \"on_target_score\": {}, \"cfd_specificity_score\": {}, \"off_target_count\": {}, \"composite_score\": {}}}{}",
            json_escape(&g.guide_id),
            g.strand,
            g.cut_site_1idx,
            json_escape(&g.spacer_seq),
            json_escape(&g.pam_seq),
            format_json_f64(g.gc_percent),
            format_json_f64(g.on_target_score),
            format_json_f64(g.cfd_specificity_score),
            g.off_target_count,
            format_json_f64(g.composite_score),
            comma
        );
    }
    let _ = writeln!(out, "  ]");
    let _ = writeln!(out, "}}");
    out
}

/// Exports a reproducible structured JSON summary for an eDNA community diversity report.
pub fn export_edna_diversity_summary_json(
    dataset_label: &str,
    report: &CommunityDiversityReport,
) -> String {
    let mut out = String::with_capacity(1024);
    let _ = writeln!(out, "{{");
    let _ = writeln!(out, "  \"schema_version\": \"osabt-edna-diversity/v1\",");
    let _ = writeln!(out, "  \"tool\": \"edna\",");
    let _ = writeln!(
        out,
        "  \"dataset_label\": \"{}\",",
        json_escape(dataset_label)
    );
    let _ = writeln!(out, "  \"clean_asv_count\": {},", report.clean_asv_count);
    let _ = writeln!(out, "  \"samples\": [");
    let n = report.alpha.len();
    for (idx, s) in report.alpha.iter().enumerate() {
        let comma = if idx + 1 < n { "," } else { "" };
        let _ = writeln!(
            out,
            "    {{\"sample_name\": \"{}\", \"total_reads_n\": {}, \"observed_richness_q0\": {}, \"chao1_asymptotic\": {}, \"hill_shannon_q1\": {}, \"hill_simpson_q2\": {}, \"hurlbert_pie\": {}, \"sample_coverage\": {}}}{}",
            json_escape(&s.sample_name),
            s.total_reads_n,
            s.observed_richness_q0,
            format_json_f64(s.chao1_asymptotic),
            format_json_f64(s.hill_shannon_q1),
            format_json_f64(s.hill_simpson_q2),
            format_json_f64(s.hurlbert_pie),
            format_json_f64(s.sample_coverage),
            comma
        );
    }
    let _ = writeln!(out, "  ]");
    let _ = writeln!(out, "}}");
    out
}

/// Exports a reproducible structured JSON summary for a PopGen FASTA diversity report.
pub fn export_popgen_summary_json(dataset_label: &str, report: &PopGenFastaReport) -> String {
    let g = &report.global_metrics;
    BatchRunSummary::new("popgen", dataset_label, report.total_sequences)
        .with_metric(
            "alignment_len_bp",
            BatchMetricValue::Usize(report.alignment_len_bp),
        )
        .with_metric(
            "population_count",
            BatchMetricValue::Usize(report.populations.len()),
        )
        .with_metric(
            "num_haplotypes",
            BatchMetricValue::Usize(g.num_haplotypes),
        )
        .with_metric(
            "haplotype_diversity_hd",
            BatchMetricValue::Float(g.haplotype_diversity_hd),
        )
        .with_metric(
            "segregating_sites_s",
            BatchMetricValue::Usize(g.segregating_sites_s),
        )
        .with_metric(
            "nucleotide_diversity_pi",
            BatchMetricValue::Float(g.nucleotide_diversity_pi),
        )
        .with_metric(
            "watterson_theta_w",
            BatchMetricValue::Float(g.watterson_theta_w),
        )
        .with_metric("tajima_d", BatchMetricValue::Float(g.tajima_d))
        .with_metric(
            "tajima_d_pvalue",
            BatchMetricValue::Float(g.tajima_d_pvalue),
        )
        .to_json()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crispr_core::{demo_crispr_target_and_offtargets, scan_crispr_candidates};
    use crate::dwc_core::{lint_and_clean_dwc, DwcCleanConfig, DwcDataset};
    use crate::dwc_taxonomy::TaxonBackbone;

    #[test]
    fn batch_summary_json_is_deterministic_and_nan_safe() {
        let a = export_statistical_summary_json(
            "flux",
            "leaf_aci_batch",
            12,
            &[
                ("vcmax", 85.73),
                ("nan_metric", f64::NAN),
                ("inf_metric", f64::INFINITY),
            ],
            &[String::from("line\"1\nline2")],
        );
        let b = export_statistical_summary_json(
            "flux",
            "leaf_aci_batch",
            12,
            &[
                ("vcmax", 85.73),
                ("nan_metric", f64::NAN),
                ("inf_metric", f64::INFINITY),
            ],
            &[String::from("line\"1\nline2")],
        );
        assert_eq!(a, b);
        assert!(a.contains("\"nan_metric\": null"));
        assert!(a.contains("\"inf_metric\": null"));
        assert!(a.contains("\"vcmax\": 85.730000"));
        assert!(a.contains("line\\\"1\\nline2"));
    }

    #[test]
    fn exports_biodiversity_and_crispr_summaries() {
        let bb = TaxonBackbone::default_curated();
        let ds = lint_and_clean_dwc(
            "PNW_Survey",
            DwcDataset::demo_pnw_field_survey_csv(),
            &bb,
            &DwcCleanConfig::default(),
        )
        .expect("demo DwC dataset should lint cleanly");
        let dwc_json = export_biodiversity_summary_json(&ds);
        assert!(dwc_json.contains("\"schema_version\": \"osabt-dwc-summary/v1\""));
        assert!(dwc_json.contains("\"dataset_name\": \"PNW_Survey\""));
        assert!(dwc_json.contains("\"issue_breakdown\":"));

        let (tgt, _) = demo_crispr_target_and_offtargets();
        let (seq, guides) =
            scan_crispr_candidates(&tgt, CrisprNuclease::SpCas9).expect("demo CRISPR scan");
        let crispr_json =
            export_crispr_summary_json("DemoLocus", CrisprNuclease::SpCas9, seq.len(), &guides, 3);
        assert!(crispr_json.contains("\"schema_version\": \"osabt-crispr-summary/v1\""));
        assert!(crispr_json.contains("\"top_guides\": ["));
        assert!(crispr_json.contains("\"candidate_count\":"));
    }
}
