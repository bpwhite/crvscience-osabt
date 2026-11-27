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

//! Command-line subcommands (`clap` v4) for Tool #3: `osabt dwc`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::dwc_core::{lint_and_clean_dwc, DwcCleanConfig, DwcDataset};
use crate::dwc_export::{
    export_dwc_geojson, export_issues_tsv, generate_dwca_eml_xml, generate_dwca_meta_xml,
    render_dwc_occurrence_map_svg, write_dwc_delimited,
};
use crate::dwc_taxonomy::{resolve_scientific_name, TaxonBackbone};

#[derive(Debug, Args)]
pub struct DwcArgs {
    #[command(subcommand)]
    pub command: Option<DwcSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum DwcSubcommand {
    /// Launch the interactive Ratatui Darwin Core & field-survey workbench.
    Tui {
        /// Optional input CSV/TSV file (launches built-in Oregon PNW demo survey if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Optional local taxonomic checklist CSV/TSV to merge into the offline backbone.
        #[arg(short, long)]
        checklist: Option<PathBuf>,
    },

    /// Audit a field CSV/TSV against TDWG Darwin Core terms, CoordinateCleaner rules, and taxonomy.
    Lint {
        /// Input CSV/TSV file (uses built-in Oregon PNW demo survey if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Optional local taxonomic checklist CSV/TSV.
        #[arg(short, long)]
        checklist: Option<PathBuf>,
        /// Optional path to save the tab-separated QA diagnostic ledger (.tsv).
        #[arg(short, long)]
        report_out: Option<PathBuf>,
        /// Optional path to save a reproducible structured JSON audit summary (.json).
        #[arg(long)]
        json_out: Option<PathBuf>,
    },

    /// Scrub headers, DMS coordinates, transposed lat/lon, ISO dates, and synonyms into clean CSV.
    Clean {
        /// Input CSV/TSV file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output scrubbed Simple Darwin Core CSV path.
        #[arg(short, long, default_value = "cleaned_occurrences.csv")]
        output: PathBuf,
        /// Replace synonymous `scientificName` values in place rather than keeping original + `acceptedNameUsage`.
        #[arg(long)]
        replace_synonyms: bool,
        /// Retain rows with unresolved spatial/schema errors instead of filtering them out.
        #[arg(long)]
        keep_errors: bool,
        /// Optional local taxonomic checklist CSV/TSV.
        #[arg(short, long)]
        checklist: Option<PathBuf>,
    },

    /// Package a complete GBIF IPT-ready Darwin Core Archive directory (`occurrence.csv`, `meta.xml`, `eml.xml`, `.geojson`, `.svg`).
    Archive {
        /// Input CSV/TSV file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output directory for the Darwin Core Archive files.
        #[arg(short, long, default_value = "osabt_dwca_bundle")]
        out_dir: PathBuf,
        /// Dataset title recorded in `eml.xml`.
        #[arg(long, default_value = "CRV Science Biodiversity Occurrence Survey")]
        title: String,
        /// Creator surname/team recorded in `eml.xml`.
        #[arg(long, default_value = "White")]
        creator: String,
        /// Optional local taxonomic checklist CSV/TSV.
        #[arg(short, long)]
        checklist: Option<PathBuf>,
    },

    /// Parse and resolve one or more scientific names (binomials, authorships, synonyms, typos, GBIF match URLs).
    Resolve {
        /// Scientific name strings to resolve (e.g. `"Pseudotsuga taxifolia (Poir.) Britton"`).
        #[arg(required = true)]
        names: Vec<String>,
        /// Optional local taxonomic checklist CSV/TSV.
        #[arg(short, long)]
        checklist: Option<PathBuf>,
    },

    /// Render a standalone vector `.svg` spatial QA distribution map (and optional `.geojson`).
    Map {
        /// Input CSV/TSV file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output `.svg` map path.
        #[arg(short, long, default_value = "dwc_spatial_map.svg")]
        output: PathBuf,
        /// Optional output `.geojson` path for QGIS / web mapping.
        #[arg(long)]
        geojson_out: Option<PathBuf>,
    },
}

pub fn run_dwc_cli(args: DwcArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "tui")]
            {
                return crate::dwc_tui::run_dwc_tui(None, None);
            }
            #[cfg(not(feature = "tui"))]
            {
                return Err(String::from(
                    "Build with `--features tui` or run `osabt dwc lint`.",
                ));
            }
        }
        Some(DwcSubcommand::Tui { input, checklist }) => {
            #[cfg(feature = "tui")]
            {
                crate::dwc_tui::run_dwc_tui(input, checklist)
            }
            #[cfg(not(feature = "tui"))]
            {
                let _ = (input, checklist);
                Err(String::from("TUI feature is disabled in this build."))
            }
        }
        Some(DwcSubcommand::Lint {
            input,
            checklist,
            report_out,
            json_out,
        }) => {
            let backbone = load_backbone(checklist.as_ref())?;
            let (name, raw) = load_input_or_demo(input.as_ref())?;
            let ds = lint_and_clean_dwc(&name, &raw, &backbone, &DwcCleanConfig::default())?;
            let s = ds.summary();

            println!("=== OSABT Darwin Core QA Audit: {} ===", ds.dataset_name);
            println!(
                "Rows: {} total | {} clean | {} warn | {} error | {} georeferenced | {} auto-fixed",
                s.total_records,
                s.clean_records,
                s.records_with_warnings,
                s.records_with_errors,
                s.georeferenced_records,
                s.auto_fixed_count
            );
            println!("---------------------------------------------------------------------");
            for iss in &ds.issues {
                println!(
                    "[Row {:>3}] {:<5} {:<22} {:<22} {}",
                    iss.row_number_1,
                    iss.severity.label(),
                    iss.code.as_str(),
                    iss.column_term,
                    iss.message
                );
            }
            if let Some(rp) = report_out {
                fs::write(&rp, export_issues_tsv(&ds))
                    .map_err(|e| format!("Failed writing report '{}': {e}", rp.display()))?;
                println!("Saved TSV audit ledger -> {}", rp.display());
            }
            if let Some(jp) = json_out {
                fs::write(&jp, crate::json_util::export_biodiversity_summary_json(&ds))
                    .map_err(|e| format!("Failed writing JSON summary '{}': {e}", jp.display()))?;
                println!("Saved JSON audit summary -> {}", jp.display());
            }
            Ok(())
        }
        Some(DwcSubcommand::Clean {
            input,
            output,
            replace_synonyms,
            keep_errors,
            checklist,
        }) => {
            let backbone = load_backbone(checklist.as_ref())?;
            let (name, raw) = load_input_or_demo(input.as_ref())?;
            let cfg = DwcCleanConfig {
                replace_synonyms_in_place: replace_synonyms,
                ..DwcCleanConfig::default()
            };
            let ds = lint_and_clean_dwc(&name, &raw, &backbone, &cfg)?;
            let csv = write_dwc_delimited(&ds, ',', !keep_errors);
            fs::write(&output, csv)
                .map_err(|e| format!("Failed writing '{}': {e}", output.display()))?;
            let s = ds.summary();
            println!(
                "Wrote scrubbed Simple Darwin Core CSV -> {} ({} auto-fixes applied, {} error rows {}).",
                output.display(),
                s.auto_fixed_count,
                s.records_with_errors,
                if keep_errors { "retained" } else { "filtered" }
            );
            Ok(())
        }
        Some(DwcSubcommand::Archive {
            input,
            out_dir,
            title,
            creator,
            checklist,
        }) => {
            let backbone = load_backbone(checklist.as_ref())?;
            let (name, raw) = load_input_or_demo(input.as_ref())?;
            let ds = lint_and_clean_dwc(&name, &raw, &backbone, &DwcCleanConfig::default())?;
            fs::create_dir_all(&out_dir)
                .map_err(|e| format!("Failed creating '{}': {e}", out_dir.display()))?;

            fs::write(
                out_dir.join("occurrence.csv"),
                write_dwc_delimited(&ds, ',', true),
            )
            .map_err(|e| format!("Failed writing occurrence.csv: {e}"))?;
            fs::write(
                out_dir.join("meta.xml"),
                generate_dwca_meta_xml(&ds, "occurrence.csv", Some("eml.xml")),
            )
            .map_err(|e| format!("Failed writing meta.xml: {e}"))?;
            fs::write(
                out_dir.join("eml.xml"),
                generate_dwca_eml_xml(&title, &creator, "CC-BY-4.0"),
            )
            .map_err(|e| format!("Failed writing eml.xml: {e}"))?;
            fs::write(out_dir.join("qa_audit.tsv"), export_issues_tsv(&ds))
                .map_err(|e| format!("Failed writing qa_audit.tsv: {e}"))?;
            fs::write(out_dir.join("occurrences.geojson"), export_dwc_geojson(&ds))
                .map_err(|e| format!("Failed writing occurrences.geojson: {e}"))?;
            fs::write(
                out_dir.join("spatial_qa_map.svg"),
                render_dwc_occurrence_map_svg(&ds),
            )
            .map_err(|e| format!("Failed writing spatial_qa_map.svg: {e}"))?;

            println!(
                "Packaged GBIF IPT-ready Darwin Core Archive bundle -> {}/ (occurrence.csv, meta.xml, eml.xml, qa_audit.tsv, occurrences.geojson, spatial_qa_map.svg)",
                out_dir.display()
            );
            Ok(())
        }
        Some(DwcSubcommand::Resolve { names, checklist }) => {
            let backbone = load_backbone(checklist.as_ref())?;
            for raw_name in &names {
                let res = resolve_scientific_name(raw_name, &backbone);
                println!("Input          : {}", res.parsed.raw_input);
                println!("Canonical Name : {}", res.parsed.canonical_name);
                println!("Authorship     : {}", res.parsed.authorship);
                println!("Rank           : {}", res.parsed.inferred_rank);
                println!("Match Kind     : {:?}", res.match_kind);
                if let Some(m) = res.matched_entry {
                    println!(
                        "Accepted Usage : {} {} [GBIF:{} | {} > {} > {}]",
                        m.accepted_name, m.authorship, m.gbif_usage_key, m.kingdom, m.family, m.genus
                    );
                }
                println!("GBIF Match URL : {}", res.gbif_lookup_url);
                println!("---");
            }
            Ok(())
        }
        Some(DwcSubcommand::Map {
            input,
            output,
            geojson_out,
        }) => {
            let backbone = TaxonBackbone::default_curated();
            let (name, raw) = load_input_or_demo(input.as_ref())?;
            let ds = lint_and_clean_dwc(&name, &raw, &backbone, &DwcCleanConfig::default())?;
            fs::write(&output, render_dwc_occurrence_map_svg(&ds))
                .map_err(|e| format!("Failed writing SVG map '{}': {e}", output.display()))?;
            println!("Saved Darwin Core spatial QA map -> {}", output.display());
            if let Some(gp) = geojson_out {
                fs::write(&gp, export_dwc_geojson(&ds))
                    .map_err(|e| format!("Failed writing GeoJSON '{}': {e}", gp.display()))?;
                println!("Saved RFC-7946 GeoJSON -> {}", gp.display());
            }
            Ok(())
        }
    }
}

fn load_backbone(checklist: Option<&PathBuf>) -> Result<TaxonBackbone, String> {
    let mut bb = TaxonBackbone::default_curated();
    if let Some(p) = checklist {
        let txt = fs::read_to_string(p)
            .map_err(|e| format!("Failed reading checklist '{}': {e}", p.display()))?;
        bb.load_checklist_delimited(&txt)?;
    }
    Ok(bb)
}

fn load_input_or_demo(input: Option<&PathBuf>) -> Result<(String, String), String> {
    match input {
        Some(p) => {
            let txt = fs::read_to_string(p)
                .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
            let stem = p
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| String::from("Survey"));
            Ok((stem, txt))
        }
        None => Ok((
            String::from("PNW_Oregon_Demo_Survey"),
            DwcDataset::demo_pnw_field_survey_csv().to_string(),
        )),
    }
}
