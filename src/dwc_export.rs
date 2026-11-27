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

//! Tool #3 (`osabt dwc`) — Exporters for scrubbed RFC-4180 Simple Darwin Core CSV/TSV,
//! TDWG Darwin Core Archive (`meta.xml` + `eml.xml`), RFC-7946 `.geojson` GIS layers,
//! audit report TSV, and standalone vector `.svg` occurrence distribution maps.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::dwc_core::{dwc_term_iri, DwcDataset};
use crate::json_util::json_escape;
use crate::svg_util::{begin_svg_canvas, xml_escape};

/// Serializes a [`DwcDataset`] to RFC-4180 CSV (or TSV when `delimiter == '\t'`).
/// If `exclude_error_rows` is true, rows with unresolved `ERROR` diagnostics (such as Null Island)
/// are omitted from the clean publication table.
pub fn write_dwc_delimited(
    dataset: &DwcDataset,
    delimiter: char,
    exclude_error_rows: bool,
) -> String {
    let mut out = String::with_capacity(dataset.records.len() * 128 + 256);

    for (idx, col) in dataset.columns.iter().enumerate() {
        if idx > 0 {
            out.push(delimiter);
        }
        push_escaped_cell(&mut out, col, delimiter);
    }
    out.push('\n');

    for rec in &dataset.records {
        if exclude_error_rows && rec.has_errors {
            continue;
        }
        for (idx, col) in dataset.columns.iter().enumerate() {
            if idx > 0 {
                out.push(delimiter);
            }
            push_escaped_cell(&mut out, rec.get(col), delimiter);
        }
        out.push('\n');
    }

    out
}

fn push_escaped_cell(buf: &mut String, val: &str, delimiter: char) {
    let needs_quotes =
        val.contains(delimiter) || val.contains('"') || val.contains('\n') || val.contains('\r');
    if needs_quotes {
        buf.push('"');
        for ch in val.chars() {
            if ch == '"' {
                buf.push_str("\"\"");
            } else {
                buf.push(ch);
            }
        }
        buf.push('"');
    } else {
        buf.push_str(val);
    }
}

/// Generates a normative TDWG Darwin Core Text Guide `meta.xml` descriptor
/// (`http://rs.tdwg.org/dwc/text/`) mapping every recognized Darwin Core column in `dataset`.
pub fn generate_dwca_meta_xml(
    dataset: &DwcDataset,
    core_filename: &str,
    eml_filename: Option<&str>,
) -> String {
    let mut xml = String::with_capacity(2048);
    xml.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    xml.push_str("<archive xmlns=\"http://rs.tdwg.org/dwc/text/\"\n");
    xml.push_str("         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
    xml.push_str("         xmlns:xs=\"http://www.w3.org/2001/XMLSchema\"\n");
    xml.push_str("         xsi:schemaLocation=\"http://rs.tdwg.org/dwc/text/ http://rs.tdwg.org/dwc/text/tdwg_dwc_text.xsd\"");
    if let Some(eml) = eml_filename {
        let _ = write!(xml, "\n         metadata=\"{}\"", xml_escape(eml));
    }
    xml.push_str(">\n");

    let id_index = dataset
        .columns
        .iter()
        .position(|c| c == "occurrenceID")
        .unwrap_or(0);

    let _ = writeln!(
        xml,
        "  <core encoding=\"UTF-8\" fieldsTerminatedBy=\",\" linesTerminatedBy=\"\\n\" fieldsEnclosedBy=\"&quot;\" ignoreHeaderLines=\"1\" rowType=\"http://rs.tdwg.org/dwc/terms/Occurrence\">"
    );
    xml.push_str("    <files>\n");
    let _ = writeln!(xml, "      <location>{}</location>", xml_escape(core_filename));
    xml.push_str("    </files>\n");
    let _ = writeln!(xml, "    <id index=\"{id_index}\" />");

    for (idx, col) in dataset.columns.iter().enumerate() {
        let iri = dwc_term_iri(col)
            .map(String::from)
            .unwrap_or_else(|| format!("https://crvscience.com/dwc/custom/{}", xml_escape(col)));
        let _ = writeln!(xml, "    <field index=\"{idx}\" term=\"{iri}\" />");
    }

    xml.push_str("  </core>\n");
    xml.push_str("</archive>\n");
    xml
}

/// Generates a minimal GBIF-compatible Ecological Metadata Language (`eml.xml`) stub
/// for packaging a Darwin Core Archive directory.
pub fn generate_dwca_eml_xml(title: &str, creator_name: &str, license_spdx: &str) -> String {
    format!(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<eml:eml packageId=\"osabt-dwca-1.0\" system=\"https://crvscience.com\"\n\
         xmlns:eml=\"eml://ecoinformatics.org/eml-2.1.1\"\n\
         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n\
  <dataset>\n\
    <title>{}</title>\n\
    <creator>\n\
      <individualName>\n\
        <surName>{}</surName>\n\
      </individualName>\n\
      <organizationName>CRV Science Network</organizationName>\n\
    </creator>\n\
    <intellectualRights>\n\
      <para>Published via OSABT (GPL-3.0). Dataset license: {}.</para>\n\
    </intellectualRights>\n\
  </dataset>\n\
</eml:eml>\n",
        xml_escape(title),
        xml_escape(creator_name),
        xml_escape(license_spdx)
    )
}

/// Exports all valid georeferenced records in `dataset` as an RFC-7946 `FeatureCollection` `.geojson`.
/// Note RFC-7946 coordinate order: `[longitude, latitude]`.
pub fn export_dwc_geojson(dataset: &DwcDataset) -> String {
    let mut out = String::from("{\n  \"type\": \"FeatureCollection\",\n  \"features\": [\n");
    let mut first = true;

    for rec in &dataset.records {
        let Some((lat, lon)) = rec.valid_lat_lon else {
            continue;
        };
        if !first {
            out.push_str(",\n");
        }
        first = false;

        let id = json_escape(rec.get("occurrenceID"));
        let sci = json_escape(rec.get("scientificName"));
        let acc = json_escape(rec.get("acceptedNameUsage"));
        let basis = json_escape(rec.get("basisOfRecord"));
        let date = json_escape(rec.get("eventDate"));
        let loc = json_escape(rec.get("locality"));
        let status = if rec.has_errors {
            "error"
        } else if rec.has_warnings {
            "warning"
        } else {
            "clean"
        };

        let _ = write!(
            out,
            "    {{\n      \"type\": \"Feature\",\n      \"geometry\": {{ \"type\": \"Point\", \"coordinates\": [{lon:.6}, {lat:.6}] }},\n      \"properties\": {{ \"occurrenceID\": \"{id}\", \"scientificName\": \"{sci}\", \"acceptedNameUsage\": \"{acc}\", \"basisOfRecord\": \"{basis}\", \"eventDate\": \"{date}\", \"locality\": \"{loc}\", \"qaStatus\": \"{status}\" }}\n    }}"
        );
    }

    out.push_str("\n  ]\n}\n");
    out
}

/// Exports the full diagnostic issue ledger as a tab-separated audit report.
pub fn export_issues_tsv(dataset: &DwcDataset) -> String {
    let mut out = String::from(
        "row\tseverity\tcode\tcolumn\tauto_applied\toriginal_value\tsuggested_value\tmessage\n",
    );
    for iss in &dataset.issues {
        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            iss.row_number_1,
            iss.severity.label(),
            iss.code.as_str(),
            iss.column_term,
            iss.auto_applied,
            iss.original_value.replace('\t', " "),
            iss.suggested_value
                .as_deref()
                .unwrap_or("")
                .replace('\t', " "),
            iss.message.replace('\t', " ")
        );
    }
    out
}

/// Renders a standalone vector `.svg` spatial scatter map of georeferenced Darwin Core occurrences
/// with automatic WGS84 bounding-box graticules and taxon color coding.
pub fn render_dwc_occurrence_map_svg(dataset: &DwcDataset) -> String {
    let width = 820.0_f64;
    let height = 520.0_f64;
    let pad_left = 74.0_f64;
    let pad_right = 240.0_f64;
    let pad_top = 58.0_f64;
    let pad_bot = 56.0_f64;
    let plot_w = width - pad_left - pad_right;
    let plot_h = height - pad_top - pad_bot;

    let pts: Vec<(f64, f64, String, bool)> = dataset
        .records
        .iter()
        .filter_map(|r| {
            r.valid_lat_lon.map(|(la, lo)| {
                let taxon = if !r.get("acceptedNameUsage").is_empty() {
                    r.get("acceptedNameUsage").to_string()
                } else {
                    r.get("scientificName").to_string()
                };
                (la, lo, taxon, r.has_warnings || r.has_errors)
            })
        })
        .collect();

    let (min_lat, max_lat, min_lon, max_lon) = if pts.is_empty() {
        (-90.0, 90.0, -180.0, 180.0)
    } else {
        let mut la0 = f64::INFINITY;
        let mut la1 = f64::NEG_INFINITY;
        let mut lo0 = f64::INFINITY;
        let mut lo1 = f64::NEG_INFINITY;
        for &(la, lo, _, _) in &pts {
            la0 = la0.min(la);
            la1 = la1.max(la);
            lo0 = lo0.min(lo);
            lo1 = lo1.max(lo);
        }
        let d_la = ((la1 - la0) * 0.15).max(0.25);
        let d_lo = ((lo1 - lo0) * 0.15).max(0.25);
        (
            (la0 - d_la).clamp(-90.0, 90.0),
            (la1 + d_la).clamp(-90.0, 90.0),
            (lo0 - d_lo).clamp(-180.0, 180.0),
            (lo1 + d_lo).clamp(-180.0, 180.0),
        )
    };

    const PALETTE: &[&str] = &[
        "#0284c7", "#16a34a", "#d97706", "#9333ea", "#e11d48", "#0d9488", "#4f46e5", "#ca8a04",
    ];
    let mut taxon_colors: BTreeMap<String, &'static str> = BTreeMap::new();
    for (_, _, t, _) in &pts {
        let next_idx = taxon_colors.len() % PALETTE.len();
        taxon_colors.entry(t.clone()).or_insert(PALETTE[next_idx]);
    }

    let mut svg = begin_svg_canvas(8192, width, height, Some("#f8fafc"));
    let _ = writeln!(
        svg,
        "  <text x=\"24\" y=\"32\" font-family=\"sans-serif\" font-size=\"16\" font-weight=\"700\" fill=\"#0f172a\">Darwin Core Spatial QA Map — {} ({} valid points)</text>",
        xml_escape(&dataset.dataset_name),
        pts.len()
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{pad_left:.1}\" y=\"{pad_top:.1}\" width=\"{plot_w:.1}\" height=\"{plot_h:.1}\" fill=\"#ecfeff\" stroke=\"#94a3b8\" stroke-width=\"1.2\"/>"
    );

    // 4x4 Graticule lines & coordinate labels
    for step in 0..=4 {
        let frac = (step as f64) / 4.0;
        let x = pad_left + frac * plot_w;
        let y = pad_top + frac * plot_h;
        let lon_v = min_lon + frac * (max_lon - min_lon);
        let lat_v = max_lat - frac * (max_lat - min_lat);

        let _ = writeln!(
            svg,
            "  <line x1=\"{x:.1}\" y1=\"{pad_top:.1}\" x2=\"{x:.1}\" y2=\"{:.1}\" stroke=\"#cbd5e1\" stroke-dasharray=\"3,3\"/>",
            pad_top + plot_h
        );
        let _ = writeln!(
            svg,
            "  <line x1=\"{pad_left:.1}\" y1=\"{y:.1}\" x2=\"{:.1}\" y2=\"{y:.1}\" stroke=\"#cbd5e1\" stroke-dasharray=\"3,3\"/>",
            pad_left + plot_w
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{x:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#334155\">{lon_v:.2}°</text>",
            pad_top + plot_h + 18.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#334155\">{lat_v:.2}°</text>",
            pad_left - 8.0,
            y + 4.0
        );
    }

    // Plot occurrence points
    let lon_span = (max_lon - min_lon).max(1e-6);
    let lat_span = (max_lat - min_lat).max(1e-6);
    for (la, lo, taxon, flagged) in &pts {
        let px = pad_left + ((*lo - min_lon) / lon_span) * plot_w;
        let py = pad_top + ((max_lat - *la) / lat_span) * plot_h;
        let fill = taxon_colors.get(taxon).copied().unwrap_or("#0284c7");
        let stroke = if *flagged { "#dc2626" } else { "#ffffff" };
        let sw = if *flagged { 2.2 } else { 1.3 };
        let _ = writeln!(
            svg,
            "  <circle cx=\"{px:.2}\" cy=\"{py:.2}\" r=\"5.5\" fill=\"{fill}\" stroke=\"{stroke}\" stroke-width=\"{sw:.1}\"/>"
        );
    }

    // Taxon legend on right sidebar
    let leg_x = pad_left + plot_w + 18.0;
    let _ = writeln!(
        svg,
        "  <text x=\"{leg_x:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"12.5\" font-weight=\"700\" fill=\"#0f172a\">Accepted Taxa</text>",
        pad_top + 14.0
    );
    for (i, (taxon, col)) in taxon_colors.iter().take(12).enumerate() {
        let ly = pad_top + 36.0 + (i as f64) * 22.0;
        let _ = writeln!(
            svg,
            "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"5\" fill=\"{col}\"/>",
            leg_x + 6.0,
            ly - 4.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{ly:.1}\" font-family=\"sans-serif\" font-style=\"italic\" font-size=\"11\" fill=\"#1e293b\">{}</text>",
            leg_x + 18.0,
            xml_escape(taxon)
        );
    }

    let _ = writeln!(
        svg,
        "  <text x=\"24\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"10.5\" fill=\"#64748b\">WGS84 Decimal Degrees • Red ring = QA warning/outlier • OSABT GPL-3.0 (CRV Science Network)</text>",
        height - 14.0
    );
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dwc_core::{lint_and_clean_dwc, DwcCleanConfig};
    use crate::dwc_taxonomy::TaxonBackbone;

    #[test]
    fn exports_valid_dwca_meta_xml_geojson_and_svg_map() {
        let bb = TaxonBackbone::default_curated();
        let ds = lint_and_clean_dwc(
            "Oregon_Survey",
            DwcDataset::demo_pnw_field_survey_csv(),
            &bb,
            &DwcCleanConfig::default(),
        )
        .unwrap();

        let meta = generate_dwca_meta_xml(&ds, "occurrence.csv", Some("eml.xml"));
        assert!(meta.contains("xmlns=\"http://rs.tdwg.org/dwc/text/\""));
        assert!(meta.contains("http://rs.tdwg.org/dwc/terms/scientificName"));

        let clean_csv = write_dwc_delimited(&ds, ',', true);
        // Row 7 (Null Island error) excluded -> 7 clean/warning rows + 1 header = 8 lines
        assert_eq!(clean_csv.lines().count(), 8);

        let geojson = export_dwc_geojson(&ds);
        assert!(geojson.contains("\"FeatureCollection\""));

        let svg = render_dwc_occurrence_map_svg(&ds);
        assert!(svg.starts_with("<svg"));
        assert!(svg.trim_end().ends_with("</svg>"));
    }
}
