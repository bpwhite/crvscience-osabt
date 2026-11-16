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

//! Tool #3 (`osabt dwc`) — Scientific binomial/trinomial parser, authorship splitter,
//! Damerau–Levenshtein fuzzy typo detector, offline taxonomic synonym resolver, custom
//! GBIF/Catalogue-of-Life checklist `.tsv` loader, and GBIF v1 Species Match API helper.

use std::collections::HashMap;

use crate::dwc_core::{detect_delimiter, parse_delimited_records};

/// Parsed structure of a scientific name string (ICZN / ICNafp).
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ParsedScientificName {
    pub raw_input: String,
    pub canonical_name: String,
    pub genus_or_uninomial: String,
    pub specific_epithet: Option<String>,
    pub infraspecific_epithet: Option<String>,
    pub authorship: String,
    pub inferred_rank: &'static str,
    pub is_hybrid: bool,
    pub open_nomenclature_qualifier: Option<String>,
}

/// Taxonomic backbone record supporting both accepted taxa and heterotypic/homotypic synonyms.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BackboneTaxon {
    pub query_canonical: String,
    pub accepted_name: String,
    pub authorship: String,
    pub taxonomic_status: String,
    pub rank: String,
    pub gbif_usage_key: u64,
    pub kingdom: String,
    pub phylum: String,
    pub class_name: String,
    pub order: String,
    pub family: String,
    pub genus: String,
}

/// Classification of how a scientific name matched the taxonomic backbone.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaxonMatchKind {
    ExactAccepted,
    Synonym,
    FuzzyOrthographic { edit_distance: usize },
    Unmatched,
}

/// Complete resolution result for a `scientificName` cell.
#[derive(Debug, Clone)]
pub struct TaxonResolution {
    pub parsed: ParsedScientificName,
    pub match_kind: TaxonMatchKind,
    pub matched_entry: Option<BackboneTaxon>,
    pub gbif_lookup_url: String,
}

/// In-memory taxonomic checklist & synonym lookup index.
#[derive(Debug, Clone)]
pub struct TaxonBackbone {
    entries_by_lower: HashMap<String, BackboneTaxon>,
}

impl TaxonBackbone {
    pub fn empty() -> Self {
        Self {
            entries_by_lower: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.entries_by_lower.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries_by_lower.is_empty()
    }

    pub fn insert(&mut self, entry: BackboneTaxon) {
        self.entries_by_lower
            .insert(entry.query_canonical.to_ascii_lowercase(), entry);
    }

    /// Loads additional taxonomic records from a local CSV/TSV checklist file with headers such as:
    /// `scientificName,acceptedNameUsage,taxonomicStatus,taxonRank,taxonID,kingdom,phylum,class,order,family,genus`.
    pub fn load_checklist_delimited(&mut self, text: &str) -> Result<usize, String> {
        let delim = detect_delimiter(text);
        let rows = parse_delimited_records(text, delim);
        if rows.len() < 2 {
            return Err(String::from("Checklist file must contain a header and >= 1 taxon row."));
        }
        let headers: Vec<String> = rows[0]
            .iter()
            .map(|h| {
                h.chars()
                    .filter(|c| c.is_ascii_alphanumeric())
                    .map(|c| c.to_ascii_lowercase())
                    .collect()
            })
            .collect();

        let col_idx = |names: &[&str]| -> Option<usize> {
            headers
                .iter()
                .position(|h| names.iter().any(|n| h == n))
        };

        let sci_col = col_idx(&["scientificname", "canonicalname", "taxon", "name"])
            .ok_or_else(|| String::from("Checklist missing 'scientificName' or 'canonicalName' column."))?;
        let acc_col = col_idx(&["acceptednameusage", "acceptedname", "acceptedscientificname"]);
        let stat_col = col_idx(&["taxonomicstatus", "status"]);
        let rank_col = col_idx(&["taxonrank", "rank"]);
        let key_col = col_idx(&["taxonid", "usagekey", "gbifkey", "id"]);
        let k_col = col_idx(&["kingdom"]);
        let p_col = col_idx(&["phylum"]);
        let c_col = col_idx(&["class"]);
        let o_col = col_idx(&["order"]);
        let f_col = col_idx(&["family"]);
        let g_col = col_idx(&["genus"]);

        let mut added = 0usize;
        for row in rows.iter().skip(1) {
            let raw_sci = row.get(sci_col).cloned().unwrap_or_default();
            if raw_sci.is_empty() {
                continue;
            }
            let parsed = parse_scientific_name(&raw_sci);
            let accepted_raw = acc_col
                .and_then(|i| row.get(i))
                .filter(|s| !s.trim().is_empty())
                .map(|s| parse_scientific_name(s).canonical_name)
                .unwrap_or_else(|| parsed.canonical_name.clone());
            let status = stat_col
                .and_then(|i| row.get(i))
                .filter(|s| !s.trim().is_empty())
                .map(|s| s.to_ascii_uppercase())
                .unwrap_or_else(|| {
                    if accepted_raw.eq_ignore_ascii_case(&parsed.canonical_name) {
                        String::from("ACCEPTED")
                    } else {
                        String::from("SYNONYM")
                    }
                });
            let usage_key = key_col
                .and_then(|i| row.get(i))
                .and_then(|s| {
                    s.trim_start_matches("GBIF:")
                        .trim()
                        .parse::<u64>()
                        .ok()
                })
                .unwrap_or(0);

            self.insert(BackboneTaxon {
                query_canonical: parsed.canonical_name.clone(),
                accepted_name: accepted_raw,
                authorship: parsed.authorship,
                taxonomic_status: status,
                rank: rank_col
                    .and_then(|i| row.get(i))
                    .cloned()
                    .filter(|s| !s.is_empty())
                    .unwrap_or_else(|| parsed.inferred_rank.to_string()),
                gbif_usage_key: usage_key,
                kingdom: k_col.and_then(|i| row.get(i)).cloned().unwrap_or_default(),
                phylum: p_col.and_then(|i| row.get(i)).cloned().unwrap_or_default(),
                class_name: c_col.and_then(|i| row.get(i)).cloned().unwrap_or_default(),
                order: o_col.and_then(|i| row.get(i)).cloned().unwrap_or_default(),
                family: f_col.and_then(|i| row.get(i)).cloned().unwrap_or_default(),
                genus: g_col
                    .and_then(|i| row.get(i))
                    .cloned()
                    .unwrap_or_else(|| parsed.genus_or_uninomial.clone()),
            });
            added += 1;
        }
        Ok(added)
    }

    /// Built-in curated offline taxonomic reference covering Pacific Northwest indicator flora,
    /// fauna, fungi, and standard biological model organisms plus classic historical synonyms.
    pub fn default_curated() -> Self {
        let mut bb = Self::empty();
        // (query, accepted, author, status, rank, gbif_key, kingdom, phylum, class, order, family, genus)
        let rows: &[(&str, &str, &str, &str, &str, u64, &str, &str, &str, &str, &str, &str)] = &[
            // Conifers & vascular plants + classic synonyms
            ("Pseudotsuga menziesii", "Pseudotsuga menziesii", "(Mirb.) Franco", "ACCEPTED", "species", 2685796, "Plantae", "Tracheophyta", "Pinopsida", "Pinales", "Pinaceae", "Pseudotsuga"),
            ("Pseudotsuga taxifolia", "Pseudotsuga menziesii", "(Mirb.) Franco", "SYNONYM", "species", 2685796, "Plantae", "Tracheophyta", "Pinopsida", "Pinales", "Pinaceae", "Pseudotsuga"),
            ("Thuja plicata", "Thuja plicata", "Donn ex D.Don", "ACCEPTED", "species", 2684174, "Plantae", "Tracheophyta", "Pinopsida", "Pinales", "Cupressaceae", "Thuja"),
            ("Tsuga heterophylla", "Tsuga heterophylla", "(Raf.) Sarg.", "ACCEPTED", "species", 2687169, "Plantae", "Tracheophyta", "Pinopsida", "Pinales", "Pinaceae", "Tsuga"),
            ("Alnus rubra", "Alnus rubra", "Bong.", "ACCEPTED", "species", 2876213, "Plantae", "Tracheophyta", "Magnoliopsida", "Fagales", "Betulaceae", "Alnus"),
            ("Acer macrophyllum", "Acer macrophyllum", "Pursh", "ACCEPTED", "species", 3189867, "Plantae", "Tracheophyta", "Magnoliopsida", "Sapindales", "Sapindaceae", "Acer"),
            ("Quercus garryana", "Quercus garryana", "Douglas ex Hook.", "ACCEPTED", "species", 2878373, "Plantae", "Tracheophyta", "Magnoliopsida", "Fagales", "Fagaceae", "Quercus"),
            ("Polystichum munitum", "Polystichum munitum", "(Kaulf.) C.Presl", "ACCEPTED", "species", 2651108, "Plantae", "Tracheophyta", "Polypodiopsida", "Polypodiales", "Dryopteridaceae", "Polystichum"),
            ("Mahonia aquifolium", "Berberis aquifolium", "Pursh", "SYNONYM", "species", 3033897, "Plantae", "Tracheophyta", "Magnoliopsida", "Ranunculales", "Berberidaceae", "Berberis"),
            ("Berberis aquifolium", "Berberis aquifolium", "Pursh", "ACCEPTED", "species", 3033897, "Plantae", "Tracheophyta", "Magnoliopsida", "Ranunculales", "Berberidaceae", "Berberis"),
            ("Trillium ovatum", "Trillium ovatum", "Pursh", "ACCEPTED", "species", 2742264, "Plantae", "Tracheophyta", "Liliopsida", "Liliales", "Melanthiaceae", "Trillium"),
            ("Camassia quamash", "Camassia quamash", "(Pursh) Greene", "ACCEPTED", "species", 2774247, "Plantae", "Tracheophyta", "Liliopsida", "Asparagales", "Asparagaceae", "Camassia"),
            ("Arabidopsis thaliana", "Arabidopsis thaliana", "(L.) Heynh.", "ACCEPTED", "species", 3042845, "Plantae", "Tracheophyta", "Magnoliopsida", "Brassicales", "Brassicaceae", "Arabidopsis"),
            ("Zea mays", "Zea mays", "L.", "ACCEPTED", "species", 5290052, "Plantae", "Tracheophyta", "Liliopsida", "Poales", "Poaceae", "Zea"),
            // Birds, Amphibians, Salmonids & Invertebrates + AOU/herp synonyms
            ("Setophaga coronata", "Setophaga coronata", "(Linnaeus, 1766)", "ACCEPTED", "species", 7341802, "Animalia", "Chordata", "Aves", "Passeriformes", "Parulidae", "Setophaga"),
            ("Dendroica coronata", "Setophaga coronata", "(Linnaeus, 1766)", "SYNONYM", "species", 7341802, "Animalia", "Chordata", "Aves", "Passeriformes", "Parulidae", "Setophaga"),
            ("Poecile atricapillus", "Poecile atricapillus", "(Linnaeus, 1766)", "ACCEPTED", "species", 2487879, "Animalia", "Chordata", "Aves", "Passeriformes", "Paridae", "Poecile"),
            ("Parus atricapillus", "Poecile atricapillus", "(Linnaeus, 1766)", "SYNONYM", "species", "2487879".parse().unwrap_or(2487879), "Animalia", "Chordata", "Aves", "Passeriformes", "Paridae", "Poecile"),
            ("Strix occidentalis", "Strix occidentalis", "(Xantus de Vesey, 1860)", "ACCEPTED", "species", 2497255, "Animalia", "Chordata", "Aves", "Strigiformes", "Strigidae", "Strix"),
            ("Rana cascadae", "Rana cascadae", "Slater, 1939", "ACCEPTED", "species", 2426797, "Animalia", "Chordata", "Amphibia", "Anura", "Ranidae", "Rana"),
            ("Rana aurora", "Rana aurora", "Baird & Girard, 1852", "ACCEPTED", "species", 2426819, "Animalia", "Chordata", "Amphibia", "Anura", "Ranidae", "Rana"),
            ("Anaxyrus boreas", "Anaxyrus boreas", "(Baird & Girard, 1852)", "ACCEPTED", "species", 2422495, "Animalia", "Chordata", "Amphibia", "Anura", "Bufonidae", "Anaxyrus"),
            ("Bufo boreas", "Anaxyrus boreas", "(Baird & Girard, 1852)", "SYNONYM", "species", 2422495, "Animalia", "Chordata", "Amphibia", "Anura", "Bufonidae", "Anaxyrus"),
            ("Taricha granulosa", "Taricha granulosa", "(Skilton, 1849)", "ACCEPTED", "species", 2431886, "Animalia", "Chordata", "Amphibia", "Caudata", "Salamandridae", "Taricha"),
            ("Oncorhynchus tshawytscha", "Oncorhynchus tshawytscha", "(Walbaum, 1792)", "ACCEPTED", "species", 2351435, "Animalia", "Chordata", "Actinopterygii", "Salmoniformes", "Salmonidae", "Oncorhynchus"),
            ("Oncorhynchus kisutch", "Oncorhynchus kisutch", "(Walbaum, 1792)", "ACCEPTED", "species", 2351477, "Animalia", "Chordata", "Actinopterygii", "Salmoniformes", "Salmonidae", "Oncorhynchus"),
            ("Oncorhynchus mykiss", "Oncorhynchus mykiss", "(Walbaum, 1792)", "ACCEPTED", "species", 2351481, "Animalia", "Chordata", "Actinopterygii", "Salmoniformes", "Salmonidae", "Oncorhynchus"),
            ("Salmo gairdneri", "Oncorhynchus mykiss", "(Walbaum, 1792)", "SYNONYM", "species", 2351481, "Animalia", "Chordata", "Actinopterygii", "Salmoniformes", "Salmonidae", "Oncorhynchus"),
            ("Drosophila melanogaster", "Drosophila melanogaster", "Meigen, 1830", "ACCEPTED", "species", 5073713, "Animalia", "Arthropoda", "Insecta", "Diptera", "Drosophilidae", "Drosophila"),
            ("Caenorhabditis elegans", "Caenorhabditis elegans", "(Maupas, 1900)", "ACCEPTED", "species", 2283683, "Animalia", "Nematoda", "Chromadorea", "Rhabditida", "Rhabditidae", "Caenorhabditis"),
            ("Mus musculus", "Mus musculus", "Linnaeus, 1758", "ACCEPTED", "species", 2438793, "Animalia", "Chordata", "Mammalia", "Rodentia", "Muridae", "Mus"),
            ("Danio rerio", "Danio rerio", "(Hamilton, 1822)", "ACCEPTED", "species", 2362651, "Animalia", "Chordata", "Actinopterygii", "Cypriniformes", "Danionidae", "Danio"),
            // Fungi & Microbes
            ("Cantharellus formosus", "Cantharellus formosus", "Corner", "ACCEPTED", "species", 5249505, "Fungi", "Basidiomycota", "Agaricomycetes", "Cantharellales", "Cantharellaceae", "Cantharellus"),
            ("Cantharellus cibarius", "Cantharellus cibarius", "Fr.", "ACCEPTED", "species", 5249504, "Fungi", "Basidiomycota", "Agaricomycetes", "Cantharellales", "Cantharellaceae", "Cantharellus"),
            ("Saccharomyces cerevisiae", "Saccharomyces cerevisiae", "Meyen ex E.C.Hansen", "ACCEPTED", "species", 2599166, "Fungi", "Ascomycota", "Saccharomycetes", "Saccharomycetales", "Saccharomycetaceae", "Saccharomyces"),
            ("Escherichia coli", "Escherichia coli", "(Migula, 1895) Castellani & Chalmers, 1919", "ACCEPTED", "species", 3221936, "Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia"),
            ("Bacillus subtilis", "Bacillus subtilis", "(Ehrenberg, 1835) Cohn, 1872", "ACCEPTED", "species", 5428117, "Bacteria", "Bacillota", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus"),
        ];

        for &(q, acc, auth, stat, rank, key, k, p, c, o, f, g) in rows {
            bb.insert(BackboneTaxon {
                query_canonical: q.to_string(),
                accepted_name: acc.to_string(),
                authorship: auth.to_string(),
                taxonomic_status: stat.to_string(),
                rank: rank.to_string(),
                gbif_usage_key: key,
                kingdom: k.to_string(),
                phylum: p.to_string(),
                class_name: c.to_string(),
                order: o.to_string(),
                family: f.to_string(),
                genus: g.to_string(),
            });
        }
        bb
    }
}

/// Parses a scientific name into canonical binomial/trinomial components and authorship.
pub fn parse_scientific_name(raw: &str) -> ParsedScientificName {
    let clean = raw.trim();
    let mut is_hybrid = false;
    let mut qualifier: Option<String> = None;

    let raw_tokens: Vec<&str> = clean.split_whitespace().collect();
    if raw_tokens.is_empty() {
        return ParsedScientificName {
            raw_input: String::new(),
            canonical_name: String::new(),
            genus_or_uninomial: String::new(),
            specific_epithet: None,
            infraspecific_epithet: None,
            authorship: String::new(),
            inferred_rank: "unranked",
            is_hybrid: false,
            open_nomenclature_qualifier: None,
        };
    }

    let mut idx = 0usize;
    if raw_tokens[idx] == "×" || raw_tokens[idx] == "x" || raw_tokens[idx] == "X" {
        is_hybrid = true;
        idx += 1;
    }

    let genus_raw = raw_tokens.get(idx).copied().unwrap_or("");
    let genus = titlecase_word(genus_raw.trim_start_matches('×'));
    if genus_raw.starts_with('×') {
        is_hybrid = true;
    }
    idx += 1;

    // Optional open-nomenclature qualifier (`cf.`, `aff.`, `sp.`, `spp.`)
    if let Some(&tok) = raw_tokens.get(idx) {
        let low = tok.to_ascii_lowercase();
        if matches!(low.as_str(), "cf." | "cf" | "aff." | "aff") {
            qualifier = Some(tok.to_string());
            idx += 1;
        } else if matches!(low.as_str(), "sp." | "sp" | "spp." | "spp") {
            return ParsedScientificName {
                raw_input: clean.to_string(),
                canonical_name: genus.clone(),
                genus_or_uninomial: genus,
                specific_epithet: None,
                infraspecific_epithet: None,
                authorship: raw_tokens[(idx + 1)..].join(" "),
                inferred_rank: "genus",
                is_hybrid,
                open_nomenclature_qualifier: Some(tok.to_string()),
            };
        }
    }

    if let Some(&tok) = raw_tokens.get(idx) {
        if tok == "×" || tok == "x" {
            is_hybrid = true;
            idx += 1;
        }
    }

    let mut specific: Option<String> = None;
    let mut infraspecific: Option<String> = None;
    let mut rank: &'static str = "uninomial";

    if let Some(&tok) = raw_tokens.get(idx) {
        if looks_like_epithet(tok) {
            specific = Some(tok.to_ascii_lowercase());
            rank = "species";
            idx += 1;
        }
    }

    // Check for infraspecific rank marker (`subsp.`, `ssp.`, `var.`, `f.`) or direct trinomial
    if specific.is_some() {
        if let Some(&tok) = raw_tokens.get(idx) {
            let low = tok.to_ascii_lowercase();
            if matches!(low.as_str(), "subsp." | "subsp" | "ssp." | "ssp") {
                rank = "subspecies";
                idx += 1;
                if let Some(&infra) = raw_tokens.get(idx) {
                    if looks_like_epithet(infra) {
                        infraspecific = Some(infra.to_ascii_lowercase());
                        idx += 1;
                    }
                }
            } else if matches!(low.as_str(), "var." | "var") {
                rank = "variety";
                idx += 1;
                if let Some(&infra) = raw_tokens.get(idx) {
                    if looks_like_epithet(infra) {
                        infraspecific = Some(infra.to_ascii_lowercase());
                        idx += 1;
                    }
                }
            } else if matches!(low.as_str(), "f." | "forma") {
                rank = "form";
                idx += 1;
                if let Some(&infra) = raw_tokens.get(idx) {
                    if looks_like_epithet(infra) {
                        infraspecific = Some(infra.to_ascii_lowercase());
                        idx += 1;
                    }
                }
            } else if looks_like_epithet(tok) {
                rank = "subspecies";
                infraspecific = Some(tok.to_ascii_lowercase());
                idx += 1;
            }
        }
    }

    let authorship = raw_tokens.get(idx..).map(|s| s.join(" ")).unwrap_or_default();
    let mut canonical = genus.clone();
    if let Some(sp) = &specific {
        canonical.push(' ');
        canonical.push_str(sp);
    }
    if let Some(inf) = &infraspecific {
        canonical.push(' ');
        canonical.push_str(inf);
    }

    ParsedScientificName {
        raw_input: clean.to_string(),
        canonical_name: canonical,
        genus_or_uninomial: genus,
        specific_epithet: specific,
        infraspecific_epithet: infraspecific,
        authorship,
        inferred_rank: rank,
        is_hybrid,
        open_nomenclature_qualifier: qualifier,
    }
}

fn looks_like_epithet(token: &str) -> bool {
    let first = match token.chars().next() {
        Some(c) => c,
        None => return false,
    };
    first.is_ascii_lowercase()
        && token
            .chars()
            .all(|c| c.is_ascii_lowercase() || c == '-')
}

fn titlecase_word(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        Some(first) => {
            let mut out = String::with_capacity(s.len());
            out.push(first.to_ascii_uppercase());
            for ch in chars {
                out.push(ch.to_ascii_lowercase());
            }
            out
        }
        None => String::new(),
    }
}

/// Resolves a scientific name against `backbone` (exact, synonym, or Damerau–Levenshtein <= 2).
pub fn resolve_scientific_name(raw_name: &str, backbone: &TaxonBackbone) -> TaxonResolution {
    let parsed = parse_scientific_name(raw_name);
    let key = parsed.canonical_name.to_ascii_lowercase();
    let gbif_url = build_gbif_species_match_url(&parsed.canonical_name);

    if let Some(entry) = backbone.entries_by_lower.get(&key) {
        let kind = if entry.taxonomic_status.eq_ignore_ascii_case("SYNONYM")
            || !entry
                .accepted_name
                .eq_ignore_ascii_case(&parsed.canonical_name)
        {
            TaxonMatchKind::Synonym
        } else {
            TaxonMatchKind::ExactAccepted
        };
        return TaxonResolution {
            parsed,
            match_kind: kind,
            matched_entry: Some(entry.clone()),
            gbif_lookup_url: gbif_url,
        };
    }

    // Fuzzy orthographic check (Damerau–Levenshtein <= 2 on binomials of length >= 8)
    if key.len() >= 8 && parsed.specific_epithet.is_some() {
        let mut best: Option<(usize, &BackboneTaxon)> = None;
        for (cand_key, entry) in &backbone.entries_by_lower {
            if cand_key.len().abs_diff(key.len()) > 2 {
                continue;
            }
            // Require matching initial character to keep taxonomic matching conservative
            if cand_key.as_bytes().first() != key.as_bytes().first() {
                continue;
            }
            let dist = damerau_levenshtein(&key, cand_key);
            if dist <= 2 && best.map(|(b, _)| dist < b).unwrap_or(true) {
                best = Some((dist, entry));
            }
        }
        if let Some((dist, entry)) = best {
            return TaxonResolution {
                parsed,
                match_kind: TaxonMatchKind::FuzzyOrthographic {
                    edit_distance: dist,
                },
                matched_entry: Some(entry.clone()),
                gbif_lookup_url: gbif_url,
            };
        }
    }

    TaxonResolution {
        parsed,
        match_kind: TaxonMatchKind::Unmatched,
        matched_entry: None,
        gbif_lookup_url: gbif_url,
    }
}

/// Optimal string alignment (restricted Damerau–Levenshtein) distance handling adjacent transpositions
/// (e.g. `menzeisii` <-> `menziesii` = edit distance 1).
pub fn damerau_levenshtein(a: &str, b: &str) -> usize {
    let ac: Vec<char> = a.chars().collect();
    let bc: Vec<char> = b.chars().collect();
    let m = ac.len();
    let n = bc.len();
    if m == 0 {
        return n;
    }
    if n == 0 {
        return m;
    }

    let mut dp = vec![vec![0usize; n + 1]; m + 1];
    for (i, row) in dp.iter_mut().enumerate().take(m + 1) {
        row[0] = i;
    }
    for j in 0..=n {
        dp[0][j] = j;
    }

    for i in 1..=m {
        for j in 1..=n {
            let cost = usize::from(ac[i - 1] != bc[j - 1]);
            let mut val = (dp[i - 1][j] + 1)
                .min(dp[i][j - 1] + 1)
                .min(dp[i - 1][j - 1] + cost);
            if i > 1 && j > 1 && ac[i - 1] == bc[j - 2] && ac[i - 2] == bc[j - 1] {
                val = val.min(dp[i - 2][j - 2] + 1);
            }
            dp[i][j] = val;
        }
    }
    dp[m][n]
}

/// Builds a URL-encoded GBIF v1 Species Match endpoint (`https://api.gbif.org/v1/species/match?name=...`).
pub fn build_gbif_species_match_url(canonical_name: &str) -> String {
    let encoded: String = canonical_name
        .bytes()
        .map(|b| match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'-' | b'_' | b'.' => {
                (b as char).to_string()
            }
            b' ' => String::from("%20"),
            other => format!("%{other:02X}"),
        })
        .collect();
    format!("https://api.gbif.org/v1/species/match?verbose=true&name={encoded}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_authorship_and_resolves_synonyms_and_transpositions() {
        let bb = TaxonBackbone::default_curated();

        let syn = resolve_scientific_name("Pseudotsuga taxifolia (Poir.) Britton", &bb);
        assert_eq!(syn.parsed.canonical_name, "Pseudotsuga taxifolia");
        assert_eq!(syn.parsed.authorship, "(Poir.) Britton");
        assert_eq!(syn.match_kind, TaxonMatchKind::Synonym);
        assert_eq!(
            syn.matched_entry.unwrap().accepted_name,
            "Pseudotsuga menziesii"
        );

        let typo = resolve_scientific_name("Pseudotsuga menzeisii", &bb);
        assert_eq!(
            typo.match_kind,
            TaxonMatchKind::FuzzyOrthographic { edit_distance: 1 }
        );
        assert_eq!(
            typo.matched_entry.unwrap().accepted_name,
            "Pseudotsuga menziesii"
        );
    }
}
