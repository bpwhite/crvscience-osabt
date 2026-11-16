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

//! Tool #3 (`osabt dwc`) — Core TDWG Darwin Core (`https://dwc.tdwg.org/terms/`) schema engine,
//! RFC-4180 CSV/TSV parser, ISO-8601 `eventDate` scrubber, and Zizka et al. (2019)
//! `CoordinateCleaner`-style spatial quality-control & DMS coordinate normalizer.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::f64::consts::PI;

use crate::dwc_taxonomy::{resolve_scientific_name, TaxonBackbone, TaxonMatchKind};

/// Severity tier for a Darwin Core lint diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum LintSeverity {
    Info,
    Warning,
    Error,
}

impl LintSeverity {
    pub fn label(self) -> &'static str {
        match self {
            Self::Info => "INFO",
            Self::Warning => "WARN",
            Self::Error => "ERROR",
        }
    }
}

/// Categorized diagnostic code modeled after TDWG Darwin Core QA & `CoordinateCleaner` tests.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DwcIssueCode {
    NonStandardHeaderMapped,
    UnrecognizedColumnHeader,
    MissingRequiredTerm,
    DuplicateOccurrenceId,
    MissingOccurrenceId,
    InvalidBasisOfRecord,
    NormalizedBasisOfRecord,
    InvalidEventDate,
    NormalizedEventDate,
    FutureEventDate,
    ConvertedDmsCoordinate,
    CoordinateOutOfBounds,
    TransposedLatLonSuspected,
    ZeroOrNullIslandCoordinate,
    EqualLatLonCoordinate,
    InstitutionOrCentroidCoordinate,
    LowCoordinatePrecision,
    MissingGeodeticDatum,
    SpeciesSpatialOutlier,
    ElevationOutOfBounds,
    InvalidIndividualCount,
    TaxonSynonymResolved,
    TaxonFuzzyTypoSuggested,
    TaxonAuthorshipExtracted,
    UnmatchedScientificName,
}

impl DwcIssueCode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NonStandardHeaderMapped => "DWC_HEADER_ALIAS",
            Self::UnrecognizedColumnHeader => "DWC_UNKNOWN_COLUMN",
            Self::MissingRequiredTerm => "DWC_MISSING_CORE_TERM",
            Self::DuplicateOccurrenceId => "DWC_DUPLICATE_ID",
            Self::MissingOccurrenceId => "DWC_MISSING_ID",
            Self::InvalidBasisOfRecord => "DWC_INVALID_BASIS",
            Self::NormalizedBasisOfRecord => "DWC_NORMALIZED_BASIS",
            Self::InvalidEventDate => "DWC_INVALID_DATE",
            Self::NormalizedEventDate => "DWC_NORMALIZED_DATE",
            Self::FutureEventDate => "DWC_FUTURE_DATE",
            Self::ConvertedDmsCoordinate => "GEO_DMS_CONVERTED",
            Self::CoordinateOutOfBounds => "GEO_OUT_OF_BOUNDS",
            Self::TransposedLatLonSuspected => "GEO_TRANSPOSED_LAT_LON",
            Self::ZeroOrNullIslandCoordinate => "GEO_NULL_ISLAND_ZERO",
            Self::EqualLatLonCoordinate => "GEO_EQUAL_LAT_LON",
            Self::InstitutionOrCentroidCoordinate => "GEO_CENTROID_OR_INST",
            Self::LowCoordinatePrecision => "GEO_LOW_PRECISION",
            Self::MissingGeodeticDatum => "GEO_MISSING_DATUM",
            Self::SpeciesSpatialOutlier => "GEO_SPECIES_OUTLIER",
            Self::ElevationOutOfBounds => "GEO_ELEVATION_RANGE",
            Self::InvalidIndividualCount => "DWC_INDIVIDUAL_COUNT",
            Self::TaxonSynonymResolved => "TAX_SYNONYM_RESOLVED",
            Self::TaxonFuzzyTypoSuggested => "TAX_FUZZY_TYPO",
            Self::TaxonAuthorshipExtracted => "TAX_AUTHORSHIP_SPLIT",
            Self::UnmatchedScientificName => "TAX_UNVERIFIED_NAME",
        }
    }
}

/// A single row- or header-level Darwin Core validation finding.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DwcLintIssue {
    /// 1-based data row index (`0` denotes a dataset/header-level issue).
    pub row_number_1: usize,
    pub column_term: String,
    pub severity: LintSeverity,
    pub code: DwcIssueCode,
    pub message: String,
    pub original_value: String,
    pub suggested_value: Option<String>,
    pub auto_applied: bool,
}

/// Canonical TDWG Darwin Core term metadata entry (`http://rs.tdwg.org/dwc/terms/`).
#[derive(Debug, Clone, Copy)]
pub struct DwcTermSpec {
    pub term: &'static str,
    pub iri: &'static str,
    pub dwc_class: &'static str,
    pub required_core: bool,
    pub aliases: &'static [&'static str],
}

/// Built-in catalog of core TDWG Darwin Core Occurrence, Event, Location, and Taxon terms.
pub static DWC_TERM_CATALOG: &[DwcTermSpec] = &[
    DwcTermSpec {
        term: "occurrenceID",
        iri: "http://rs.tdwg.org/dwc/terms/occurrenceID",
        dwc_class: "Occurrence",
        required_core: true,
        aliases: &["id", "occ_id", "occurrence_id", "guid", "record_id"],
    },
    DwcTermSpec {
        term: "basisOfRecord",
        iri: "http://rs.tdwg.org/dwc/terms/basisOfRecord",
        dwc_class: "Occurrence",
        required_core: true,
        aliases: &["basis", "basis_of_record", "record_type", "specimen_type"],
    },
    DwcTermSpec {
        term: "scientificName",
        iri: "http://rs.tdwg.org/dwc/terms/scientificName",
        dwc_class: "Taxon",
        required_core: true,
        aliases: &[
            "species",
            "taxon",
            "sci_name",
            "scientific_name",
            "species_name",
            "binomial",
        ],
    },
    DwcTermSpec {
        term: "scientificNameAuthorship",
        iri: "http://rs.tdwg.org/dwc/terms/scientificNameAuthorship",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["author", "authorship", "authority", "taxon_author"],
    },
    DwcTermSpec {
        term: "taxonRank",
        iri: "http://rs.tdwg.org/dwc/terms/taxonRank",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["rank", "taxon_rank"],
    },
    DwcTermSpec {
        term: "taxonomicStatus",
        iri: "http://rs.tdwg.org/dwc/terms/taxonomicStatus",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["taxonomic_status", "name_status"],
    },
    DwcTermSpec {
        term: "acceptedNameUsage",
        iri: "http://rs.tdwg.org/dwc/terms/acceptedNameUsage",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["accepted_name", "accepted_scientific_name"],
    },
    DwcTermSpec {
        term: "taxonID",
        iri: "http://rs.tdwg.org/dwc/terms/taxonID",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["taxon_id", "gbif_key", "usage_key"],
    },
    DwcTermSpec {
        term: "kingdom",
        iri: "http://rs.tdwg.org/dwc/terms/kingdom",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["kingdom_name"],
    },
    DwcTermSpec {
        term: "phylum",
        iri: "http://rs.tdwg.org/dwc/terms/phylum",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["division", "phylum_name"],
    },
    DwcTermSpec {
        term: "class",
        iri: "http://rs.tdwg.org/dwc/terms/class",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["class_name"],
    },
    DwcTermSpec {
        term: "order",
        iri: "http://rs.tdwg.org/dwc/terms/order",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["order_name"],
    },
    DwcTermSpec {
        term: "family",
        iri: "http://rs.tdwg.org/dwc/terms/family",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["family_name"],
    },
    DwcTermSpec {
        term: "genus",
        iri: "http://rs.tdwg.org/dwc/terms/genus",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["genus_name"],
    },
    DwcTermSpec {
        term: "vernacularName",
        iri: "http://rs.tdwg.org/dwc/terms/vernacularName",
        dwc_class: "Taxon",
        required_core: false,
        aliases: &["common_name", "vernacular_name"],
    },
    DwcTermSpec {
        term: "eventDate",
        iri: "http://rs.tdwg.org/dwc/terms/eventDate",
        dwc_class: "Event",
        required_core: true,
        aliases: &[
            "date",
            "collection_date",
            "observed_on",
            "event_date",
            "date_collected",
            "survey_date",
        ],
    },
    DwcTermSpec {
        term: "year",
        iri: "http://rs.tdwg.org/dwc/terms/year",
        dwc_class: "Event",
        required_core: false,
        aliases: &["collect_year"],
    },
    DwcTermSpec {
        term: "month",
        iri: "http://rs.tdwg.org/dwc/terms/month",
        dwc_class: "Event",
        required_core: false,
        aliases: &["collect_month"],
    },
    DwcTermSpec {
        term: "day",
        iri: "http://rs.tdwg.org/dwc/terms/day",
        dwc_class: "Event",
        required_core: false,
        aliases: &["collect_day"],
    },
    DwcTermSpec {
        term: "decimalLatitude",
        iri: "http://rs.tdwg.org/dwc/terms/decimalLatitude",
        dwc_class: "Location",
        required_core: true,
        aliases: &["lat", "latitude", "decimallatitude", "lat_dd", "y"],
    },
    DwcTermSpec {
        term: "decimalLongitude",
        iri: "http://rs.tdwg.org/dwc/terms/decimalLongitude",
        dwc_class: "Location",
        required_core: true,
        aliases: &[
            "lon",
            "lng",
            "long",
            "longitude",
            "decimallongitude",
            "lon_dd",
            "x",
        ],
    },
    DwcTermSpec {
        term: "geodeticDatum",
        iri: "http://rs.tdwg.org/dwc/terms/geodeticDatum",
        dwc_class: "Location",
        required_core: false,
        aliases: &["datum", "crs", "srs", "geodetic_datum"],
    },
    DwcTermSpec {
        term: "coordinateUncertaintyInMeters",
        iri: "http://rs.tdwg.org/dwc/terms/coordinateUncertaintyInMeters",
        dwc_class: "Location",
        required_core: false,
        aliases: &[
            "uncertainty_m",
            "accuracy_m",
            "positional_accuracy",
            "coord_uncertainty",
        ],
    },
    DwcTermSpec {
        term: "countryCode",
        iri: "http://rs.tdwg.org/dwc/terms/countryCode",
        dwc_class: "Location",
        required_core: false,
        aliases: &["country", "iso_country", "country_code"],
    },
    DwcTermSpec {
        term: "stateProvince",
        iri: "http://rs.tdwg.org/dwc/terms/stateProvince",
        dwc_class: "Location",
        required_core: false,
        aliases: &["state", "province", "state_province"],
    },
    DwcTermSpec {
        term: "locality",
        iri: "http://rs.tdwg.org/dwc/terms/locality",
        dwc_class: "Location",
        required_core: false,
        aliases: &["site", "site_name", "location", "place_guess"],
    },
    DwcTermSpec {
        term: "minimumElevationInMeters",
        iri: "http://rs.tdwg.org/dwc/terms/minimumElevationInMeters",
        dwc_class: "Location",
        required_core: false,
        aliases: &["elevation", "elev_m", "elevation_m", "altitude_m"],
    },
    DwcTermSpec {
        term: "recordedBy",
        iri: "http://rs.tdwg.org/dwc/terms/recordedBy",
        dwc_class: "Occurrence",
        required_core: false,
        aliases: &["collector", "observer", "recorded_by", "collected_by"],
    },
    DwcTermSpec {
        term: "individualCount",
        iri: "http://rs.tdwg.org/dwc/terms/individualCount",
        dwc_class: "Occurrence",
        required_core: false,
        aliases: &["count", "abundance", "individual_count", "n_individuals"],
    },
    DwcTermSpec {
        term: "occurrenceStatus",
        iri: "http://rs.tdwg.org/dwc/terms/occurrenceStatus",
        dwc_class: "Occurrence",
        required_core: false,
        aliases: &["presence", "occurrence_status"],
    },
    DwcTermSpec {
        term: "institutionCode",
        iri: "http://rs.tdwg.org/dwc/terms/institutionCode",
        dwc_class: "RecordLevel",
        required_core: false,
        aliases: &["institution", "institution_code", "herbarium", "museum"],
    },
    DwcTermSpec {
        term: "catalogNumber",
        iri: "http://rs.tdwg.org/dwc/terms/catalogNumber",
        dwc_class: "Occurrence",
        required_core: false,
        aliases: &["catalog_no", "catalog_number", "voucher_id", "accession"],
    },
    DwcTermSpec {
        term: "occurrenceRemarks",
        iri: "http://rs.tdwg.org/dwc/terms/occurrenceRemarks",
        dwc_class: "Occurrence",
        required_core: false,
        aliases: &["notes", "remarks", "comments", "field_notes"],
    },
];

/// Known biodiversity institution headquarters & national/regional centroids
/// frequently misassigned by automated geocoders (`CoordinateCleaner` reference points).
const KNOWN_CENTROIDS_AND_INSTITUTIONS: &[(&str, f64, f64, f64)] = &[
    ("Null Island (0,0)", 0.0, 0.0, 5.0),
    ("GBIF Secretariat (Copenhagen)", 55.7069, 12.5594, 1.5),
    ("Smithsonian NMNH (Washington, DC)", 38.8913, -77.0261, 1.0),
    ("Natural History Museum (London)", 51.4967, -0.1764, 1.0),
    ("Royal Botanic Gardens, Kew", 51.4787, -0.2956, 1.5),
    ("Geographic Center of Contiguous USA (Kansas)", 39.8283, -98.5795, 4.0),
    ("Geographic Center of Canada (Arviat)", 62.3933, -96.8167, 5.0),
    ("Geographic Center of Oregon (Crook County)", 43.9336, -120.5583, 2.5),
];

/// Configuration options for scrubbing and validating a Darwin Core occurrence table.
#[derive(Debug, Clone)]
pub struct DwcCleanConfig {
    /// Automatically apply safe fixes (header mapping, DMS->DD, transposed lat/lon swap, ISO date,
    /// `basisOfRecord` camelCase, `WGS84` datum, missing `occurrenceID` minting, synonym update).
    pub auto_fix: bool,
    /// Replace synonyms in `scientificName` with accepted names while storing the original in
    /// `taxonRemarks` / `previousIdentifications`, or preserve `scientificName` and populate `acceptedNameUsage`.
    pub replace_synonyms_in_place: bool,
    /// Reference year used to flag impossible future `eventDate` records.
    pub max_valid_year: i32,
    /// Minimum decimal places expected on `decimalLatitude`/`decimalLongitude` before warning.
    pub min_coord_decimals: usize,
    /// Flag conspecific spatial outliers exceeding `outlier_mad_multiplier` * MAD (minimum 250 km).
    pub detect_species_outliers: bool,
    /// Prefix used when minting deterministic `occurrenceID`s for blank rows.
    pub mint_id_prefix: String,
}

impl Default for DwcCleanConfig {
    fn default() -> Self {
        Self {
            auto_fix: true,
            replace_synonyms_in_place: false,
            max_valid_year: 2026,
            min_coord_decimals: 3,
            detect_species_outliers: true,
            mint_id_prefix: String::from("urn:catalog:CRV:OSABT"),
        }
    }
}

/// Parsed occurrence row storing Darwin Core term -> cell value mappings.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DwcRecord {
    pub row_number_1: usize,
    pub fields: BTreeMap<String, String>,
    pub valid_lat_lon: Option<(f64, f64)>,
    pub has_errors: bool,
    pub has_warnings: bool,
}

impl DwcRecord {
    pub fn get(&self, term: &str) -> &str {
        self.fields.get(term).map(String::as_str).unwrap_or("")
    }

    pub fn set(&mut self, term: &str, value: impl Into<String>) {
        self.fields.insert(term.to_string(), value.into());
    }
}

/// Complete Darwin Core dataset plus audit report produced by [`lint_and_clean_dwc`].
#[derive(Debug, Clone)]
pub struct DwcDataset {
    pub dataset_name: String,
    pub delimiter: char,
    /// Ordered list of canonical/retained column headers for export.
    pub columns: Vec<String>,
    /// Mapping of original file header -> mapped Darwin Core term (when renamed).
    pub header_mappings: Vec<(String, String)>,
    pub records: Vec<DwcRecord>,
    pub issues: Vec<DwcLintIssue>,
}

/// Summary counts across a linted `DwcDataset`.
#[derive(Debug, Clone, Default)]
pub struct DwcAuditSummary {
    pub total_records: usize,
    pub clean_records: usize,
    pub records_with_warnings: usize,
    pub records_with_errors: usize,
    pub georeferenced_records: usize,
    pub error_count: usize,
    pub warning_count: usize,
    pub info_count: usize,
    pub auto_fixed_count: usize,
    pub unique_taxa: usize,
    pub synonyms_resolved: usize,
}

impl DwcDataset {
    /// Computes aggregate quality-control metrics for CLI/TUI summaries.
    pub fn summary(&self) -> DwcAuditSummary {
        let mut s = DwcAuditSummary {
            total_records: self.records.len(),
            ..Default::default()
        };
        let mut taxa = HashSet::new();
        for rec in &self.records {
            if rec.has_errors {
                s.records_with_errors += 1;
            } else if rec.has_warnings {
                s.records_with_warnings += 1;
            } else {
                s.clean_records += 1;
            }
            if rec.valid_lat_lon.is_some() {
                s.georeferenced_records += 1;
            }
            let sci = rec.get("scientificName").trim();
            if !sci.is_empty() {
                taxa.insert(sci.to_string());
            }
        }
        s.unique_taxa = taxa.len();
        for iss in &self.issues {
            match iss.severity {
                LintSeverity::Error => s.error_count += 1,
                LintSeverity::Warning => s.warning_count += 1,
                LintSeverity::Info => s.info_count += 1,
            }
            if iss.auto_applied {
                s.auto_fixed_count += 1;
            }
            if iss.code == DwcIssueCode::TaxonSynonymResolved {
                s.synonyms_resolved += 1;
            }
        }
        s
    }

    /// Built-in Pacific Northwest (Oregon) field survey dataset demonstrating realistic
    /// collector sheet pitfalls: non-standard headers, DMS strings, swapped lat/lon, Null Island,
    /// taxonomic synonyms (*Pseudotsuga taxifolia*, *Berberis aquifolium*, *Dendroica coronata*),
    /// slash dates, and missing `occurrenceID`s.
    pub fn demo_pnw_field_survey_csv() -> &'static str {
        "\
occ_id,basis,species,date_collected,lat,lon,datum,uncertainty_m,elev_m,locality,collector,count
CRV-OR-001,human observation,\"Pseudotsuga taxifolia (Poir.) Britton\",07/14/2026,\"44° 30' 14.4\"\" N\",\"123° 33' 03.6\"\" W\",,15,1248,\"Marys Peak Meadow, Benton Co., OR\",B. P. White,4
CRV-OR-002,PreservedSpecimen,Polystichum munitum,2026-07-14,44.5646,-123.2620,WGS84,8,86,\"McDonald-Dunn Forest, Corvallis, OR\",B. P. White,12
,HumanObservation,Berberis aquifolium,2026/07/15,-123.2794,44.5582,WGS84,25,92,\"Bald Hill Natural Area, Corvallis, OR\",B. P. White,3
CRV-OR-004,MachineObservation,Dendroica coronata,2026-05-19,44.6215,-123.2351,WGS84,30,68,\"William L. Finley NWR, OR\",CRV Acoustic Node 2,1
CRV-OR-005,HumanObservation,Cantharellus formosus,2026-10-22,44.41, -123.58,WGS84,,410,\"Alsea Falls Trail, OR\",B. P. White,9
CRV-OR-006,PreservedSpecimen,Rana cascadae,2026-08-03,44.3984,-121.8492,WGS84,12,1540,\"Three Sisters Wilderness, OR\",B. P. White,2
CRV-OR-007,HumanObservation,Pseudotsuga menzeisii,2026-08-04,0.0,0.0,WGS84,,0,\"Unlogged GPS Fix\",Field Intern,1
CRV-OR-002,HumanObservation,Trillium ovatum,2026-04-11,44.5712,-123.2745,WGS84,10,140,\"Oak Creek Trail, Corvallis, OR\",B. P. White,5
"
    }
}

/// Looks up the canonical TDWG Darwin Core term spec for a raw spreadsheet column header.
/// Returns `(Some(spec), exact_match)` when recognized either verbatim or via common field aliases.
pub fn match_dwc_term(raw_header: &str) -> Option<(&'static DwcTermSpec, bool)> {
    let trimmed = raw_header.trim().trim_start_matches("dwc:");
    for spec in DWC_TERM_CATALOG {
        if spec.term == trimmed {
            return Some((spec, true));
        }
    }
    let norm: String = trimmed
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect();
    for spec in DWC_TERM_CATALOG {
        let spec_norm: String = spec
            .term
            .chars()
            .map(|c| c.to_ascii_lowercase())
            .collect();
        if norm == spec_norm {
            return Some((spec, false));
        }
        for alias in spec.aliases {
            let alias_norm: String = alias
                .chars()
                .filter(|c| c.is_ascii_alphanumeric())
                .map(|c| c.to_ascii_lowercase())
                .collect();
            if norm == alias_norm {
                return Some((spec, false));
            }
        }
    }
    None
}

/// Returns the TDWG canonical IRI (`http://rs.tdwg.org/dwc/terms/...`) for a standard term.
pub fn dwc_term_iri(term: &str) -> Option<&'static str> {
    DWC_TERM_CATALOG
        .iter()
        .find(|s| s.term == term)
        .map(|s| s.iri)
}

/// Auto-detects the field delimiter (`,`, `\t`, `;`, or `|`) from the header line.
pub fn detect_delimiter(text: &str) -> char {
    let first_line = text
        .trim_start_matches('\u{feff}')
        .lines()
        .next()
        .unwrap_or("");
    let candidates = [('\t', 0usize), (',', 0usize), (';', 0usize), ('|', 0usize)];
    let mut best = (',', 0usize);
    for (delim, _) in candidates {
        let count = split_delimited_line(first_line, delim).len();
        if count > best.1 {
            best = (delim, count);
        }
    }
    best.0
}

/// Parses an RFC-4180 delimited text stream (CSV/TSV) supporting quoted newlines and `""` escapes.
pub fn parse_delimited_records(text: &str, delimiter: char) -> Vec<Vec<String>> {
    let clean = text.trim_start_matches('\u{feff}');
    let mut rows: Vec<Vec<String>> = Vec::new();
    let mut current_row: Vec<String> = Vec::new();
    let mut field = String::new();
    let mut in_quotes = false;
    let mut chars = clean.chars().peekable();

    while let Some(ch) = chars.next() {
        if in_quotes {
            if ch == '"' {
                if chars.peek() == Some(&'"') {
                    field.push('"');
                    chars.next();
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(ch);
            }
        } else if ch == '"' && field.trim().is_empty() {
            field.clear();
            in_quotes = true;
        } else if ch == delimiter {
            current_row.push(field.trim().to_string());
            field.clear();
        } else if ch == '\r' {
            if chars.peek() == Some(&'\n') {
                chars.next();
            }
            current_row.push(field.trim().to_string());
            field.clear();
            if current_row.iter().any(|c| !c.is_empty()) {
                rows.push(current_row);
            }
            current_row = Vec::new();
        } else if ch == '\n' {
            current_row.push(field.trim().to_string());
            field.clear();
            if current_row.iter().any(|c| !c.is_empty()) {
                rows.push(current_row);
            }
            current_row = Vec::new();
        } else {
            field.push(ch);
        }
    }

    if !field.is_empty() || !current_row.is_empty() {
        current_row.push(field.trim().to_string());
        if current_row.iter().any(|c| !c.is_empty()) {
            rows.push(current_row);
        }
    }

    rows
}

fn split_delimited_line(line: &str, delimiter: char) -> Vec<String> {
    parse_delimited_records(line, delimiter)
        .into_iter()
        .next()
        .unwrap_or_default()
}

/// Parses and validates a Darwin Core occurrence CSV/TSV string against TDWG vocabularies,
/// `CoordinateCleaner` spatial rules, and the supplied [`TaxonBackbone`].
pub fn lint_and_clean_dwc(
    dataset_name: &str,
    raw_text: &str,
    backbone: &TaxonBackbone,
    config: &DwcCleanConfig,
) -> Result<DwcDataset, String> {
    let delimiter = detect_delimiter(raw_text);
    let raw_rows = parse_delimited_records(raw_text, delimiter);
    if raw_rows.is_empty() {
        return Err(String::from(
            "Input Darwin Core file is empty or contains no header row.",
        ));
    }

    let raw_headers = &raw_rows[0];
    let mut issues: Vec<DwcLintIssue> = Vec::new();
    let mut columns: Vec<String> = Vec::with_capacity(raw_headers.len() + 6);
    let mut header_mappings: Vec<(String, String)> = Vec::new();

    for h in raw_headers {
        match match_dwc_term(h) {
            Some((spec, true)) => {
                columns.push(spec.term.to_string());
            }
            Some((spec, false)) => {
                let target = if config.auto_fix {
                    spec.term.to_string()
                } else {
                    h.clone()
                };
                header_mappings.push((h.clone(), spec.term.to_string()));
                issues.push(DwcLintIssue {
                    row_number_1: 0,
                    column_term: spec.term.to_string(),
                    severity: LintSeverity::Info,
                    code: DwcIssueCode::NonStandardHeaderMapped,
                    message: format!(
                        "Mapped non-standard column header '{}' -> canonical Darwin Core term '{}'.",
                        h, spec.term
                    ),
                    original_value: h.clone(),
                    suggested_value: Some(spec.term.to_string()),
                    auto_applied: config.auto_fix,
                });
                columns.push(target);
            }
            None => {
                issues.push(DwcLintIssue {
                    row_number_1: 0,
                    column_term: h.clone(),
                    severity: LintSeverity::Warning,
                    code: DwcIssueCode::UnrecognizedColumnHeader,
                    message: format!(
                        "Column '{}' is not a standard Darwin Core term; retained as custom field.",
                        h
                    ),
                    original_value: h.clone(),
                    suggested_value: None,
                    auto_applied: false,
                });
                columns.push(h.clone());
            }
        }
    }

    // Check core required Darwin Core Occurrence terms
    for spec in DWC_TERM_CATALOG.iter().filter(|s| s.required_core) {
        if !columns.iter().any(|c| c == spec.term) {
            if config.auto_fix && spec.term == "occurrenceID" {
                columns.insert(0, String::from("occurrenceID"));
                issues.push(DwcLintIssue {
                    row_number_1: 0,
                    column_term: String::from("occurrenceID"),
                    severity: LintSeverity::Warning,
                    code: DwcIssueCode::MissingRequiredTerm,
                    message: String::from(
                        "Missing required 'occurrenceID' column; inserted and minting deterministic IDs.",
                    ),
                    original_value: String::new(),
                    suggested_value: Some(String::from("occurrenceID")),
                    auto_applied: true,
                });
            } else {
                issues.push(DwcLintIssue {
                    row_number_1: 0,
                    column_term: spec.term.to_string(),
                    severity: LintSeverity::Error,
                    code: DwcIssueCode::MissingRequiredTerm,
                    message: format!(
                        "Dataset is missing recommended/required Darwin Core column '{}'.",
                        spec.term
                    ),
                    original_value: String::new(),
                    suggested_value: None,
                    auto_applied: false,
                });
            }
        }
    }

    let mut seen_ids: HashMap<String, usize> = HashMap::new();
    let mut records: Vec<DwcRecord> = Vec::with_capacity(raw_rows.len().saturating_sub(1));

    for (idx, cells) in raw_rows.iter().skip(1).enumerate() {
        let row_num = idx + 1;
        let mut fields = BTreeMap::new();
        for (col_idx, raw_col_name) in raw_headers.iter().enumerate() {
            let mapped_col = match match_dwc_term(raw_col_name) {
                Some((spec, _)) if config.auto_fix => spec.term.to_string(),
                _ => raw_col_name.clone(),
            };
            let val = cells.get(col_idx).cloned().unwrap_or_default();
            fields.insert(mapped_col, val);
        }

        let mut rec = DwcRecord {
            row_number_1: row_num,
            fields,
            valid_lat_lon: None,
            has_errors: false,
            has_warnings: false,
        };

        let row_issue_start = issues.len();

        // 1. Validate / Mint occurrenceID
        validate_occurrence_id(&mut rec, &mut seen_ids, config, &mut issues);

        // 2. Validate / Normalize basisOfRecord
        validate_basis_of_record(&mut rec, config, &mut issues);

        // 3. Validate / Scrub ISO-8601 eventDate
        validate_event_date(&mut rec, config, &mut columns, &mut issues);

        // 4. Validate / Scrub Coordinates (DMS, transposed, Null Island, precision, datum)
        validate_coordinates(&mut rec, config, &mut columns, &mut issues);

        // 5. Validate elevation & individualCount
        validate_numeric_measurements(&rec, &mut issues);

        // 6. Validate & Resolve Scientific Name + Higher Taxonomy
        validate_taxonomy(&mut rec, backbone, config, &mut columns, &mut issues);

        for iss in &issues[row_issue_start..] {
            if iss.severity == LintSeverity::Error && !iss.auto_applied {
                rec.has_errors = true;
            } else if iss.severity == LintSeverity::Warning && !iss.auto_applied {
                rec.has_warnings = true;
            }
        }

        records.push(rec);
    }

    // 7. Dataset-wide per-species spatial outlier check (Zizka et al. 2019 MAD filter)
    if config.detect_species_outliers {
        detect_conspecific_spatial_outliers(&mut records, &mut issues);
    }

    Ok(DwcDataset {
        dataset_name: dataset_name.to_string(),
        delimiter,
        columns,
        header_mappings,
        records,
        issues,
    })
}

fn ensure_column(columns: &mut Vec<String>, col: &str) {
    if !columns.iter().any(|c| c == col) {
        columns.push(col.to_string());
    }
}

fn validate_occurrence_id(
    rec: &mut DwcRecord,
    seen_ids: &mut HashMap<String, usize>,
    config: &DwcCleanConfig,
    issues: &mut Vec<DwcLintIssue>,
) {
    let raw_id = rec.get("occurrenceID").trim().to_string();
    if raw_id.is_empty() {
        let minted = mint_deterministic_occurrence_id(&config.mint_id_prefix, rec);
        if config.auto_fix {
            rec.set("occurrenceID", minted.clone());
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("occurrenceID"),
            severity: LintSeverity::Warning,
            code: DwcIssueCode::MissingOccurrenceId,
            message: format!("Blank occurrenceID; minted deterministic ID '{minted}'."),
            original_value: String::new(),
            suggested_value: Some(minted),
            auto_applied: config.auto_fix,
        });
        return;
    }

    if let Some(first_row) = seen_ids.insert(raw_id.clone(), rec.row_number_1) {
        let dedup = format!("{}-dup-r{}", raw_id, rec.row_number_1);
        if config.auto_fix {
            rec.set("occurrenceID", dedup.clone());
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("occurrenceID"),
            severity: LintSeverity::Error,
            code: DwcIssueCode::DuplicateOccurrenceId,
            message: format!(
                "Duplicate occurrenceID '{}' already appeared on row #{}.",
                raw_id, first_row
            ),
            original_value: raw_id,
            suggested_value: Some(dedup),
            auto_applied: config.auto_fix,
        });
    }
}

fn mint_deterministic_occurrence_id(prefix: &str, rec: &DwcRecord) -> String {
    // FNV-1a 64-bit deterministic hash across row number + key content fields
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    let payload = format!(
        "{}|{}|{}|{}|{}",
        rec.row_number_1,
        rec.get("scientificName"),
        rec.get("eventDate"),
        rec.get("decimalLatitude"),
        rec.get("decimalLongitude")
    );
    for b in payload.bytes() {
        hash ^= b as u64;
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{}:{:04}:{:08x}", prefix, rec.row_number_1, (hash & 0xffff_ffff))
}

fn validate_basis_of_record(
    rec: &mut DwcRecord,
    config: &DwcCleanConfig,
    issues: &mut Vec<DwcLintIssue>,
) {
    const CANONICAL_BASIS: &[&str] = &[
        "HumanObservation",
        "PreservedSpecimen",
        "MachineObservation",
        "MaterialSample",
        "LivingSpecimen",
        "FossilSpecimen",
        "MaterialCitation",
        "Occurrence",
    ];
    let raw = rec.get("basisOfRecord").trim().to_string();
    if CANONICAL_BASIS.contains(&raw.as_str()) {
        return;
    }

    let key: String = raw
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect();

    let mapped = match key.as_str() {
        "humanobservation" | "observation" | "obs" | "visual" | "fieldobservation" | "sighting" => {
            Some("HumanObservation")
        }
        "preservedspecimen" | "specimen" | "voucher" | "herbarium" | "museumspecimen" => {
            Some("PreservedSpecimen")
        }
        "machineobservation" | "camera" | "cameratrap" | "acoustic" | "audiomoth" | "sensor" => {
            Some("MachineObservation")
        }
        "materialsample" | "dna" | "edna" | "tissue" | "blood" | "soil" | "sample" => {
            Some("MaterialSample")
        }
        "livingspecimen" | "living" | "culture" | "captive" => Some("LivingSpecimen"),
        "fossilspecimen" | "fossil" => Some("FossilSpecimen"),
        "materialcitation" | "literature" | "citation" => Some("MaterialCitation"),
        "occurrence" => Some("Occurrence"),
        _ => None,
    };

    if let Some(canon) = mapped {
        if config.auto_fix {
            rec.set("basisOfRecord", canon);
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("basisOfRecord"),
            severity: LintSeverity::Info,
            code: DwcIssueCode::NormalizedBasisOfRecord,
            message: format!(
                "Normalized basisOfRecord '{}' -> TDWG UpperCamelCase '{}'.",
                raw, canon
            ),
            original_value: raw,
            suggested_value: Some(canon.to_string()),
            auto_applied: config.auto_fix,
        });
    } else {
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("basisOfRecord"),
            severity: LintSeverity::Error,
            code: DwcIssueCode::InvalidBasisOfRecord,
            message: format!(
                "Invalid basisOfRecord '{}'; expected TDWG controlled vocabulary (e.g. HumanObservation, PreservedSpecimen).",
                raw
            ),
            original_value: raw,
            suggested_value: Some(String::from("HumanObservation")),
            auto_applied: false,
        });
    }
}

/// Parses and normalizes dates into ISO-8601 (`YYYY-MM-DD`, `YYYY-MM`, `YYYY`, or interval).
pub fn normalize_iso_event_date(raw: &str) -> Result<(String, i32, Option<u32>, Option<u32>), String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(String::from("Blank eventDate."));
    }
    // Handle ISO interval YYYY-MM-DD/YYYY-MM-DD by validating both endpoints
    if let Some((left, right)) = trimmed.split_once('/') {
        if left.len() >= 4 && right.len() >= 4 && left.contains('-') && right.contains('-') {
            let (l_iso, y, m, d) = normalize_iso_event_date(left)?;
            let (r_iso, _, _, _) = normalize_iso_event_date(right)?;
            return Ok((format!("{l_iso}/{r_iso}"), y, m, d));
        }
    }

    // Strip optional ISO time suffix `T14:22:00Z` while preserving it if valid
    let (date_part, time_suffix) = match trimmed.split_once('T') {
        Some((d, t)) if !t.is_empty() => (d, Some(t)),
        _ => (trimmed, None),
    };

    // Year only (`YYYY`)
    if date_part.len() == 4 && date_part.chars().all(|c| c.is_ascii_digit()) {
        let y: i32 = date_part.parse().map_err(|_| "Invalid year")?;
        return Ok((format!("{y:04}"), y, None, None));
    }

    // `YYYY-MM`
    let sep = if date_part.contains('-') {
        '-'
    } else if date_part.contains('/') {
        '/'
    } else {
        return Err(format!("Unrecognized date format '{trimmed}'."));
    };

    let parts: Vec<&str> = date_part.split(sep).collect();
    if parts.len() == 2 && parts[0].len() == 4 {
        let y: i32 = parts[0].parse().map_err(|_| "Invalid year")?;
        let m: u32 = parts[1].parse().map_err(|_| "Invalid month")?;
        if !(1..=12).contains(&m) {
            return Err(format!("Month {m} out of range 1..=12."));
        }
        return Ok((format!("{y:04}-{m:02}"), y, Some(m), None));
    }

    if parts.len() == 3 {
        let (y, m, d): (i32, u32, u32) = if parts[0].len() == 4 {
            // YYYY-MM-DD or YYYY/MM/DD
            (
                parts[0].parse().map_err(|_| "Invalid year")?,
                parts[1].parse().map_err(|_| "Invalid month")?,
                parts[2].parse().map_err(|_| "Invalid day")?,
            )
        } else if parts[2].len() == 4 {
            // Ambiguous slash/dash date: if first token > 12 -> DD/MM/YYYY, else US field sheet MM/DD/YYYY
            let a: u32 = parts[0].parse().map_err(|_| "Invalid month/day")?;
            let b: u32 = parts[1].parse().map_err(|_| "Invalid day/month")?;
            let y: i32 = parts[2].parse().map_err(|_| "Invalid year")?;
            if a > 12 {
                (y, b, a)
            } else {
                (y, a, b)
            }
        } else {
            return Err(format!("Two-digit years are ambiguous in '{trimmed}'; use 4-digit YYYY."));
        };

        validate_calendar_ymd(y, m, d)?;
        let iso = match time_suffix {
            Some(ts) => format!("{y:04}-{m:02}-{d:02}T{ts}"),
            None => format!("{y:04}-{m:02}-{d:02}"),
        };
        return Ok((iso, y, Some(m), Some(d)));
    }

    Err(format!("Could not parse eventDate '{trimmed}' as ISO-8601."))
}

fn validate_calendar_ymd(year: i32, month: u32, day: u32) -> Result<(), String> {
    if !(1500..=2100).contains(&year) {
        return Err(format!("Year {year} outside plausible biological survey range (1500..=2100)."));
    }
    if !(1..=12).contains(&month) {
        return Err(format!("Month {month} out of range 1..=12."));
    }
    let leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    let max_d = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if leap {
                29
            } else {
                28
            }
        }
        _ => 31,
    };
    if day < 1 || day > max_d {
        return Err(format!(
            "Invalid calendar day {year:04}-{month:02}-{day:02} (max {max_d} days in month)."
        ));
    }
    Ok(())
}

fn validate_event_date(
    rec: &mut DwcRecord,
    config: &DwcCleanConfig,
    columns: &mut Vec<String>,
    issues: &mut Vec<DwcLintIssue>,
) {
    let raw = rec.get("eventDate").trim().to_string();
    match normalize_iso_event_date(&raw) {
        Ok((iso, year, month, day)) => {
            if year > config.max_valid_year {
                issues.push(DwcLintIssue {
                    row_number_1: rec.row_number_1,
                    column_term: String::from("eventDate"),
                    severity: LintSeverity::Error,
                    code: DwcIssueCode::FutureEventDate,
                    message: format!(
                        "eventDate '{}' lies beyond maximum configured survey year ({}).",
                        iso, config.max_valid_year
                    ),
                    original_value: raw.clone(),
                    suggested_value: None,
                    auto_applied: false,
                });
            }
            if iso != raw {
                if config.auto_fix {
                    rec.set("eventDate", iso.clone());
                }
                issues.push(DwcLintIssue {
                    row_number_1: rec.row_number_1,
                    column_term: String::from("eventDate"),
                    severity: LintSeverity::Info,
                    code: DwcIssueCode::NormalizedEventDate,
                    message: format!("Reformatted eventDate '{}' -> ISO-8601 '{}'.", raw, iso),
                    original_value: raw,
                    suggested_value: Some(iso),
                    auto_applied: config.auto_fix,
                });
            }
            if config.auto_fix {
                ensure_column(columns, "year");
                rec.set("year", year.to_string());
                if let Some(m) = month {
                    ensure_column(columns, "month");
                    rec.set("month", m.to_string());
                }
                if let Some(d) = day {
                    ensure_column(columns, "day");
                    rec.set("day", d.to_string());
                }
            }
        }
        Err(reason) => {
            issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("eventDate"),
                severity: LintSeverity::Error,
                code: DwcIssueCode::InvalidEventDate,
                message: reason,
                original_value: raw,
                suggested_value: None,
                auto_applied: false,
            });
        }
    }
}

/// Parses either a decimal coordinate (`-123.2620`) or a DMS string (`44° 30' 14.4" N`, `123d33m03.6sW`).
/// Returns `(decimal_degrees, was_dms_converted, decimal_places)`.
pub fn parse_coordinate_or_dms(raw: &str) -> Result<(f64, bool, usize), String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err(String::from("Missing coordinate value."));
    }
    if let Ok(dd) = trimmed.parse::<f64>() {
        if !dd.is_finite() {
            return Err(String::from("Non-finite coordinate."));
        }
        let decimals = trimmed
            .split_once('.')
            .map(|(_, frac)| frac.chars().take_while(|c| c.is_ascii_digit()).count())
            .unwrap_or(0);
        return Ok((dd, false, decimals));
    }

    // Check for DMS tokens + compass hemisphere (N/S/E/W)
    let upper = trimmed.to_ascii_uppercase();
    let has_south_or_west = upper.contains('S') || upper.contains('W') || trimmed.starts_with('-');
    let mut nums: Vec<f64> = Vec::new();
    let mut cur = String::new();
    for ch in upper.chars() {
        if ch.is_ascii_digit() || ch == '.' {
            cur.push(ch);
        } else if !cur.is_empty() {
            let v: f64 = cur
                .parse()
                .map_err(|_| format!("Malformed DMS number in '{trimmed}'."))?;
            nums.push(v);
            cur.clear();
        }
    }
    if !cur.is_empty() {
        let v: f64 = cur
            .parse()
            .map_err(|_| format!("Malformed DMS number in '{trimmed}'."))?;
        nums.push(v);
    }

    if nums.is_empty() || nums.len() > 3 {
        return Err(format!("Cannot parse coordinate '{trimmed}' as decimal degrees or DMS."));
    }

    let deg = nums[0];
    let min = nums.get(1).copied().unwrap_or(0.0);
    let sec = nums.get(2).copied().unwrap_or(0.0);
    if min >= 60.0 || sec >= 60.0 {
        return Err(format!(
            "DMS minutes/seconds >= 60 in '{trimmed}' ({deg}° {min}' {sec}\")."
        ));
    }

    let mut dd = deg + (min / 60.0) + (sec / 3600.0);
    if has_south_or_west {
        dd = -dd.abs();
    }
    Ok((dd, true, 6))
}

fn validate_coordinates(
    rec: &mut DwcRecord,
    config: &DwcCleanConfig,
    columns: &mut Vec<String>,
    issues: &mut Vec<DwcLintIssue>,
) {
    let raw_lat = rec.get("decimalLatitude").trim().to_string();
    let raw_lon = rec.get("decimalLongitude").trim().to_string();

    let (mut lat, lat_dms, lat_dec) = match parse_coordinate_or_dms(&raw_lat) {
        Ok(v) => v,
        Err(e) => {
            issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("decimalLatitude"),
                severity: LintSeverity::Error,
                code: DwcIssueCode::CoordinateOutOfBounds,
                message: e,
                original_value: raw_lat,
                suggested_value: None,
                auto_applied: false,
            });
            return;
        }
    };

    let (mut lon, lon_dms, lon_dec) = match parse_coordinate_or_dms(&raw_lon) {
        Ok(v) => v,
        Err(e) => {
            issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("decimalLongitude"),
                severity: LintSeverity::Error,
                code: DwcIssueCode::CoordinateOutOfBounds,
                message: e,
                original_value: raw_lon,
                suggested_value: None,
                auto_applied: false,
            });
            return;
        }
    };

    if lat_dms {
        let formatted = format!("{lat:.6}");
        if config.auto_fix {
            rec.set("decimalLatitude", formatted.clone());
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("decimalLatitude"),
            severity: LintSeverity::Info,
            code: DwcIssueCode::ConvertedDmsCoordinate,
            message: format!("Converted DMS latitude '{}' -> WGS84 decimal degrees {}.", raw_lat, formatted),
            original_value: raw_lat.clone(),
            suggested_value: Some(formatted),
            auto_applied: config.auto_fix,
        });
    }

    if lon_dms {
        let formatted = format!("{lon:.6}");
        if config.auto_fix {
            rec.set("decimalLongitude", formatted.clone());
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("decimalLongitude"),
            severity: LintSeverity::Info,
            code: DwcIssueCode::ConvertedDmsCoordinate,
            message: format!(
                "Converted DMS longitude '{}' -> WGS84 decimal degrees {}.",
                raw_lon, formatted
            ),
            original_value: raw_lon.clone(),
            suggested_value: Some(formatted),
            auto_applied: config.auto_fix,
        });
    }

    // Detect transposed Latitude / Longitude (|lat| > 90 && |lat| <= 180 && |lon| <= 90)
    if lat.abs() > 90.0 && lat.abs() <= 180.0 && lon.abs() <= 90.0 {
        let swapped = format!("lat={lon:.6}, lon={lat:.6}");
        if config.auto_fix {
            std::mem::swap(&mut lat, &mut lon);
            rec.set("decimalLatitude", format!("{lat:.6}"));
            rec.set("decimalLongitude", format!("{lon:.6}"));
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("decimalLatitude"),
            severity: LintSeverity::Error,
            code: DwcIssueCode::TransposedLatLonSuspected,
            message: format!(
                "Latitude ({raw_lat}) exceeds +-90 deg while longitude ({raw_lon}) is within +-90 deg; transposed lat/lon columns."
            ),
            original_value: format!("{raw_lat}, {raw_lon}"),
            suggested_value: Some(swapped),
            auto_applied: config.auto_fix,
        });
    }

    if !(-90.0..=90.0).contains(&lat) || !(-180.0..=180.0).contains(&lon) {
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("decimalLatitude"),
            severity: LintSeverity::Error,
            code: DwcIssueCode::CoordinateOutOfBounds,
            message: format!(
                "Coordinates ({lat}, {lon}) fall outside WGS84 bounds (lat -90..=90, lon -180..=180)."
            ),
            original_value: format!("{raw_lat}, {raw_lon}"),
            suggested_value: None,
            auto_applied: false,
        });
        return;
    }

    // Null Island / exact zero checks
    if haversine_km(lat, lon, 0.0, 0.0) <= 5.0 || lat.abs() < 1e-7 || lon.abs() < 1e-7 {
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("decimalLatitude"),
            severity: LintSeverity::Error,
            code: DwcIssueCode::ZeroOrNullIslandCoordinate,
            message: format!(
                "Suspicious zero / Null Island coordinate ({lat:.5}, {lon:.5}) flagged by CoordinateCleaner zero-check."
            ),
            original_value: format!("{lat}, {lon}"),
            suggested_value: None,
            auto_applied: false,
        });
        return;
    }

    // Equal |lat| == |lon| artifact check
    if (lat - lon).abs() < 1e-6 || (lat.abs() - lon.abs()).abs() < 1e-6 {
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("decimalLatitude"),
            severity: LintSeverity::Warning,
            code: DwcIssueCode::EqualLatLonCoordinate,
            message: format!(
                "Identical |latitude| and |longitude| ({lat:.5}, {lon:.5}); verify copy-paste error."
            ),
            original_value: format!("{lat}, {lon}"),
            suggested_value: None,
            auto_applied: false,
        });
    }

    // Check against known centroids & biodiversity institutions
    for &(label, c_lat, c_lon, radius_km) in KNOWN_CENTROIDS_AND_INSTITUTIONS.iter().skip(1) {
        let d_km = haversine_km(lat, lon, c_lat, c_lon);
        if d_km <= radius_km {
            issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("decimalLatitude"),
                severity: LintSeverity::Warning,
                code: DwcIssueCode::InstitutionOrCentroidCoordinate,
                message: format!(
                    "Coordinate lies {:.2} km from '{}' — likely automated centroid/institution geocoding artifact.",
                    d_km, label
                ),
                original_value: format!("{lat}, {lon}"),
                suggested_value: None,
                auto_applied: false,
            });
        }
    }

    // Precision & uncertainty check
    let min_dec = lat_dec.min(lon_dec);
    if min_dec < config.min_coord_decimals
        && rec.get("coordinateUncertaintyInMeters").trim().is_empty()
    {
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("coordinateUncertaintyInMeters"),
            severity: LintSeverity::Warning,
            code: DwcIssueCode::LowCoordinatePrecision,
            message: format!(
                "Low coordinate precision ({min_dec} decimal places, ~{:.1} km grid) with blank coordinateUncertaintyInMeters.",
                111.32 / 10_f64.powi(min_dec as i32)
            ),
            original_value: format!("{lat}, {lon}"),
            suggested_value: Some(format!(
                "{:.0}",
                (111_320.0 / 10_f64.powi(min_dec as i32)).round()
            )),
            auto_applied: false,
        });
    }

    // Populate WGS84 datum if empty
    if rec.get("geodeticDatum").trim().is_empty() {
        if config.auto_fix {
            ensure_column(columns, "geodeticDatum");
            rec.set("geodeticDatum", "WGS84");
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("geodeticDatum"),
            severity: LintSeverity::Info,
            code: DwcIssueCode::MissingGeodeticDatum,
            message: String::from("Missing geodeticDatum on decimal coordinates; defaulted to 'WGS84'."),
            original_value: String::new(),
            suggested_value: Some(String::from("WGS84")),
            auto_applied: config.auto_fix,
        });
    }

    rec.valid_lat_lon = Some((lat, lon));
}

fn validate_numeric_measurements(rec: &DwcRecord, issues: &mut Vec<DwcLintIssue>) {
    let elev_raw = rec.get("minimumElevationInMeters").trim();
    if !elev_raw.is_empty() {
        match elev_raw.parse::<f64>() {
            Ok(el) if (-430.0..=8850.0).contains(&el) => {}
            _ => issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("minimumElevationInMeters"),
                severity: LintSeverity::Warning,
                code: DwcIssueCode::ElevationOutOfBounds,
                message: format!(
                    "minimumElevationInMeters '{}' is non-numeric or outside terrestrial bounds (-430..=8850 m).",
                    elev_raw
                ),
                original_value: elev_raw.to_string(),
                suggested_value: None,
                auto_applied: false,
            }),
        }
    }

    let cnt_raw = rec.get("individualCount").trim();
    if !cnt_raw.is_empty() {
        match cnt_raw.parse::<i64>() {
            Ok(n) if n >= 0 => {}
            _ => issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("individualCount"),
                severity: LintSeverity::Warning,
                code: DwcIssueCode::InvalidIndividualCount,
                message: format!("individualCount '{}' must be a non-negative integer.", cnt_raw),
                original_value: cnt_raw.to_string(),
                suggested_value: None,
                auto_applied: false,
            }),
        }
    }
}

fn validate_taxonomy(
    rec: &mut DwcRecord,
    backbone: &TaxonBackbone,
    config: &DwcCleanConfig,
    columns: &mut Vec<String>,
    issues: &mut Vec<DwcLintIssue>,
) {
    let raw_sci = rec.get("scientificName").trim().to_string();
    if raw_sci.is_empty() {
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("scientificName"),
            severity: LintSeverity::Error,
            code: DwcIssueCode::UnmatchedScientificName,
            message: String::from("Blank scientificName."),
            original_value: String::new(),
            suggested_value: None,
            auto_applied: false,
        });
        return;
    }

    let res = resolve_scientific_name(&raw_sci, backbone);

    // Split embedded authorship into `scientificNameAuthorship` when missing
    if !res.parsed.authorship.is_empty() && rec.get("scientificNameAuthorship").trim().is_empty() {
        if config.auto_fix {
            ensure_column(columns, "scientificNameAuthorship");
            rec.set("scientificNameAuthorship", res.parsed.authorship.clone());
        }
        issues.push(DwcLintIssue {
            row_number_1: rec.row_number_1,
            column_term: String::from("scientificNameAuthorship"),
            severity: LintSeverity::Info,
            code: DwcIssueCode::TaxonAuthorshipExtracted,
            message: format!(
                "Extracted authorship '{}' from scientificName '{}'.",
                res.parsed.authorship, raw_sci
            ),
            original_value: raw_sci.clone(),
            suggested_value: Some(res.parsed.authorship.clone()),
            auto_applied: config.auto_fix,
        });
    }

    match res.match_kind {
        TaxonMatchKind::ExactAccepted => {
            if let Some(entry) = res.matched_entry {
                if config.auto_fix {
                    populate_higher_taxonomy(rec, &entry, columns);
                }
            }
        }
        TaxonMatchKind::Synonym => {
            if let Some(entry) = res.matched_entry {
                if config.auto_fix {
                    ensure_column(columns, "acceptedNameUsage");
                    ensure_column(columns, "taxonomicStatus");
                    rec.set("acceptedNameUsage", entry.accepted_name.clone());
                    rec.set("taxonomicStatus", if config.replace_synonyms_in_place {
                        "ACCEPTED"
                    } else {
                        "SYNONYM"
                    });
                    if config.replace_synonyms_in_place {
                        rec.set("scientificName", entry.accepted_name.clone());
                    }
                    populate_higher_taxonomy(rec, &entry, columns);
                }
                issues.push(DwcLintIssue {
                    row_number_1: rec.row_number_1,
                    column_term: String::from("scientificName"),
                    severity: LintSeverity::Warning,
                    code: DwcIssueCode::TaxonSynonymResolved,
                    message: format!(
                        "Taxonomic synonym '{}' resolves to accepted name '{}' (GBIF key {}).",
                        res.parsed.canonical_name, entry.accepted_name, entry.gbif_usage_key
                    ),
                    original_value: raw_sci,
                    suggested_value: Some(entry.accepted_name),
                    auto_applied: config.auto_fix,
                });
            }
        }
        TaxonMatchKind::FuzzyOrthographic { edit_distance } => {
            if let Some(entry) = res.matched_entry {
                if config.auto_fix {
                    rec.set("scientificName", entry.accepted_name.clone());
                    populate_higher_taxonomy(rec, &entry, columns);
                }
                issues.push(DwcLintIssue {
                    row_number_1: rec.row_number_1,
                    column_term: String::from("scientificName"),
                    severity: LintSeverity::Warning,
                    code: DwcIssueCode::TaxonFuzzyTypoSuggested,
                    message: format!(
                        "Probable spelling typo '{}' (edit distance {}) -> '{}'.",
                        res.parsed.canonical_name, edit_distance, entry.accepted_name
                    ),
                    original_value: raw_sci,
                    suggested_value: Some(entry.accepted_name),
                    auto_applied: config.auto_fix,
                });
            }
        }
        TaxonMatchKind::Unmatched => {
            if config.auto_fix && rec.get("taxonRank").trim().is_empty() {
                ensure_column(columns, "taxonRank");
                rec.set("taxonRank", res.parsed.inferred_rank);
            }
            issues.push(DwcLintIssue {
                row_number_1: rec.row_number_1,
                column_term: String::from("scientificName"),
                severity: LintSeverity::Info,
                code: DwcIssueCode::UnmatchedScientificName,
                message: format!(
                    "Taxon '{}' not in local offline backbone; verify via GBIF /v1/species/match or load regional checklist TSV.",
                    res.parsed.canonical_name
                ),
                original_value: raw_sci,
                suggested_value: Some(res.gbif_lookup_url),
                auto_applied: false,
            });
        }
    }
}

fn populate_higher_taxonomy(
    rec: &mut DwcRecord,
    entry: &crate::dwc_taxonomy::BackboneTaxon,
    columns: &mut Vec<String>,
) {
    for (term, val) in [
        ("taxonID", format!("GBIF:{}", entry.gbif_usage_key)),
        ("kingdom", entry.kingdom.clone()),
        ("phylum", entry.phylum.clone()),
        ("class", entry.class_name.clone()),
        ("order", entry.order.clone()),
        ("family", entry.family.clone()),
        ("genus", entry.genus.clone()),
        ("taxonRank", entry.rank.clone()),
    ] {
        if rec.get(term).trim().is_empty() && !val.is_empty() {
            ensure_column(columns, term);
            rec.set(term, val);
        }
    }
}

/// Great-circle Haversine distance in kilometers on WGS84 mean spherical radius (6371.0088 km).
pub fn haversine_km(lat1_deg: f64, lon1_deg: f64, lat2_deg: f64, lon2_deg: f64) -> f64 {
    const R_KM: f64 = 6371.0088;
    let dlat = (lat2_deg - lat1_deg) * (PI / 180.0);
    let dlon = (lon2_deg - lon1_deg) * (PI / 180.0);
    let rlat1 = lat1_deg * (PI / 180.0);
    let rlat2 = lat2_deg * (PI / 180.0);
    let a = (dlat * 0.5).sin().powi(2) + rlat1.cos() * rlat2.cos() * (dlon * 0.5).sin().powi(2);
    2.0 * R_KM * a.sqrt().min(1.0).asin()
}

fn detect_conspecific_spatial_outliers(
    records: &mut [DwcRecord],
    issues: &mut Vec<DwcLintIssue>,
) {
    let mut by_taxon: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, rec) in records.iter().enumerate() {
        if rec.valid_lat_lon.is_some() {
            let key = if !rec.get("acceptedNameUsage").is_empty() {
                rec.get("acceptedNameUsage").to_string()
            } else {
                rec.get("scientificName").to_string()
            };
            if !key.is_empty() {
                by_taxon.entry(key).or_default().push(i);
            }
        }
    }

    for (taxon, indices) in by_taxon {
        if indices.len() < 4 {
            continue;
        }
        let mut lats: Vec<f64> = indices
            .iter()
            .filter_map(|&i| records[i].valid_lat_lon.map(|p| p.0))
            .collect();
        let mut lons: Vec<f64> = indices
            .iter()
            .filter_map(|&i| records[i].valid_lat_lon.map(|p| p.1))
            .collect();
        lats.sort_by(|a, b| a.partial_cmp(b).unwrap());
        lons.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med_lat = lats[lats.len() / 2];
        let med_lon = lons[lons.len() / 2];

        let mut dists: Vec<f64> = indices
            .iter()
            .map(|&i| {
                let (la, lo) = records[i].valid_lat_lon.unwrap();
                haversine_km(la, lo, med_lat, med_lon)
            })
            .collect();
        let mut sorted_d = dists.clone();
        sorted_d.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med_d = sorted_d[sorted_d.len() / 2];

        let mut abs_dev: Vec<f64> = sorted_d.iter().map(|d| (d - med_d).abs()).collect();
        abs_dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mad_km = (abs_dev[abs_dev.len() / 2] * 1.4826).max(50.0);
        let threshold_km = (med_d + 5.0 * mad_km).max(350.0);

        for (pos, &rec_idx) in indices.iter().enumerate() {
            let d_km = dists[pos];
            if d_km > threshold_km {
                records[rec_idx].has_warnings = true;
                let (la, lo) = records[rec_idx].valid_lat_lon.unwrap();
                issues.push(DwcLintIssue {
                    row_number_1: records[rec_idx].row_number_1,
                    column_term: String::from("decimalLatitude"),
                    severity: LintSeverity::Warning,
                    code: DwcIssueCode::SpeciesSpatialOutlier,
                    message: format!(
                        "Geographic outlier for '{}': {:.0} km from conspecific median centroid ({:.3}, {:.3}).",
                        taxon, d_km, med_lat, med_lon
                    ),
                    original_value: format!("{la:.5}, {lo:.5}"),
                    suggested_value: None,
                    auto_applied: false,
                });
            }
        }
        let _ = &mut dists;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cleans_demo_pnw_survey_dms_transposed_synonyms_and_null_island() {
        let bb = TaxonBackbone::default_curated();
        let cfg = DwcCleanConfig::default();
        let ds = lint_and_clean_dwc(
            "PNW_Demo",
            DwcDataset::demo_pnw_field_survey_csv(),
            &bb,
            &cfg,
        )
        .unwrap();

        assert_eq!(ds.records.len(), 8);
        // Row 1: DMS converted + slash date normalized + Pseudotsuga taxifolia synonym mapped
        let r1 = &ds.records[0];
        assert_eq!(r1.get("eventDate"), "2026-07-14");
        assert_eq!(r1.get("basisOfRecord"), "HumanObservation");
        assert_eq!(r1.get("acceptedNameUsage"), "Pseudotsuga menziesii");
        let (lat1, lon1) = r1.valid_lat_lon.unwrap();
        assert!((lat1 - 44.504).abs() < 0.01);
        assert!((lon1 - (-123.551)).abs() < 0.01);

        // Row 3: Transposed lat/lon swapped into Oregon (~44.5582, -123.2794) + minted occurrenceID
        let r3 = &ds.records[2];
        assert!(!r3.get("occurrenceID").is_empty());
        let (lat3, lon3) = r3.valid_lat_lon.unwrap();
        assert!((lat3 - 44.5582).abs() < 1e-3);
        assert!((lon3 - (-123.2794)).abs() < 1e-3);

        // Row 7: Null Island flagged as error & fuzzy typo fixed -> Pseudotsuga menziesii
        let r7 = &ds.records[6];
        assert!(r7.has_errors);
        assert!(r7.valid_lat_lon.is_none());
        assert_eq!(r7.get("scientificName"), "Pseudotsuga menziesii");
    }
}
