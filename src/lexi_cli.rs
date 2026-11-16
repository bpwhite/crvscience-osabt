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

//! Command-line subcommands (`clap` v4) for Tool #6: `osabt lexi`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::lexi_engine::{
    build_cognate_tree, decode_biological_phrase, demo_etymology_terms, export_anki_flashcards_tsv,
    export_flashcards_markdown, search_morpheme_lexicon, DecodedPhrase,
};

#[derive(Debug, Args)]
pub struct LexiArgs {
    #[command(subcommand)]
    pub command: Option<LexiSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum LexiSubcommand {
    /// Launch the interactive 4-tab Ratatui Greek/Latin etymology & flashcard quiz studio.
    Tui {
        /// Optional initial binomial or anatomical term to inspect.
        #[arg(short, long)]
        term: Option<String>,
    },

    /// Decode one or more taxonomic binomials or anatomical terms into Greek/Latin morphemes.
    Decode {
        /// Scientific names or anatomical terms (uses built-in showcase list if omitted).
        terms: Vec<String>,
    },

    /// Render the shared-morpheme cognate tree for a biological or anatomical term.
    Cognates {
        /// Term or binomial to expand into related cognates (e.g. `"Actinopterygii"`).
        #[arg(default_value = "Drosophila melanogaster")]
        term: String,
    },

    /// Search the offline Classical Greek, Latin, and ICZN/ICNafp morpheme lexicon.
    SearchRoot {
        /// Substring matching root form, English gloss, or cognate example (e.g. `"wing"` or `"pter"`).
        query: String,
    },

    /// Generate an Anki-importable `.tsv` deck (`#separator:tab`) and/or Markdown study sheet.
    Flashcards {
        /// Optional newline-delimited text file of scientific/anatomical terms (uses showcase set if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output Anki `.tsv` file path.
        #[arg(short, long, default_value = "osabt_etymology_anki.tsv")]
        anki_out: PathBuf,
        /// Optional Markdown study sheet output path.
        #[arg(long)]
        md_out: Option<PathBuf>,
        /// Anki deck tag column string.
        #[arg(long, default_value = "OSABT_Biology_Etymology")]
        tag: String,
    },
}

pub fn run_lexi_cli(args: LexiArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "tui")]
            {
                return crate::lexi_tui::run_lexi_tui(None);
            }
            #[cfg(not(feature = "tui"))]
            {
                return Err(String::from(
                    "Build with `--features tui` or run `osabt lexi decode`.",
                ));
            }
        }
        Some(LexiSubcommand::Tui { term }) => {
            #[cfg(feature = "tui")]
            {
                crate::lexi_tui::run_lexi_tui(term)
            }
            #[cfg(not(feature = "tui"))]
            {
                let _ = term;
                Err(String::from("TUI feature is disabled in this build."))
            }
        }
        Some(LexiSubcommand::Decode { terms }) => {
            let list: Vec<String> = if terms.is_empty() {
                demo_etymology_terms()
                    .into_iter()
                    .map(String::from)
                    .collect()
            } else {
                terms
            };

            for raw in &list {
                let dec = decode_biological_phrase(raw);
                println!("=== {} ===", dec.raw_input);
                println!("Literal Gloss : {}", dec.combined_gloss);
                for w in &dec.words {
                    println!("  Word '{}' -> {}", w.original_word, w.hyphenated);
                    for seg in &w.segments {
                        println!(
                            "    • {:<10} [{:<17} | {:<24}] -> {}",
                            seg.surface_slice, seg.origin_label, seg.classical_lemma, seg.gloss
                        );
                    }
                }
                println!();
            }
            Ok(())
        }
        Some(LexiSubcommand::Cognates { term }) => {
            let dec = decode_biological_phrase(&term);
            let nodes = build_cognate_tree(&dec);
            println!("=== Cognate Tree for '{}' ===", dec.raw_input);
            for node in &nodes {
                println!(
                    "┌─ {} ({} — {}) = \"{}\"",
                    node.morpheme_id, node.origin, node.classical_lemma, node.gloss
                );
                for (i, c) in node.cognates.iter().enumerate() {
                    let b = if i + 1 == node.cognates.len() {
                        "└──"
                    } else {
                        "├──"
                    };
                    println!("   {b} {c}");
                }
            }
            Ok(())
        }
        Some(LexiSubcommand::SearchRoot { query }) => {
            let hits = search_morpheme_lexicon(&query);
            println!(
                "Found {} morpheme entries matching '{}':",
                hits.len(),
                query
            );
            for e in hits {
                println!(
                    "  • {:<10} ({:<9} | {:<17} | {}) -> {}  [e.g. {}]",
                    e.canonical_id,
                    e.role.label(),
                    e.origin.label(),
                    e.classical_lemma,
                    e.gloss,
                    e.cognates.join(", ")
                );
            }
            Ok(())
        }
        Some(LexiSubcommand::Flashcards {
            input,
            anki_out,
            md_out,
            tag,
        }) => {
            let phrases: Vec<DecodedPhrase> = match input {
                Some(p) => {
                    let txt = fs::read_to_string(&p)
                        .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                    txt.lines()
                        .map(str::trim)
                        .filter(|l| !l.is_empty() && !l.starts_with('#'))
                        .map(decode_biological_phrase)
                        .collect()
                }
                None => demo_etymology_terms()
                    .into_iter()
                    .map(decode_biological_phrase)
                    .collect(),
            };

            fs::write(&anki_out, export_anki_flashcards_tsv(&phrases, &tag))
                .map_err(|e| format!("Failed writing Anki TSV '{}': {e}", anki_out.display()))?;
            println!(
                "Wrote {} Anki flashcards -> {}",
                phrases.len(),
                anki_out.display()
            );

            if let Some(mp) = md_out {
                fs::write(&mp, export_flashcards_markdown(&phrases))
                    .map_err(|e| format!("Failed writing Markdown '{}': {e}", mp.display()))?;
                println!("Wrote Markdown study sheet -> {}", mp.display());
            }
            Ok(())
        }
    }
}