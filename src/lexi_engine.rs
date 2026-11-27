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

//! Tool #6 (`osabt lexi`) — Dynamic-programming Greek/Latin morpheme lattice decoder,
//! patronymic/toponymic heuristic recognizer, cognate tree builder, and Anki TSV / Markdown
//! study flashcard exporter.

use std::fmt::Write as _;

use crate::lexi_roots::{
    MorphemeEntry, MorphemeOrigin, MorphemeRole, LEXI_MORPHEME_CATALOG,
};

/// Single segmented piece of a decoded biological or anatomical word.
#[derive(Debug, Clone)]
pub struct DecodedSegment {
    pub surface_slice: String,
    pub canonical_id: String,
    pub role_label: String,
    pub origin_label: String,
    pub classical_lemma: String,
    pub gloss: String,
    pub cognates: Vec<String>,
    pub is_connecting_vowel: bool,
}

/// Complete etymological breakdown of one word token.
#[derive(Debug, Clone)]
pub struct DecodedWord {
    pub original_word: String,
    pub hyphenated: String,
    pub literal_gloss: String,
    pub segments: Vec<DecodedSegment>,
}

/// Complete breakdown of a single or multi-word binomial/anatomical phrase.
#[derive(Debug, Clone)]
pub struct DecodedPhrase {
    pub raw_input: String,
    pub words: Vec<DecodedWord>,
    pub combined_gloss: String,
}

/// Cognate cluster rooted on a shared Classical Greek/Latin morpheme.
#[derive(Debug, Clone)]
pub struct CognateNode {
    pub morpheme_id: String,
    pub classical_lemma: String,
    pub gloss: String,
    pub origin: String,
    pub cognates: Vec<String>,
}

/// Decodes a single word or multi-word scientific/anatomical phrase (e.g. `"Drosophila melanogaster"`
/// or `"Sternocleidomastoideus"`) into constituent Greek/Latin morphemes and literal meaning.
pub fn decode_biological_phrase(raw_phrase: &str) -> DecodedPhrase {
    let words: Vec<DecodedWord> = raw_phrase
        .split_whitespace()
        .filter(|t| !t.is_empty())
        .map(|tok| {
            let clean: String = tok
                .chars()
                .filter(|c| c.is_ascii_alphabetic() || *c == '-')
                .collect();
            decode_single_word(&clean)
        })
        .collect();

    let combined = words
        .iter()
        .map(|w| format!("{} (\"{}\")", w.original_word, w.literal_gloss))
        .collect::<Vec<_>>()
        .join(" + ");

    DecodedPhrase {
        raw_input: raw_phrase.trim().to_string(),
        words,
        combined_gloss: combined,
    }
}

/// Segments a single biological/anatomical word via dynamic programming over `LEXI_MORPHEME_CATALOG`,
/// handling `-o-`/`-i-` connecting vowels plus `-ii`/`-ae` patronymics and `-ensis` toponyms.
pub fn decode_single_word(word: &str) -> DecodedWord {
    let low = word.trim().to_ascii_lowercase().replace('-', "");
    let n = low.len();
    if n == 0 {
        return DecodedWord {
            original_word: word.to_string(),
            hyphenated: String::new(),
            literal_gloss: String::new(),
            segments: Vec::new(),
        };
    }

    // Special fast-path check for Neo-Latin eponymous genitives (`menziesii`, `cascadae`) or
    // geographic toponyms (`mississippiensis`, `oregonensis`) when no Classical stem covers the prefix
    if let Some(patronym) = try_patronym_or_toponym(word, &low) {
        return patronym;
    }

    #[derive(Clone)]
    struct Edge {
        prev: usize,
        cost: i64,
        seg: DecodedSegment,
    }

    let mut dp: Vec<Option<Edge>> = vec![None; n + 1];
    dp[0] = Some(Edge {
        prev: 0,
        cost: 0,
        seg: DecodedSegment {
            surface_slice: String::new(),
            canonical_id: String::new(),
            role_label: String::new(),
            origin_label: String::new(),
            classical_lemma: String::new(),
            gloss: String::new(),
            cognates: Vec::new(),
            is_connecting_vowel: false,
        },
    });

    for i in 0..n {
        let Some(cur_edge) = dp[i].clone() else {
            continue;
        };
        let base_cost = cur_edge.cost;

        // 1. Try matching every catalog morpheme form starting at byte `i`
        for entry in LEXI_MORPHEME_CATALOG {
            if entry.role == MorphemeRole::Prefix && i > 0 {
                continue;
            }
            for &form in entry.forms {
                let flen = form.len();
                if flen == 0 || i + flen > n {
                    continue;
                }
                if &low[i..(i + flen)] != form {
                    continue;
                }
                if entry.role == MorphemeRole::Suffix && i + flen != n {
                    continue;
                }
                if entry.role == MorphemeRole::Suffix && i < 2 {
                    continue;
                }

                // Reward longer verified Classical stems strongly (`10 - 14 * flen`) so single long
                // roots beat choppy 1-2 char accidental matches.
                let role_adj: i64 = match entry.role {
                    MorphemeRole::Stem => -8,
                    MorphemeRole::Prefix => -4,
                    MorphemeRole::Suffix => {
                        if flen >= 4 {
                            -10
                        } else {
                            6 // Keep 1-2 char inflectional endings (`-a`, `-us`) lower priority than real stems
                        }
                    }
                };
                let step_cost = 18 - (15 * (flen as i64)) + role_adj;
                let cand_cost = base_cost + step_cost;
                let j = i + flen;

                if dp[j].as_ref().map(|e| cand_cost < e.cost).unwrap_or(true) {
                    dp[j] = Some(Edge {
                        prev: i,
                        cost: cand_cost,
                        seg: DecodedSegment {
                            surface_slice: low[i..j].to_string(),
                            canonical_id: entry.canonical_id.to_string(),
                            role_label: entry.role.label().to_string(),
                            origin_label: entry.origin.label().to_string(),
                            classical_lemma: entry.classical_lemma.to_string(),
                            gloss: entry.gloss.to_string(),
                            cognates: entry.cognates.iter().map(|s| (*s).to_string()).collect(),
                            is_connecting_vowel: false,
                        },
                    });
                }
            }
        }

        // 2. Single connecting vowel (`-o-` or `-i-`) internally between morphemes
        let ch = low.as_bytes()[i] as char;
        if i > 1 && i + 1 < n && (ch == 'o' || ch == 'i') {
            let cv_cost = base_cost + 4;
            if dp[i + 1].as_ref().map(|e| cv_cost < e.cost).unwrap_or(true) {
                dp[i + 1] = Some(Edge {
                    prev: i,
                    cost: cv_cost,
                    seg: DecodedSegment {
                        surface_slice: ch.to_string(),
                        canonical_id: format!("-{ch}-"),
                        role_label: String::from("linker"),
                        origin_label: if ch == 'o' {
                            String::from("Greek linker")
                        } else {
                            String::from("Latin linker")
                        },
                        classical_lemma: format!("-{ch}-"),
                        gloss: String::from("connecting vowel"),
                        cognates: Vec::new(),
                        is_connecting_vowel: true,
                    },
                });
            }
        }

        // 3. Fallback unknown 1-char advance (high penalty = +35 per char)
        let unk_cost = base_cost + 35;
        if dp[i + 1]
            .as_ref()
            .map(|e| unk_cost < e.cost)
            .unwrap_or(true)
        {
            dp[i + 1] = Some(Edge {
                prev: i,
                cost: unk_cost,
                seg: DecodedSegment {
                    surface_slice: low[i..(i + 1)].to_string(),
                    canonical_id: String::from("?"),
                    role_label: String::from("stem"),
                    origin_label: String::from("Neo-Latin/Vernacular"),
                    classical_lemma: String::from("(unindexed stem)"),
                    gloss: low[i..(i + 1)].to_string(),
                    cognates: Vec::new(),
                    is_connecting_vowel: false,
                },
            });
        }
    }

    // Backtrack from `n` to `0` and coalesce consecutive unknown `?` characters
    let mut raw_segs = Vec::new();
    let mut cur = n;
    while cur > 0 {
        let Some(edge) = dp[cur].clone() else {
            break;
        };
        raw_segs.push(edge.seg);
        cur = edge.prev;
    }
    raw_segs.reverse();

    let mut segments: Vec<DecodedSegment> = Vec::new();
    for seg in raw_segs {
        if seg.canonical_id == "?" {
            if let Some(last) = segments.last_mut() {
                if last.canonical_id == "?" {
                    last.surface_slice.push_str(&seg.surface_slice);
                    last.gloss = format!("proper/vernacular stem '{}'", last.surface_slice);
                    continue;
                }
            }
            let mut s = seg;
            s.gloss = format!("proper/vernacular stem '{}'", s.surface_slice);
            segments.push(s);
        } else {
            segments.push(seg);
        }
    }

    let hyphenated = segments
        .iter()
        .map(|s| s.surface_slice.as_str())
        .collect::<Vec<_>>()
        .join("-");

    let meaningful_glosses: Vec<&str> = segments
        .iter()
        .filter(|s| !s.is_connecting_vowel)
        .map(|s| s.gloss.as_str())
        .collect();
    let literal_gloss = meaningful_glosses.join(" + ");

    DecodedWord {
        original_word: word.to_string(),
        hyphenated,
        literal_gloss,
        segments,
    }
}

fn try_patronym_or_toponym(original: &str, low: &str) -> Option<DecodedWord> {
    // Toponym ending `-ensis` or `-ense` where stem isn't a Classical root (e.g. `mississippiensis`, `corvallisensis`)
    for suf in ["ensis", "ense"] {
        if let Some(stem) = low.strip_suffix(suf) {
            if stem.len() >= 4 && !LEXI_MORPHEME_CATALOG.iter().any(|e| e.forms.contains(&stem)) {
                let place = titlecase(stem);
                return Some(DecodedWord {
                    original_word: original.to_string(),
                    hyphenated: format!("{stem}-{suf}"),
                    literal_gloss: format!("native to / originating from {place}"),
                    segments: vec![
                        DecodedSegment {
                            surface_slice: stem.to_string(),
                            canonical_id: stem.to_string(),
                            role_label: String::from("toponym"),
                            origin_label: String::from("Geographic"),
                            classical_lemma: place.clone(),
                            gloss: format!("locality '{place}'"),
                            cognates: Vec::new(),
                            is_connecting_vowel: false,
                        },
                        DecodedSegment {
                            surface_slice: suf.to_string(),
                            canonical_id: String::from("-ensis"),
                            role_label: String::from("suffix"),
                            origin_label: String::from("Latin"),
                            classical_lemma: String::from("-ensis / -ense"),
                            gloss: String::from("originating in, native to"),
                            cognates: vec![
                                String::from("Alligator mississippiensis"),
                                String::from("Homo neanderthalensis"),
                            ],
                            is_connecting_vowel: false,
                        },
                    ],
                });
            }
        }
    }

    // Eponymous male genitive `-ii` (e.g. `menziesii`, `townsendii`, `clarkii`) when stem is not in Classical catalog
    if let Some(stem) = low.strip_suffix("ii") {
        let known_classical = LEXI_MORPHEME_CATALOG
            .iter()
            .any(|e| e.forms.iter().any(|f| low.starts_with(f) && f.len() >= 4));
        if stem.len() >= 4 && !known_classical {
            let person = titlecase(stem);
            return Some(DecodedWord {
                original_word: original.to_string(),
                hyphenated: format!("{stem}-ii"),
                literal_gloss: format!("named in honor of {person}"),
                segments: vec![
                    DecodedSegment {
                        surface_slice: stem.to_string(),
                        canonical_id: stem.to_string(),
                        role_label: String::from("eponym"),
                        origin_label: String::from("Neo-Latin Eponym"),
                        classical_lemma: person.clone(),
                        gloss: format!("naturalist/person '{person}'"),
                        cognates: Vec::new(),
                        is_connecting_vowel: false,
                    },
                    DecodedSegment {
                        surface_slice: String::from("ii"),
                        canonical_id: String::from("-ii"),
                        role_label: String::from("suffix"),
                        origin_label: MorphemeOrigin::NeoLatinTaxonomic.label().to_string(),
                        classical_lemma: String::from("-ii (ICZN/ICNafp genitive)"),
                        gloss: String::from("named in honor of"),
                        cognates: vec![String::from("Pseudotsuga menziesii")],
                        is_connecting_vowel: false,
                    },
                ],
            });
        }
    }

    None
}

fn titlecase(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        Some(f) => f.to_ascii_uppercase().to_string() + c.as_str(),
        None => String::new(),
    }
}

/// Searches `LEXI_MORPHEME_CATALOG` by root form, English gloss, or cognate substring.
pub fn search_morpheme_lexicon(query: &str) -> Vec<&'static MorphemeEntry> {
    let q = query.trim().trim_matches('-').to_ascii_lowercase();
    if q.is_empty() {
        return LEXI_MORPHEME_CATALOG.iter().collect();
    }
    LEXI_MORPHEME_CATALOG
        .iter()
        .filter(|e| {
            e.canonical_id.to_ascii_lowercase().contains(&q)
                || e.forms.iter().any(|f| f.contains(&q))
                || e.gloss.to_ascii_lowercase().contains(&q)
                || e.cognates
                    .iter()
                    .any(|c| c.to_ascii_lowercase().contains(&q))
        })
        .collect()
}

/// Builds the cognate network for every Classical root discovered inside `decoded`.
pub fn build_cognate_tree(decoded: &DecodedPhrase) -> Vec<CognateNode> {
    let mut nodes = Vec::new();
    let mut seen = Vec::new();

    for w in &decoded.words {
        for seg in &w.segments {
            if seg.is_connecting_vowel || seg.cognates.is_empty() {
                continue;
            }
            if seen.contains(&seg.canonical_id) {
                continue;
            }
            seen.push(seg.canonical_id.clone());
            nodes.push(CognateNode {
                morpheme_id: seg.canonical_id.clone(),
                classical_lemma: seg.classical_lemma.clone(),
                gloss: seg.gloss.clone(),
                origin: seg.origin_label.clone(),
                cognates: seg.cognates.clone(),
            });
        }
    }
    nodes
}

/// Exports decoded terms as an Anki-ready tab-separated deck (`#separator:tab`, `#html:true`)
/// conforming to `https://docs.ankiweb.net/manual/importing/text-files`.
pub fn export_anki_flashcards_tsv(phrases: &[DecodedPhrase], tag: &str) -> String {
    let mut out = String::from("#separator:tab\n#html:true\n#tags column:3\n");
    let clean_tag = tag.replace(char::is_whitespace, "_");

    for p in phrases {
        let front = format!("<b><i>{}</i></b> — Etymology &amp; Literal Meaning?", p.raw_input);
        let mut back = String::new();
        for w in &p.words {
            let _ = write!(
                back,
                "<div><b>{}</b> (<code>{}</code>): {}</div><ul>",
                w.original_word, w.hyphenated, w.literal_gloss
            );
            for seg in &w.segments {
                if seg.is_connecting_vowel {
                    continue;
                }
                let _ = write!(
                    back,
                    "<li><code>{}</code> [{} — {}]: <b>{}</b></li>",
                    seg.surface_slice, seg.origin_label, seg.classical_lemma, seg.gloss
                );
            }
            back.push_str("</ul>");
        }
        let _ = writeln!(
            out,
            "{}\t{}\t{}",
            front.replace('\t', " "),
            back.replace('\t', " "),
            clean_tag
        );
    }
    out
}

/// Exports decoded terms as a clean Markdown study sheet with morpheme breakdowns and cognates.
pub fn export_flashcards_markdown(phrases: &[DecodedPhrase]) -> String {
    let mut md = String::from("# OSABT Biological & Anatomical Etymology Study Sheet\n\n");
    for p in phrases {
        let _ = writeln!(md, "## *{}*\n", p.raw_input);
        let _ = writeln!(md, "- **Combined Gloss:** {}\n", p.combined_gloss);
        md.push_str("| Word | Morpheme | Origin | Classical Lemma | Meaning | Cognates |\n");
        md.push_str("| :--- | :--- | :--- | :--- | :--- | :--- |\n");
        for w in &p.words {
            for seg in &w.segments {
                let cogs = if seg.cognates.is_empty() {
                    String::from("—")
                } else {
                    seg.cognates.join(", ")
                };
                let _ = writeln!(
                    md,
                    "| {} | `{}` | {} | {} | {} | {} |",
                    w.original_word,
                    seg.surface_slice,
                    seg.origin_label,
                    seg.classical_lemma,
                    seg.gloss,
                    cogs
                );
            }
        }
        md.push('\n');
    }
    md
}

/// Built-in showcase list of classic taxonomic & anatomical compounds for demos and study decks.
pub fn demo_etymology_terms() -> Vec<&'static str> {
    vec![
        "Drosophila melanogaster",
        "Ornithorhynchus anatinus",
        "Pseudotsuga menziesii",
        "Actinopterygii",
        "Sternocleidomastoideus",
        "Saccharomyces cerevisiae",
        "Enterobacteriaceae",
        "Rhododendron macrophyllum",
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_classic_binomials_anatomy_and_iczn_icnafp_suffixes() {
        let dros = decode_biological_phrase("Drosophila melanogaster");
        assert_eq!(dros.words.len(), 2);
        assert!(dros.words[0].literal_gloss.contains("dew"));
        assert!(dros.words[1].literal_gloss.contains("black"));
        assert!(dros.words[1].literal_gloss.contains("stomach"));

        let ps = decode_biological_phrase("Pseudotsuga menziesii");
        assert!(ps.words[0].literal_gloss.contains("false"));
        assert!(ps.words[0].literal_gloss.contains("hemlock"));
        assert!(ps.words[1].literal_gloss.contains("Menzies"));

        let scm = decode_single_word("Sternocleidomastoideus");
        assert!(scm.literal_gloss.contains("sternum") || scm.literal_gloss.contains("breastbone"));
        assert!(scm.literal_gloss.contains("clavicle"));

        let anki = export_anki_flashcards_tsv(&[dros], "OSABT_Biology");
        assert!(anki.starts_with("#separator:tab"));
    }
}
