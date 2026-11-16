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

//! Shared zero-dependency vector `.svg` canvas, coordinate formatting, and XML/SVG escaping
//! primitives used across all 20 OSABT scientific visualizers.

use std::fmt::Write as _;

/// Escapes reserved XML/SVG text and attribute characters (`&`, `<`, `>`, `"`, `'`).
pub fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

/// Appends an opening `<svg xmlns=...>` root element and optional solid background `<rect>`
/// to an existing SVG buffer.
pub fn write_svg_open(svg: &mut String, width: f64, height: f64, background_fill: Option<&str>) {
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    if let Some(fill) = background_fill {
        let _ = writeln!(
            svg,
            "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"{}\"/>",
            xml_escape(fill)
        );
    }
}

/// Allocates a new SVG string buffer with the given initial capacity, opening `<svg>` tag,
/// and optional full-bleed background `<rect>`.
pub fn begin_svg_canvas(
    capacity: usize,
    width: f64,
    height: f64,
    background_fill: Option<&str>,
) -> String {
    let mut svg = String::with_capacity(capacity);
    write_svg_open(&mut svg, width, height, background_fill);
    svg
}

/// Formats `(x, y)` tuples into a space-delimited SVG `<polyline points="...">` attribute string.
pub fn format_polyline_points(points: &[(f64, f64)]) -> String {
    let mut poly = String::with_capacity(points.len().saturating_mul(16));
    for (idx, &(x, y)) in points.iter().enumerate() {
        if idx > 0 {
            poly.push(' ');
        }
        let _ = write!(poly, "{x:.2},{y:.2}");
    }
    poly
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_xml_opens_svg_and_formats_polyline() {
        assert_eq!(
            xml_escape("A&B <C> \"D\" 'E'"),
            "A&amp;B &lt;C&gt; &quot;D&quot; &apos;E&apos;"
        );
        let svg = begin_svg_canvas(256, 800.0, 420.0, Some("#ffffff"));
        assert!(svg.contains("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 800 420\" width=\"800\" height=\"420\">"));
        assert!(svg.contains("<rect x=\"0\" y=\"0\" width=\"800\" height=\"420\" fill=\"#ffffff\"/>"));
        assert_eq!(format_polyline_points(&[(1.0, 2.5), (3.25, 4.125)]), "1.00,2.50 3.25,4.12");
    }
}
