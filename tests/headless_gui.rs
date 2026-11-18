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

#![cfg(feature = "gui")]

use eframe::egui;
use osabt::brand::draw_standalone_footer;
use osabt::ecosound_gui::EcosoundStudioApp;
use osabt::flow_gui::FlowGuiApp;
use osabt::flux_gui::FluxGuiApp;
use osabt::forge3d_gui::Forge3dStudioApp;
use osabt::kinetics_gui::KineticsStudioApp;
use osabt::morpho_gui::MorphoStudioApp;
use osabt::plasmid_gui::PlasmidStudioApp;
use osabt::popgen_gui::PopGenGuiApp;
use osabt::workbench_gui::{OsabtWorkbenchApp, ToolId};

fn run_frames(
    ctx: &egui::Context,
    size: [f32; 2],
    frames: usize,
    mut draw: impl FnMut(&egui::Context, usize),
) {
    for frame_idx in 0..frames {
        let mut raw_input = egui::RawInput::default();
        raw_input.screen_rect = Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(size[0], size[1]),
        ));
        raw_input
            .viewports
            .entry(egui::ViewportId::ROOT)
            .and_modify(|v| v.native_pixels_per_point = Some(1.0));
        raw_input.time = Some(frame_idx as f64 / 60.0);
        ctx.set_pixels_per_point(1.0);
        let out = ctx.run(raw_input, |c| draw(c, frame_idx));
        assert!(
            !out.shapes.is_empty(),
            "headless egui frame {frame_idx} emitted zero clipped shapes"
        );
    }
}

#[test]
fn workbench_renders_launcher_and_all_20_tools_without_panic() {
    const ALL_TOOLS: [ToolId; 20] = [
        ToolId::Bench,
        ToolId::Plasmid,
        ToolId::Dwc,
        ToolId::Phylo,
        ToolId::Epi,
        ToolId::Lexi,
        ToolId::Forge3d,
        ToolId::Morpho,
        ToolId::Plate,
        ToolId::Ecosound,
        ToolId::Infoevo,
        ToolId::Qpcr,
        ToolId::Kinetics,
        ToolId::Edna,
        ToolId::Wildlife,
        ToolId::Popgen,
        ToolId::Crispr,
        ToolId::Flow,
        ToolId::Survival,
        ToolId::Flux,
    ];
    let mut app = OsabtWorkbenchApp::new();
    let ctx = egui::Context::default();

    // Home Launcher at desktop & minimum sizes
    run_frames(&ctx, [1480.0, 920.0], 2, |c, _| app.draw_panels(c));
    run_frames(&ctx, [960.0, 640.0], 1, |c, _| app.draw_panels(c));

    // All 20 tools (12 workbench panels + 8 dedicated studios, including ToolId::Plate)
    for tool in ALL_TOOLS {
        app.select_tool(tool);
        run_frames(&ctx, [1480.0, 920.0], 2, |c, _| app.draw_panels(c));
        run_frames(&ctx, [960.0, 640.0], 1, |c, _| app.draw_panels(c));
    }
}

#[test]
fn all_8_standalone_studios_render_with_footer_and_about_modal() {
    let mut plasmid = PlasmidStudioApp::new(None);
    let mut forge3d = Forge3dStudioApp::new(None);
    let mut morpho = MorphoStudioApp::new(None);
    let mut ecosound = EcosoundStudioApp::new(None);
    let mut kinetics = KineticsStudioApp::new(None);
    let mut popgen = PopGenGuiApp::new(None, None).expect("demo popgen");
    let mut flow = FlowGuiApp::new(None).expect("demo flow");
    let mut flux = FluxGuiApp::new(None, None).expect("demo flux");

    let studios: [(&str, &mut dyn FnMut(&egui::Context)); 8] = [
        ("plasmid", &mut |c| plasmid.draw_panels(c)),
        ("forge3d", &mut |c| forge3d.draw_panels(c)),
        ("morpho", &mut |c| morpho.draw_panels(c)),
        ("ecosound", &mut |c| ecosound.draw_panels(c)),
        ("kinetics", &mut |c| kinetics.draw_panels(c)),
        ("popgen", &mut |c| popgen.draw_panels(c)),
        ("flow", &mut |c| flow.draw_panels(c)),
        ("flux", &mut |c| flux.draw_panels(c)),
    ];

    for (_name, draw_fn) in studios {
        let ctx = egui::Context::default();
        for about_open in [false, true] {
            let mut open = about_open;
            run_frames(&ctx, [1280.0, 800.0], 2, |c, _| {
                draw_standalone_footer(c, &mut open);
                draw_fn(c);
            });
            run_frames(&ctx, [960.0, 640.0], 1, |c, _| {
                draw_standalone_footer(c, &mut open);
                draw_fn(c);
            });
        }
    }
}
