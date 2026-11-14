use osabt::infoevo_core::{
    compute_alignment_information, demo_ecoli_promoter_and_stemloop_fasta, parse_fasta_alignment,
};
use osabt::phylo_core::{compare_robinson_foulds, parse_phylo_trees};
use osabt::plate_curve::{
    demo_measured_96well_plate, evaluate_plate_qc_and_normalize, fit_all_plate_compounds,
};
use osabt::plate_layout::{
    export_plate_pipetting_csv, generate_randomized_plate, render_plate_map_svg, PlatePlanConfig,
};
use osabt::popgen_core::PopSequenceRecord;
use osabt::popgen_network::{build_median_joining_network, render_popgen_network_svg};

#[test]
fn plate_design_qc_4pl_and_exports_end_to_end() {
    let cfg = PlatePlanConfig::default();
    let planned = generate_randomized_plate(&cfg).expect("plate planning succeeds");
    assert_eq!(planned.wells.len(), 96);
    let csv = export_plate_pipetting_csv(&planned);
    let svg = render_plate_map_svg(&planned);
    assert!(csv.contains("well_id,row,col,role"));
    assert!(svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>"));

    let mut measured = demo_measured_96well_plate();
    let qc = evaluate_plate_qc_and_normalize(&mut measured).expect("QC succeeds");
    assert!(qc.z_prime_factor > 0.5);
    let fits = fit_all_plate_compounds(&measured);
    assert_eq!(fits.len(), 3);
    assert!(fits.iter().all(|f| f.relative_ec50_nm > 0.0 && f.r_squared > 0.90));
}

#[test]
fn fasta_to_infoevo_popgen_network_and_phylo_newick_round_trip() {
    let aln = parse_fasta_alignment("demo", demo_ecoli_promoter_and_stemloop_fasta())
        .expect("FASTA parses");
    let info = compute_alignment_information(&aln, true, None);
    assert!(info.alignment_length > 0 && info.total_r_sequence_bits > 0.0);

    let pop_records: Vec<PopSequenceRecord> = aln
        .records
        .iter()
        .enumerate()
        .map(|(idx, rec)| PopSequenceRecord {
            sample_id: rec.id.clone(),
            population: format!("Clade_{}", idx % 2),
            sequence: String::from_utf8_lossy(&rec.residues).into_owned(),
        })
        .collect();
    let net = build_median_joining_network(&pop_records, 0).expect("MJN succeeds");
    let net_svg = render_popgen_network_svg(&net, "E. coli MJN");
    assert!(net_svg.starts_with("<svg") && net_svg.trim_end().ends_with("</svg>"));

    let nwk = "(((Taxon_A:0.12,Taxon_B:0.15):0.08,Taxon_C:0.20):0.05,(Taxon_D:0.18,Taxon_E:0.22):0.10);";
    let trees = parse_phylo_trees("t", nwk).expect("Newick parses");
    let tree_a = &trees[0];
    let mut tree_b = tree_a.clone();
    tree_b.reroot_by_outgroup("Taxon_E").expect("reroot succeeds");
    let reparsed = parse_phylo_trees("b", &tree_b.to_newick()).expect("round-trip Newick parses");
    let rf = compare_robinson_foulds(tree_a, &reparsed[0]).expect("RF comparison succeeds");
    assert_eq!(rf.rf_distance, 0);
}

#[cfg(feature = "cli")]
#[test]
fn cli_plate_plan_and_phylo_svg_generate_outputs_in_temp_dir() {
    use clap::Parser;
    use osabt::phylo_cli::{run_phylo_cli, PhyloArgs};
    use osabt::plate_cli::{run_plate_cli, PlateArgs};

    #[derive(Parser)]
    struct TestPlateCli {
        #[command(flatten)]
        args: PlateArgs,
    }
    #[derive(Parser)]
    struct TestPhyloCli {
        #[command(flatten)]
        args: PhyloArgs,
    }

    let tmp = std::env::temp_dir().join(format!("osabt_it_{}", std::process::id()));
    std::fs::create_dir_all(&tmp).expect("create temp dir");
    let prefix = tmp.join("screen_p1");
    let prefix_str = prefix.to_string_lossy().into_owned();

    let parsed_plate = TestPlateCli::try_parse_from([
        "plate",
        "plan",
        "--wells",
        "96",
        "--moat",
        "1",
        "--seed",
        "42",
        "--out-prefix",
        &prefix_str,
    ])
    .expect("parse plate plan CLI args");
    run_plate_cli(parsed_plate.args).expect("run_plate_cli succeeds");

    let csv_path = tmp.join("screen_p1_layout.csv");
    let md_path = tmp.join("screen_p1_protocol.md");
    let svg_path = tmp.join("screen_p1_map.svg");
    assert!(csv_path.exists() && md_path.exists() && svg_path.exists());

    let phylo_svg = tmp.join("tree.svg");
    let phylo_svg_str = phylo_svg.to_string_lossy().into_owned();
    let parsed_phylo = TestPhyloCli::try_parse_from(["phylo", "svg", "-o", &phylo_svg_str])
        .expect("parse phylo svg CLI args");
    run_phylo_cli(parsed_phylo.args).expect("run_phylo_cli succeeds");
    let tree_svg_text = std::fs::read_to_string(&phylo_svg).expect("read tree svg");
    assert!(tree_svg_text.starts_with("<svg"));

    let _ = std::fs::remove_dir_all(&tmp);
}
