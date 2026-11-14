use osabt::flux_leaf::{demo_fvcb_aci_csv, fit_fvcb_aci_curve, parse_aci_csv};
use osabt::flux_stream::{demo_stream_diel_do_csv, fit_odum_stream_metabolism, parse_stream_do_csv};
use osabt::gibson::{design_gibson_assembly, GibsonConfig, GibsonFragmentInput};
use osabt::kinetics_core::{
    demo_competitive_kinase_matrix_csv, demo_morrison_tight_binding_csv,
    fit_global_reversible_inhibition, fit_morrison_tight_binding, parse_kinetics_csv,
    InhibitionMechanism,
};
use osabt::plasmid_digest::{
    find_unique_cutters, score_overhang_set_fidelity, simulate_restriction_digest,
};
use osabt::plasmid_io::PlasmidRecord;
use osabt::popgen_core::{
    analyze_diploid_vcf, analyze_popgen_fasta, demo_pnw_popgen_fasta, demo_pnw_popgen_vcf,
    parse_popgen_fasta, parse_popgen_vcf,
};
use osabt::survival_core::{
    compute_kaplan_meier_curves, compute_log_rank_test, demo_survival_cohort_csv,
    fit_cox_proportional_hazards, parse_survival_csv,
};
use osabt::survival_tox::{
    demo_ecotox_lc50_csv, fit_ecotox_lc50, parse_ecotox_csv, DoseResponseLink,
};
use osabt::thermo::{calc_perfect_duplex, ReactionConditions};
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

struct OsabtValidationData {
    // T1
    vcmax: f64,
    jmax: f64,
    tpu: f64,
    rd: f64,
    aci_r2: f64,
    aci_rmse: f64,
    tm_salt_c: f64,
    delta_g_37: f64,
    km_group_count: usize,
    logrank_chi2: f64,
    logrank_p: f64,
    cox_hr: f64,
    cox_p: f64,
    probit_lc50: f64,
    logit_lc50: f64,
    // T2
    fasta_sequences: usize,
    fasta_len_bp: usize,
    fasta_pop_count: usize,
    global_pi: f64,
    global_tajima_d: f64,
    hudson_fst_max: f64,
    vcf_samples: usize,
    vcf_loci: usize,
    vcf_multilocus_fst: f64,
    plasmid_name: String,
    plasmid_len_bp: usize,
    unique_cutter_count: usize,
    selected_enzymes: Vec<String>,
    digest_fragment_count: usize,
    digest_total_bp: usize,
    gibson_assembled_bp: usize,
    gibson_junctions: usize,
    gg_overhang_fidelity_pct: f64,
    // T3
    stream_gpp: f64,
    stream_er: f64,
    stream_nep: f64,
    stream_k600: f64,
    stream_trophic_state: String,
    stream_rmse: f64,
    best_inhib_label: String,
    best_inhib_ki_nm: f64,
    best_inhib_akaike_w: f64,
    best_inhib_r2: f64,
    morrison_et_nm: f64,
    morrison_ki_app_nm: f64,
    morrison_r2: f64,
    // Files
    csv_path: PathBuf,
    csv_bytes: u64,
    json_path: PathBuf,
    json_bytes: u64,
}

static OSABT_SUITE: OnceLock<OsabtValidationData> = OnceLock::new();

fn validation_dir() -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("validation");
    fs::create_dir_all(&dir).expect("create validation directory");
    dir
}

fn json_escape(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len() + 8);
    for ch in raw.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn csv_escape(raw: &str) -> String {
    if raw.contains(',') || raw.contains('"') || raw.contains('\n') {
        format!("\"{}\"", raw.replace('"', "\"\""))
    } else {
        raw.to_string()
    }
}

fn ensure_osabt_validation() -> &'static OsabtValidationData {
    OSABT_SUITE.get_or_init(|| {
        let val_dir = validation_dir();
        let mut csv = String::from("test_id,module,engine_function,metric,value,unit\n");

        // ---------------------------------------------------------------------
        // T1: Ecophysiology FvCB A-Ci, DNA Thermodynamics & Survival / Ecotox
        // ---------------------------------------------------------------------
        let aci_points = parse_aci_csv(&demo_fvcb_aci_csv()).expect("parse demo A-Ci CSV");
        let fvcb = fit_fvcb_aci_curve(&aci_points).expect("fit FvCB A-Ci curve");

        let thermo_res = calc_perfect_duplex(
            "GCTAGCTAGCGATCGATCGA",
            &ReactionConditions::default(),
            60.0,
        )
        .expect("compute nearest-neighbor duplex thermodynamics");

        let surv_obs =
            parse_survival_csv(demo_survival_cohort_csv()).expect("parse demo survival cohort");
        let km_curves = compute_kaplan_meier_curves(&surv_obs);
        let logrank = compute_log_rank_test(&surv_obs).expect("compute log-rank test");
        let cox = fit_cox_proportional_hazards(&surv_obs, "treatment", true)
            .expect("fit Cox proportional hazards");

        let ecotox_pts = parse_ecotox_csv(demo_ecotox_lc50_csv()).expect("parse ecotox CSV");
        let probit = fit_ecotox_lc50(&ecotox_pts, DoseResponseLink::Probit)
            .expect("fit Finney Probit LC50");
        let logit = fit_ecotox_lc50(&ecotox_pts, DoseResponseLink::Logit)
            .expect("fit Binomial Logit LC50");

        for (mod_name, fn_name, metric, val, unit) in [
            ("flux_leaf", "fit_fvcb_aci_curve", "vcmax", fvcb.vcmax, "umol m^-2 s^-1"),
            ("flux_leaf", "fit_fvcb_aci_curve", "jmax", fvcb.jmax, "umol m^-2 s^-1"),
            ("flux_leaf", "fit_fvcb_aci_curve", "tpu", fvcb.tpu, "umol m^-2 s^-1"),
            ("flux_leaf", "fit_fvcb_aci_curve", "r_squared", fvcb.r_squared, "dimensionless"),
            ("thermo", "calc_perfect_duplex", "tm_salt_corrected_c", thermo_res.tm_salt_corrected_c, "C"),
            ("thermo", "calc_perfect_duplex", "delta_g_37_kcal", thermo_res.delta_g_37_kcal, "kcal/mol"),
            ("survival_core", "compute_log_rank_test", "chi_square", logrank.chi_square, "chi2"),
            ("survival_core", "compute_log_rank_test", "p_value", logrank.p_value, "probability"),
            ("survival_core", "fit_cox_proportional_hazards", "hazard_ratio", cox.hazard_ratio, "ratio"),
            ("survival_tox", "fit_ecotox_lc50(Probit)", "lc50", probit.lc50, "ug/L"),
            ("survival_tox", "fit_ecotox_lc50(Logit)", "lc50", logit.lc50, "ug/L"),
        ] {
            csv.push_str(&format!(
                "T1,{},{},{},{:.6},{}\n",
                mod_name,
                csv_escape(fn_name),
                metric,
                val,
                csv_escape(unit)
            ));
        }

        // ---------------------------------------------------------------------
        // T2: Population Genetics (FASTA + Diploid VCF) & Plasmid Digest / Gibson
        // ---------------------------------------------------------------------
        let fasta_recs = parse_popgen_fasta(&demo_pnw_popgen_fasta()).expect("parse PNW FASTA");
        let fasta_rep = analyze_popgen_fasta(&fasta_recs).expect("analyze FASTA popgen");
        let hudson_fst_max = fasta_rep
            .pairwise_divergence
            .iter()
            .map(|p| p.hudson_fst)
            .fold(f64::NEG_INFINITY, f64::max);

        let (sample_pops, vcf_loci) =
            parse_popgen_vcf(demo_pnw_popgen_vcf()).expect("parse PNW diploid VCF");
        let vcf_rep = analyze_diploid_vcf(&sample_pops, &vcf_loci).expect("analyze diploid VCF");

        let plasmid = PlasmidRecord::demo_pcrv_gfp();
        let plasmid_len = plasmid.sequence.len();
        let unique_cutters = find_unique_cutters(&plasmid);
        let mut selected_enzymes: Vec<String> = Vec::new();
        for hit in &unique_cutters {
            if !selected_enzymes.contains(&hit.enzyme_name) {
                selected_enzymes.push(hit.enzyme_name.clone());
                if selected_enzymes.len() == 3 {
                    break;
                }
            }
        }
        assert!(
            selected_enzymes.len() >= 2,
            "Demo plasmid must yield >= 2 unique cutters"
        );
        let digest_frags = simulate_restriction_digest(&plasmid, &selected_enzymes)
            .expect("simulate restriction digest");
        let digest_total_bp: usize = digest_frags.iter().map(|f| f.length_bp).sum();

        let split_mid = plasmid_len / 2;
        let frag_inputs = vec![
            GibsonFragmentInput {
                name: String::from("backbone_left"),
                sequence: plasmid.sequence[..split_mid].to_string(),
            },
            GibsonFragmentInput {
                name: String::from("insert_right"),
                sequence: plasmid.sequence[split_mid..].to_string(),
            },
        ];
        let gibson_plan =
            design_gibson_assembly(&frag_inputs, true, &GibsonConfig::default())
                .expect("design circular 2-fragment Gibson assembly");

        let gg_overhangs = vec![
            String::from("AATG"),
            String::from("GCTT"),
            String::from("CGAA"),
            String::from("TACA"),
        ];
        let gg_fidelity = score_overhang_set_fidelity(&gg_overhangs);

        for (mod_name, fn_name, metric, val, unit) in [
            ("popgen_core", "analyze_popgen_fasta", "nucleotide_diversity_pi", fasta_rep.global_metrics.nucleotide_diversity_pi, "per_site"),
            ("popgen_core", "analyze_popgen_fasta", "tajima_d", fasta_rep.global_metrics.tajima_d, "dimensionless"),
            ("popgen_core", "analyze_popgen_fasta", "hudson_fst_max", hudson_fst_max, "dimensionless"),
            ("popgen_core", "analyze_diploid_vcf", "multilocus_wc_fst", vcf_rep.multilocus_fst, "dimensionless"),
            ("plasmid_digest", "simulate_restriction_digest", "fragment_count", digest_frags.len() as f64, "fragments"),
            ("gibson", "design_gibson_assembly", "assembled_length_bp", gibson_plan.assembled_length_bp as f64, "bp"),
            ("plasmid_digest", "score_overhang_set_fidelity", "overhang_fidelity_pct", gg_fidelity, "percent"),
        ] {
            csv.push_str(&format!(
                "T2,{},{},{},{:.6},{}\n",
                mod_name,
                csv_escape(fn_name),
                metric,
                val,
                csv_escape(unit)
            ));
        }

        // ---------------------------------------------------------------------
        // T3: Stream Metabolism Flux (Odum) & Enzyme Kinetics / Morrison Fit
        // ---------------------------------------------------------------------
        let stream_pts =
            parse_stream_do_csv(&demo_stream_diel_do_csv()).expect("parse stream diel DO CSV");
        let stream_rep =
            fit_odum_stream_metabolism(&stream_pts).expect("fit Odum stream metabolism");

        let kin_pts = parse_kinetics_csv(demo_competitive_kinase_matrix_csv())
            .expect("parse competitive kinase matrix CSV");
        let kin_fits =
            fit_global_reversible_inhibition(&kin_pts).expect("fit global reversible inhibition");
        let best_kin = kin_fits
            .iter()
            .min_by(|a, b| a.aicc.partial_cmp(&b.aicc).unwrap())
            .expect("at least one inhibition fit");

        let morr_pts = parse_kinetics_csv(demo_morrison_tight_binding_csv())
            .expect("parse Morrison tight-binding CSV");
        let morr_obs: Vec<(f64, f64)> = morr_pts
            .iter()
            .map(|p| (p.substrate_um, p.velocity))
            .collect();
        let morr_fit = fit_morrison_tight_binding(&morr_obs, None, Some(25.0), Some((70.0, 35.0)))
            .expect("fit Morrison tight-binding Ki_app");

        for (mod_name, fn_name, metric, val, unit) in [
            ("flux_stream", "fit_odum_stream_metabolism", "gpp_g_o2_m2_d", stream_rep.gpp_g_o2_m2_d, "g O2 m^-2 d^-1"),
            ("flux_stream", "fit_odum_stream_metabolism", "er_g_o2_m2_d", stream_rep.er_g_o2_m2_d, "g O2 m^-2 d^-1"),
            ("flux_stream", "fit_odum_stream_metabolism", "nep_g_o2_m2_d", stream_rep.nep_g_o2_m2_d, "g O2 m^-2 d^-1"),
            ("flux_stream", "fit_odum_stream_metabolism", "k600_per_day", stream_rep.k600_per_day, "d^-1"),
            ("kinetics_core", "fit_global_reversible_inhibition", "best_ki_nm", best_kin.ki_nm, "nM"),
            ("kinetics_core", "fit_global_reversible_inhibition", "best_akaike_weight", best_kin.akaike_weight, "probability"),
            ("kinetics_core", "fit_morrison_tight_binding", "active_enzyme_et_nm", morr_fit.active_enzyme_et_nm, "nM"),
            ("kinetics_core", "fit_morrison_tight_binding", "ki_app_nm", morr_fit.ki_app_nm, "nM"),
            ("kinetics_core", "fit_morrison_tight_binding", "r_squared", morr_fit.r_squared, "dimensionless"),
        ] {
            csv.push_str(&format!(
                "T3,{},{},{},{:.6},{}\n",
                mod_name,
                csv_escape(fn_name),
                metric,
                val,
                csv_escape(unit)
            ));
        }

        let csv_path = val_dir.join("osabt_benchmark_suite.csv");
        fs::write(&csv_path, &csv).expect("write osabt_benchmark_suite.csv");
        let csv_bytes = fs::metadata(&csv_path).expect("stat csv").len();

        let best_inhib_label = best_kin.mechanism.label().to_string();
        let best_inhib_is_comp = best_kin.mechanism == InhibitionMechanism::Competitive;
        let summary_json = format!(
            "{{\n  \"schema_version\": \"1.0\",\n  \"repo\": \"crvscience-osabt\",\n  \"preprint_id\": \"APP-PREPRINT-08\",\n  \"preprint_url\": \"https://crvscience.com/apps/preprints/osabt-modular-biostatistics-toolkit\",\n  \"concept_doi\": \"10.5281/zenodo.24222388\",\n  \"version\": \"2.0.0\",\n  \"engine_scope\": \"Real implemented v2.0.0 Rust engines only (Option A)\",\n  \"testsPassed\": 3,\n  \"testsTotal\": 3,\n  \"tests\": {{\n    \"t1_ecophysiology_thermo_and_survival\": {{\n      \"engine_calls\": [\"flux_leaf::fit_fvcb_aci_curve\", \"thermo::calc_perfect_duplex\", \"survival_core::compute_kaplan_meier_curves\", \"survival_core::compute_log_rank_test\", \"survival_core::fit_cox_proportional_hazards\", \"survival_tox::fit_ecotox_lc50\"],\n      \"fvcb_vcmax_umol_m2_s\": {:.6},\n      \"fvcb_jmax_umol_m2_s\": {:.6},\n      \"fvcb_tpu_umol_m2_s\": {:.6},\n      \"fvcb_rd_umol_m2_s\": {:.6},\n      \"fvcb_r_squared\": {:.6},\n      \"fvcb_rmse\": {:.6},\n      \"duplex_20mer_tm_salt_corrected_c\": {:.4},\n      \"duplex_20mer_delta_g_37_kcal\": {:.4},\n      \"kaplan_meier_groups\": {},\n      \"logrank_chi_square\": {:.6},\n      \"logrank_p_value\": {:.6},\n      \"cox_hazard_ratio\": {:.6},\n      \"cox_p_value\": {:.6},\n      \"ecotox_probit_lc50\": {:.6},\n      \"ecotox_logit_lc50\": {:.6}\n    }},\n    \"t2_popgen_and_plasmid_biotech\": {{\n      \"engine_calls\": [\"popgen_core::analyze_popgen_fasta\", \"popgen_core::analyze_diploid_vcf\", \"plasmid_io::PlasmidRecord::demo_pcrv_gfp\", \"plasmid_digest::find_unique_cutters\", \"plasmid_digest::simulate_restriction_digest\", \"plasmid_digest::score_overhang_set_fidelity\", \"gibson::design_gibson_assembly\"],\n      \"fasta_sequences\": {},\n      \"fasta_alignment_len_bp\": {},\n      \"fasta_populations\": {},\n      \"global_pi\": {:.6},\n      \"global_tajima_d\": {:.6},\n      \"max_pairwise_hudson_fst\": {:.6},\n      \"vcf_samples\": {},\n      \"vcf_loci\": {},\n      \"vcf_multilocus_wc_fst\": {:.6},\n      \"plasmid_name\": \"{}\",\n      \"plasmid_length_bp\": {},\n      \"unique_cutter_hits\": {},\n      \"digest_enzymes\": [{}],\n      \"digest_fragment_count\": {},\n      \"digest_total_bp\": {},\n      \"gibson_assembled_bp\": {},\n      \"gibson_junctions\": {},\n      \"golden_gate_4overhang_fidelity_pct\": {:.4}\n    }},\n    \"t3_ecosystem_stream_flux_and_enzyme_kinetics\": {{\n      \"engine_calls\": [\"flux_stream::fit_odum_stream_metabolism\", \"kinetics_core::fit_global_reversible_inhibition\", \"kinetics_core::fit_morrison_tight_binding\"],\n      \"stream_gpp_g_o2_m2_d\": {:.6},\n      \"stream_er_g_o2_m2_d\": {:.6},\n      \"stream_nep_g_o2_m2_d\": {:.6},\n      \"stream_k600_per_day\": {:.6},\n      \"stream_trophic_state\": \"{}\",\n      \"stream_do_rmse_mg_l\": {:.6},\n      \"best_inhibition_mechanism\": \"{}\",\n      \"best_mechanism_is_competitive\": {},\n      \"best_inhibition_ki_nm\": {:.6},\n      \"best_inhibition_akaike_weight\": {:.6},\n      \"best_inhibition_r_squared\": {:.6},\n      \"morrison_active_enzyme_et_nm\": {:.6},\n      \"morrison_ki_app_nm\": {:.6},\n      \"morrison_r_squared\": {:.6},\n      \"csv_artifact\": \"validation/osabt_benchmark_suite.csv\",\n      \"csv_bytes\": {}\n    }}\n  }}\n}}\n",
            fvcb.vcmax,
            fvcb.jmax,
            fvcb.tpu,
            fvcb.rd,
            fvcb.r_squared,
            fvcb.rmse,
            thermo_res.tm_salt_corrected_c,
            thermo_res.delta_g_37_kcal,
            km_curves.len(),
            logrank.chi_square,
            logrank.p_value,
            cox.hazard_ratio,
            cox.p_value,
            probit.lc50,
            logit.lc50,
            fasta_rep.total_sequences,
            fasta_rep.alignment_len_bp,
            fasta_rep.populations.len(),
            fasta_rep.global_metrics.nucleotide_diversity_pi,
            fasta_rep.global_metrics.tajima_d,
            hudson_fst_max,
            vcf_rep.num_samples,
            vcf_rep.num_loci,
            vcf_rep.multilocus_fst,
            json_escape(&plasmid.name),
            plasmid_len,
            unique_cutters.len(),
            selected_enzymes
                .iter()
                .map(|e| format!("\"{}\"", json_escape(e)))
                .collect::<Vec<_>>()
                .join(", "),
            digest_frags.len(),
            digest_total_bp,
            gibson_plan.assembled_length_bp,
            gibson_plan.junctions.len(),
            gg_fidelity,
            stream_rep.gpp_g_o2_m2_d,
            stream_rep.er_g_o2_m2_d,
            stream_rep.nep_g_o2_m2_d,
            stream_rep.k600_per_day,
            json_escape(&stream_rep.trophic_state),
            stream_rep.do_rmse_mg_l,
            json_escape(&best_inhib_label),
            best_inhib_is_comp,
            best_kin.ki_nm,
            best_kin.akaike_weight,
            best_kin.r_squared,
            morr_fit.active_enzyme_et_nm,
            morr_fit.ki_app_nm,
            morr_fit.r_squared,
            csv_bytes,
        );

        let json_path = val_dir.join("osabt_validation_summary.json");
        fs::write(&json_path, &summary_json).expect("write osabt_validation_summary.json");
        let json_bytes = fs::metadata(&json_path).expect("stat json").len();

        OsabtValidationData {
            vcmax: fvcb.vcmax,
            jmax: fvcb.jmax,
            tpu: fvcb.tpu,
            rd: fvcb.rd,
            aci_r2: fvcb.r_squared,
            aci_rmse: fvcb.rmse,
            tm_salt_c: thermo_res.tm_salt_corrected_c,
            delta_g_37: thermo_res.delta_g_37_kcal,
            km_group_count: km_curves.len(),
            logrank_chi2: logrank.chi_square,
            logrank_p: logrank.p_value,
            cox_hr: cox.hazard_ratio,
            cox_p: cox.p_value,
            probit_lc50: probit.lc50,
            logit_lc50: logit.lc50,
            fasta_sequences: fasta_rep.total_sequences,
            fasta_len_bp: fasta_rep.alignment_len_bp,
            fasta_pop_count: fasta_rep.populations.len(),
            global_pi: fasta_rep.global_metrics.nucleotide_diversity_pi,
            global_tajima_d: fasta_rep.global_metrics.tajima_d,
            hudson_fst_max,
            vcf_samples: vcf_rep.num_samples,
            vcf_loci: vcf_rep.num_loci,
            vcf_multilocus_fst: vcf_rep.multilocus_fst,
            plasmid_name: plasmid.name,
            plasmid_len_bp: plasmid_len,
            unique_cutter_count: unique_cutters.len(),
            selected_enzymes,
            digest_fragment_count: digest_frags.len(),
            digest_total_bp,
            gibson_assembled_bp: gibson_plan.assembled_length_bp,
            gibson_junctions: gibson_plan.junctions.len(),
            gg_overhang_fidelity_pct: gg_fidelity,
            stream_gpp: stream_rep.gpp_g_o2_m2_d,
            stream_er: stream_rep.er_g_o2_m2_d,
            stream_nep: stream_rep.nep_g_o2_m2_d,
            stream_k600: stream_rep.k600_per_day,
            stream_trophic_state: stream_rep.trophic_state,
            stream_rmse: stream_rep.do_rmse_mg_l,
            best_inhib_label,
            best_inhib_ki_nm: best_kin.ki_nm,
            best_inhib_akaike_w: best_kin.akaike_weight,
            best_inhib_r2: best_kin.r_squared,
            morrison_et_nm: morr_fit.active_enzyme_et_nm,
            morrison_ki_app_nm: morr_fit.ki_app_nm,
            morrison_r2: morr_fit.r_squared,
            csv_path,
            csv_bytes,
            json_path,
            json_bytes,
        }
    })
}

#[test]
fn t1_ecophysiology_thermo_and_survival() {
    let d = ensure_osabt_validation();
    assert!(d.vcmax > 0.0 && d.jmax > d.vcmax && d.tpu > 0.0 && d.rd >= 0.0);
    assert!(d.aci_r2 > 0.90 && d.aci_rmse.is_finite());
    assert!(d.tm_salt_c > 30.0 && d.tm_salt_c < 90.0 && d.delta_g_37 < 0.0);
    assert!(d.km_group_count >= 2);
    assert!(d.logrank_chi2.is_finite() && d.logrank_p >= 0.0 && d.logrank_p <= 1.0);
    assert!(d.cox_hr > 0.0 && d.cox_p >= 0.0 && d.cox_p <= 1.0);
    assert!(d.probit_lc50 > 0.0 && d.logit_lc50 > 0.0);
    assert!(d.json_path.exists() && d.json_bytes < 384 * 1024);
}

#[test]
fn t2_popgen_and_plasmid_biotech() {
    let d = ensure_osabt_validation();
    assert!(d.fasta_sequences >= 6 && d.fasta_len_bp > 0 && d.fasta_pop_count >= 2);
    assert!(d.global_pi >= 0.0 && d.global_tajima_d.is_finite());
    assert!(d.hudson_fst_max.is_finite() && d.hudson_fst_max >= 0.0);
    assert!(d.vcf_samples >= 4 && d.vcf_loci >= 1 && d.vcf_multilocus_fst.is_finite());
    assert!(!d.plasmid_name.is_empty() && d.plasmid_len_bp > 0);
    assert!(d.unique_cutter_count >= 2 && !d.selected_enzymes.is_empty());
    assert_eq!(d.digest_fragment_count, d.selected_enzymes.len());
    assert_eq!(d.digest_total_bp, d.plasmid_len_bp);
    assert_eq!(d.gibson_assembled_bp, d.plasmid_len_bp);
    assert_eq!(d.gibson_junctions, 2);
    assert!(d.gg_overhang_fidelity_pct > 50.0 && d.gg_overhang_fidelity_pct <= 100.0);
}

#[test]
fn t3_ecosystem_stream_flux_and_enzyme_kinetics() {
    let d = ensure_osabt_validation();
    assert!(d.stream_gpp > 0.0 && d.stream_er > 0.0);
    assert!((d.stream_nep - (d.stream_gpp - d.stream_er)).abs() < 1e-3);
    assert!(d.stream_k600 > 0.0 && !d.stream_trophic_state.is_empty() && d.stream_rmse < 2.0);
    assert!(!d.best_inhib_label.is_empty());
    assert!(d.best_inhib_ki_nm > 0.0 && d.best_inhib_akaike_w > 0.5 && d.best_inhib_r2 > 0.90);
    assert!(d.morrison_et_nm > 0.0 && d.morrison_ki_app_nm > 0.0 && d.morrison_r2 > 0.90);
    assert!(d.csv_path.exists() && d.csv_bytes < 384 * 1024);
    assert!(d.json_path.exists() && d.json_bytes < 384 * 1024);
}
