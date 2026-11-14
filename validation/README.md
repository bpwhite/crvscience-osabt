# Headless Validation Artifacts (`crvscience-osabt`)

Generated directly by running the compiled Rust integration test suite against real implemented `osabt` (`v2.0.0`) modules:

```bash
cargo test --test headless_validation
```

All outputs stay under `384 KB` (`393,216` bytes) per file:

1. **`osabt_benchmark_suite.csv`** (`T1`–`T3`): Quantitative benchmark metrics emitted across `flux_leaf::fit_fvcb_aci_curve`, `thermo::calc_perfect_duplex`, `survival_core::{compute_kaplan_meier_curves, compute_log_rank_test, fit_cox_proportional_hazards}`, `survival_tox::fit_ecotox_lc50`, `popgen_core::{analyze_popgen_fasta, analyze_diploid_vcf}`, `plasmid_digest::{find_unique_cutters, simulate_restriction_digest, score_overhang_set_fidelity}`, `gibson::design_gibson_assembly`, `flux_stream::fit_odum_stream_metabolism`, and `kinetics_core::{fit_global_reversible_inhibition, fit_morrison_tight_binding}`.
2. **`osabt_validation_summary.json`**: Structured validation summary containing empirical parameter estimates, statistical test outputs, fragment assembly conservation checks, and artifact byte sizes.
