# Preprints for OSABT (`crvscience-osabt`)

**Author:** Bryan P. White (ORCID: [0000-0003-4191-3511](https://orcid.org/0000-0003-4191-3511))  
**Hub:** https://crvscience.com/apps

This directory contains the self-contained LaTeX source (`.tex`) and single-column PDF (`.pdf`) editions for both:

1. **Companion Application Preprint (`APP-PREPRINT-08`)**
   - **Title:** OSABT: A 20-Module Biology and Biostatistics Workbench for Molecular Labs, Genetics, Kinetics, and Ecology
   - **Files:** `osabt-modular-biostatistics-toolkit.tex`, `osabt-modular-biostatistics-toolkit.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24222388
   - **Figshare Mirror:** https://figshare.com/articles/preprint/OSABT_All-in-One_Toolkit_for_Plant_Physiology_Survival_Curves_Genetics_and_Enzyme_Kinetics/35013763
   - **Web Preprint:** https://crvscience.com/apps/preprints/osabt-modular-biostatistics-toolkit

2. **CRV Science Public Apps Network Overview (`APP-PREPRINT-01`)**
   - **Title:** The CRV Science App Network: Twenty Open-Source Rust Programs for Science, Data, and Simulation
   - **Files:** `crv-science-public-apps-network.tex`, `crv-science-public-apps-network.pdf`
   - **Canonical Zenodo DOI:** https://doi.org/10.5281/zenodo.24222398
   - **Figshare Mirror:** https://figshare.com/articles/preprint/The_CRV_Science_App_Network_Seventeen_Open-Source_Rust_Programs_for_Science_Data_and_Simulation/35013736
   - **Web Preprint:** https://crvscience.com/apps/preprints/crv-science-public-apps-network

## Compiling Standalone LaTeX

Both `.tex` files inline their bibliography so each compiles standalone with standard `pdflatex`:

```bash
pdflatex osabt-modular-biostatistics-toolkit.tex
pdflatex crv-science-public-apps-network.tex
```
