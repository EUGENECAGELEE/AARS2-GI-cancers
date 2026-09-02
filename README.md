# AARS2 in Gastrointestinal Cancers — Analysis Code and Data

This repository contains the analysis scripts and processed data underlying the manuscript:

> Wang X, Li H. *Mitochondrial Alanyl-tRNA Synthetase AARS2 Links Heat Shock Response to Immune Microenvironment Remodeling in Gastrointestinal Cancers: A Pan-Cancer Multi-Omics Analysis.* Submitted to PLOS ONE.

## Overview

All analyses are purely bioinformatics and use publicly available datasets:
- **TCGA** pan-cancer bulk transcriptomic data accessed via the **cBioPortal API** (https://www.cbioportal.org)
- **GEO GSE149614** single-cell RNA-seq of hepatocellular carcinoma (71,913 cells)

No private or restricted data are used.

## Repository structure

```
scripts/   # Python scripts corresponding to the Methods section
data/      # Processed result tables (.xlsx) and summary statistics (.json)
```

### Scripts (Methods mapping)

| Script | Method section | Description |
|---|---|---|
| `immune_infiltration.py` | 2.3 | Immune infiltration scores (IMMUNE_SETS: CD8_T, Treg, M2, etc.) |
| `revision_core.py` | 2.2, 2.5 | Pathway correlations + partial correlation (controlling for glycolysis / HIF1A) |
| `revision_cox.py` | 2.5 | Continuous Cox regression for OS |
| `revision_rcs.py` | 2.5 | Restricted cubic spline (RCS) models |
| `revision_control_immune.py` | 2.3 | Control-gene specificity analysis (24 mitochondrial genes) |
| `revision_hsf1_analysis.py` | 2.3 | AARS2 vs HSF1 / CD274 (PD-L1) / IDO1 |
| `revision_immune_crossval.py` | 2.3 | ssGSEA + marker-based NNLS cross-validation |
| `revision_sc_malignant.py` | 2.4 | Single-cell malignant-hepatocyte validation |
| `revision_sc_umi_threshold.py` | 2.4 | Single-cell UMI threshold sensitivity |
| `revision_partial_fdr.py` | 2.5 | FDR correction for partial correlations |
| `revision_fetch_control.py` | 2.2 | Fetch control-gene expression via cBioPortal API |
| `revision_fetch_hsf1.py` | 2.3 | Fetch HSF1 / CD274 / IDO1 expression via cBioPortal API |

### Data

- `data/*.xlsx` — Supplementary Tables S1–S5 (sample sizes, pathway correlations, immune infiltration, control-gene analysis, single-cell statistics)
- `data/*.json` — correlation, survival, and immune summary statistics used to generate the figures

## Requirements

- Python 3.13
- See individual scripts for imports: `numpy`, `pandas`, `scipy`, `statsmodels`, `lifelines`, `sklearn`, `gseapy`, `scanpy`, `requests`

## Reproducibility note

Raw TCGA/GEO data are re-downloaded at runtime from public APIs (cBioPortal, GEO).
This repository provides the processed outputs and the exact analysis code that
generated the manuscript results, satisfying the PLOS ONE data-availability policy.

## License

CC-BY 4.0 — code; data derived from TCGA (cBioPortal) and GEO GSE149614 under their respective access policies.
