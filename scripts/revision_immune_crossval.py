"""R3: Immune infiltration cross-validation (Reviewer 2, comment 3)

Method 1: gseapy.ssgsea (proper single-sample GSEA ranking algorithm)
Method 2: NNLS deconvolution (CIBERSORT-style) with a curated immune signature matrix

Compare with the original marker-mean method (immune_infiltration_corr.json):
- If AARS2-immune correlations are consistent across methods, robustness is established.

Outputs:
- results/revision_ssgsea_scores.json (per patient ssGSEA enrichment scores)
- results/revision_ssgsea_corr.json (AARS2 vs ssGSEA scores)
- results/revision_nnls_corr.json (AARS2 vs NNLS cell fractions)
"""
import json, os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
DATA = os.path.join(BASE, 'data', 'pancancer')

def load_json(p):
    with open(p, encoding='utf-8') as f:
        return json.load(f)

# Immune cell marker gene sets (same as original analysis)
IMMUNE_SETS = {
    'CD8_T': ['CD8A','CD8B','GZMA','GZMB','PRF1','IFNG','TBX21','CXCR3'],
    'CD4_T': ['CD4','IL7R','CCR7','LEF1','TCF7','FOXP3'],
    'Treg': ['FOXP3','IL2RA','CTLA4','IKZF2','CCR8','TNFRSF18'],
    'NK': ['NKG7','KLRD1','KLRK1','NCR1','GNLY','KIR2DL3'],
    'B': ['CD19','MS4A1','CD79A','CD79B','BANK1'],
    'Macrophage': ['CD68','CD163','CSF1R','MRC1','ITGAM','CD14'],
    'M1': ['NOS2','IL1B','TNF','CXCL9','CXCL10','IL12A'],
    'M2': ['CD163','MRC1','MSR1','IL10','TGFB1','ARG1'],
    'DC': ['ITGAX','CD1C','CLEC9A','BATF3','FLT3'],
    'Neutrophil': ['FCGR3B','CSF3R','S100A8','S100A9','CEACAM8'],
    'Monocyte': ['CD14','FCGR3A','LYZ','S100A12'],
}

GI_CANCERS = ['STAD', 'ESCA', 'LIHC', 'PAAD', 'COADREAD']

# ---------------------------------------------------------------
# Load and prepare expression matrices
# ---------------------------------------------------------------
expr = load_json(os.path.join(DATA, 'expr_matrix.json'))
immune_data = load_json(os.path.join(DATA, 'immune_expr_matrix.json'))

# Build per-cancer DataFrame: rows=patients, cols=genes (AARS2 + immune genes)
def build_df(cancer):
    """Return (df, aars2_series) with patients as rows."""
    pats = immune_data[cancer].get('patients', {})
    a2_pat = expr.get(cancer, {}).get('patients', {})
    genes = set()
    for p in pats.values():
        genes.update(p.keys())
    for p in a2_pat.values():
        genes.update(p.keys())
    genes = sorted(genes)
    rows = []
    pids = []
    aars2 = []
    for pid, p in pats.items():
        row = {g: p.get(g) for g in genes}
        a2 = a2_pat.get(pid, {}).get('AARS2')
        if a2 is None:
            continue
        rows.append(row)
        pids.append(pid)
        aars2.append(a2)
    df = pd.DataFrame(rows, index=pids)
    df = df.apply(pd.to_numeric, errors='coerce')
    return df, pd.Series(aars2, index=pids)

# ---------------------------------------------------------------
# Method 1: gseapy.ssgsea per cancer
# ---------------------------------------------------------------
import gseapy as gp

ssgsea_results = {}
ssgsea_corr = {}
for cancer in GI_CANCERS:
    try:
        df, a2 = build_df(cancer)
        if df.shape[0] < 30:
            continue
        # Drop genes with >50% missing
        keep = df.notna().sum(axis=0) > df.shape[0] * 0.5
        df = df.loc[:, keep]
        # Fill remaining NA with gene mean
        df = df.fillna(df.mean())
        # ssGSEA needs genes as rows, samples as columns
        expr_t = df.T
        # Run ssGSEA
        ss = gp.ssgsea(expr_t, gene_sets=IMMUNE_SETS, outdir=None, min_size=3, max_size=500,
                       sample_norm_method='rank', no_plot=True, threads=2, verbose=False)
        ss_df = ss.res2d.pivot(index='Term', columns='Name', values='ES')
        # Correlate AARS2 with each cell type score
        row = {'cancer': cancer, 'n': len(a2)}
        for cell in ss_df.index:
            scores = ss_df.loc[cell].astype(float)
            common = a2.index.intersection(scores.index)
            if len(common) < 30:
                continue
            rho, p = stats.spearmanr(a2[common], scores[common])
            row[f'{cell}_rho'] = float(rho)
            row[f'{cell}_p'] = float(p)
        ssgsea_corr[cancer] = row
        ssgsea_results[cancer] = ss_df.to_dict()
        print(f'{cancer}: ssGSEA done ({len(a2)} patients)')
    except Exception as e:
        print(f'{cancer}: ssGSEA FAIL {e}')

with open(os.path.join(BASE, 'results', 'revision_ssgsea_corr.json'), 'w', encoding='utf-8') as f:
    json.dump(ssgsea_corr, f, ensure_ascii=False, indent=2)
print('Saved: revision_ssgsea_corr.json')

# FDR for ssGSEA correlations
for cancer in ssgsea_corr:
    tests = [(k, v) for k, v in ssgsea_corr[cancer].items() if k.endswith('_p')]
    if tests:
        _, qs, _, _ = multipletests([v for _, v in tests], method='fdr_bh')
        for (k, _), q in zip(tests, qs):
            ssgsea_corr[cancer][k.replace('_p', '_q')] = float(q)
with open(os.path.join(BASE, 'results', 'revision_ssgsea_corr.json'), 'w', encoding='utf-8') as f:
    json.dump(ssgsea_corr, f, ensure_ascii=False, indent=2)
print('FDR applied to ssGSEA corr')

# ---------------------------------------------------------------
# Method 2: NNLS deconvolution (CIBERSORT-style)
# ---------------------------------------------------------------
from scipy.optimize import nnls

# Build signature matrix: immune cell types x marker genes (mean expression per gene across cancers)
# Use the curated marker sets; genes are shared across cell types
def build_signature():
    """Signature matrix: rows=cell types, cols=genes (mean marker expression across all GI patients)."""
    cell_genes = {}
    for cell, genes in IMMUNE_SETS.items():
        cell_genes[cell] = genes
    # Compute mean expression per gene per cell (average across all GI cancer patients)
    gene_means = {}
    for cancer in GI_CANCERS:
        pats = immune_data.get(cancer, {}).get('patients', {})
        for p in pats.values():
            for g, v in p.items():
                if v is not None:
                    gene_means.setdefault(g, []).append(v)
    sig_rows = {}
    for cell, genes in cell_genes.items():
        vals = [np.mean(gene_means[g]) for g in genes if g in gene_means]
        sig_rows[cell] = {g: np.mean(gene_means[g]) for g in genes if g in gene_means}
    return sig_rows

sig = build_signature()
all_sig_genes = sorted(set(g for genes in sig.values() for g in genes))
sig_matrix = pd.DataFrame(0.0, index=list(sig.keys()), columns=all_sig_genes)
for cell, gdict in sig.items():
    for g, v in gdict.items():
        sig_matrix.loc[cell, g] = v

nnls_corr = {}
for cancer in GI_CANCERS:
    try:
        df, a2 = build_df(cancer)
        if df.shape[0] < 30:
            continue
        # Genes in signature present in data
        sig_genes = [g for g in all_sig_genes if g in df.columns]
        if len(sig_genes) < 10:
            continue
        X = sig_matrix[sig_genes].values  # (cell x gene)
        # For each patient, solve NNLS: min ||X^T f - expr||  (f = cell fractions)
        fractions = []
        for _, row in df.iterrows():
            y = row[sig_genes].fillna(0).values
            if np.sum(y) == 0:
                fractions.append(np.zeros(X.shape[0]))
                continue
            f, _ = nnls(X.T, y)
            fractions.append(f)
        frac_df = pd.DataFrame(fractions, columns=sig_matrix.index, index=df.index)
        # Normalize to sum 1
        frac_df = frac_df.div(frac_df.sum(axis=1), axis=0).replace([np.inf, -np.inf], 0).fillna(0)
        # Correlate AARS2 with each cell fraction
        row = {'cancer': cancer, 'n': len(a2)}
        for cell in frac_df.columns:
            rho, p = stats.spearmanr(a2, frac_df[cell])
            row[f'{cell}_rho'] = float(rho)
            row[f'{cell}_p'] = float(p)
        nnls_corr[cancer] = row
        print(f'{cancer}: NNLS done ({len(a2)} patients)')
    except Exception as e:
        print(f'{cancer}: NNLS FAIL {e}')

with open(os.path.join(BASE, 'results', 'revision_nnls_corr.json'), 'w', encoding='utf-8') as f:
    json.dump(nnls_corr, f, ensure_ascii=False, indent=2)

# FDR for NNLS correlations
for cancer in nnls_corr:
    tests = [(k, v) for k, v in nnls_corr[cancer].items() if k.endswith('_p')]
    if tests:
        _, qs, _, _ = multipletests([v for _, v in tests], method='fdr_bh')
        for (k, _), q in zip(tests, qs):
            nnls_corr[cancer][k.replace('_p', '_q')] = float(q)
with open(os.path.join(BASE, 'results', 'revision_nnls_corr.json'), 'w', encoding='utf-8') as f:
    json.dump(nnls_corr, f, ensure_ascii=False, indent=2)
print('Saved: revision_nnls_corr.json (with FDR)')

# ---------------------------------------------------------------
# Comparison table
# ---------------------------------------------------------------
print('\n' + '='*110)
print('AARS2 vs CD8_T correlation: 3 methods comparison')
print('='*110)
orig = load_json(os.path.join(BASE, 'results', 'immune_infiltration_corr.json'))
orig_by_cancer = {r['cancer']: r for r in orig}
print(f'{"Cancer":8} {"marker-mean rho":15} {"ssGSEA rho":12} {"NNLS rho":10}')
for cancer in GI_CANCERS:
    o = orig_by_cancer.get(cancer, {}).get('CD8_T_rho')
    s = ssgsea_corr.get(cancer, {}).get('CD8_T_rho')
    n = nnls_corr.get(cancer, {}).get('CD8_T_rho')
    print(f'{cancer:8} {o if o is not None else "NA":>15} {s if s is not None else "NA":>12} {n if n is not None else "NA":>10}')

print('\nAARS2 vs M2 correlation: 3 methods comparison')
print(f'{"Cancer":8} {"marker-mean rho":15} {"ssGSEA rho":12} {"NNLS rho":10}')
for cancer in GI_CANCERS:
    o = orig_by_cancer.get(cancer, {}).get('M2_rho')
    s = ssgsea_corr.get(cancer, {}).get('M2_rho')
    n = nnls_corr.get(cancer, {}).get('M2_rho')
    print(f'{cancer:8} {o if o is not None else "NA":>15} {s if s is not None else "NA":>12} {n if n is not None else "NA":>10}')
