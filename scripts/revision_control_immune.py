"""R4: Control gene specificity (Reviewer 2 comment 4)

Compare AARS2 vs control genes (mitochondrial markers, mitochondrial aaRS,
TCA enzymes, mito TF/HSP) for immune infiltration correlations, using the
SAME immune scoring pipeline as immune_infiltration.py.

Question: is AARS2's negative correlation with CD8_T/NK/M2/Treg specific,
or shared by all mitochondrial genes?
"""
import json, os, statistics
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
DATA = os.path.join(BASE, 'data', 'pancancer')

def load_json(p):
    with open(p, encoding='utf-8') as f:
        return json.load(f)

# Immune cell marker gene sets (identical to immune_infiltration.py)
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

# Control gene categories
GENE_CAT = {}
for g in ['COX4I1','COX5A','TOMM20','TOMM22']: GENE_CAT[g] = 'Mito marker'
for g in ['DARS2','EARS2','LARS2','MARS2','TARS2','YARS2','SARS2','GARS1']: GENE_CAT[g] = 'Mito aaRS'
for g in ['PDHA1','SDHB','MDH2','CS','IDH2','SUCLA2','OGDH','DLST','FH']: GENE_CAT[g] = 'TCA enzyme'
for g in ['TFAM']: GENE_CAT[g] = 'Mito TF'
for g in ['HSPA9','HSPD1']: GENE_CAT[g] = 'Mito HSP'

ALL_CONTROLS = list(GENE_CAT.keys())
ALL_GENES = ['AARS2'] + ALL_CONTROLS

GI_CANCERS = ['STAD', 'ESCA', 'LIHC', 'PAAD', 'COADREAD']
IMMUNE_TYPES = ['CD8_T', 'CD4_T', 'Treg', 'NK', 'B', 'Macrophage', 'M1', 'M2', 'DC', 'Neutrophil', 'Monocyte']

# ---------------------------------------------------------------
# Load data
# ---------------------------------------------------------------
expr = load_json(os.path.join(DATA, 'expr_matrix.json'))
immune_data = load_json(os.path.join(DATA, 'immune_expr_matrix.json'))
control = load_json(os.path.join(DATA, 'control_expr_matrix.json'))

# Build merged patient dicts: gene expr + immune gene expr
merged = {}
for cancer in immune_data:
    if 'patients' not in immune_data[cancer]:
        continue
    a2_pat = expr.get(cancer, {}).get('patients', {})
    ctrl_pat = control.get(cancer, {}).get('patients', {})
    m = {}
    for pid, p in immune_data[cancer]['patients'].items():
        d = dict(p)
        if pid in a2_pat and 'AARS2' in a2_pat[pid]:
            d['AARS2'] = a2_pat[pid]['AARS2']
        if pid in ctrl_pat:
            for g in ALL_CONTROLS:
                if g in ctrl_pat[pid]:
                    d[g] = ctrl_pat[pid][g]
        m[pid] = d
    merged[cancer] = {'n': immune_data[cancer]['n'], 'patients': m}

# ---------------------------------------------------------------
# Compute immune scores per patient (same method as original)
# ---------------------------------------------------------------
def immune_scores(patient_dict):
    """Compute mean marker expression per immune cell type."""
    scores = {}
    for s, genes in IMMUNE_SETS.items():
        present = [patient_dict[g] for g in genes
                   if g in patient_dict and patient_dict[g] is not None]
        scores[s] = statistics.mean(present) if len(present) >= 3 else None
    return scores

# ---------------------------------------------------------------
# Correlations: each gene x each immune type x each GI cancer
# ---------------------------------------------------------------
results = {}
for cancer in GI_CANCERS:
    if cancer not in merged:
        continue
    pts = merged[cancer]['patients']
    # Precompute immune scores per patient
    imm_scores = {pid: immune_scores(p) for pid, p in pts.items()}
    cancer_res = {}
    for gene in ALL_GENES:
        gene_res = {}
        for s in IMMUNE_TYPES:
            pairs = []
            for pid, p in pts.items():
                if gene not in p or p[gene] is None:
                    continue
                isc = imm_scores[pid][s]
                if isc is None:
                    continue
                pairs.append((float(p[gene]), float(isc)))
            if len(pairs) >= 20:
                rho, pval = stats.spearmanr([x for x, _ in pairs], [y for _, y in pairs])
                gene_res[s] = {'rho': round(float(rho), 4), 'p': float(pval), 'n': len(pairs)}
            else:
                gene_res[s] = None
        cancer_res[gene] = gene_res
    results[cancer] = cancer_res
    print(f'{cancer}: done ({len(pts)} patients)')

# ---------------------------------------------------------------
# FDR across all tests per cancer
# ---------------------------------------------------------------
for cancer in results:
    tests = []
    for gene in ALL_GENES:
        for s in IMMUNE_TYPES:
            r = results[cancer][gene][s]
            if r is not None:
                tests.append((gene, s, r))
    if tests:
        _, qs, _, _ = multipletests([r['p'] for _, _, r in tests], method='fdr_bh')
        for (gene, s, r), q in zip(tests, qs):
            results[cancer][gene][s]['q'] = float(q)

# Save full
with open(os.path.join(BASE, 'results', 'revision_control_immune.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print('Saved: revision_control_immune.json')

# ---------------------------------------------------------------
# Summary table: AARS2 vs controls for CD8_T / NK / M2 / Treg
# ---------------------------------------------------------------
print('\n' + '='*110)
print(f'{"Cancer":7} {"Gene":8} {"Cat":12} | ' + ' | '.join(f'{s}' for s in ['CD8_T', 'NK', 'M2', 'Treg']))
print('='*110)
for cancer in GI_CANCERS:
    if cancer not in results:
        continue
    for gene in ALL_GENES:
        r = results[cancer][gene]
        parts = []
        for s in ['CD8_T', 'NK', 'M2', 'Treg']:
            v = r[s]
            if v is None:
                parts.append(f'{"-":>10}')
            else:
                sig = '*' if v['q'] < 0.05 else ''
                parts.append(f'{v["rho"]:+.3f}{sig:1}(q={v["q"]:.2e})')
        cat = GENE_CAT.get(gene, 'AARS2')
        print(f'{cancer:7} {gene:8} {cat:12} | ' + ' | '.join(parts))
    print('-'*110)

# ---------------------------------------------------------------
# Key statistic: mean |rho| of AARS2 vs controls for CD8_T
# ---------------------------------------------------------------
print('\n=== Effect size comparison (|rho| for CD8_T, 5 GI cancers) ===')
for s in ['CD8_T', 'NK', 'M2', 'Treg']:
    a2_rhos = []
    ctrl_rhos = []
    for cancer in GI_CANCERS:
        if cancer not in results:
            continue
        for gene in ALL_GENES:
            v = results[cancer][gene][s]
            if v is None:
                continue
            if gene == 'AARS2':
                a2_rhos.append(abs(v['rho']))
            else:
                ctrl_rhos.append(abs(v['rho']))
    if a2_rhos and ctrl_rhos:
        print(f'{s}: AARS2 mean|rho|={np.mean(a2_rhos):.3f} (n={len(a2_rhos)}), '
              f'controls mean|rho|={np.mean(ctrl_rhos):.3f} (n={len(ctrl_rhos)})')

# Mean AARS2 vs mean control rho per cancer for CD8_T
print('\n=== CD8_T: AARS2 rho vs control range per cancer ===')
for cancer in GI_CANCERS:
    if cancer not in results:
        continue
    a2 = results[cancer]['AARS2']['CD8_T']
    ctrl_vals = []
    for g in ALL_CONTROLS:
        v = results[cancer][g]['CD8_T']
        if v is not None:
            ctrl_vals.append(v['rho'])
    if a2 is not None and ctrl_vals:
        print(f'{cancer}: AARS2 rho={a2["rho"]:+.3f} | controls min={min(ctrl_vals):+.3f} '
              f'max={max(ctrl_vals):+.3f} mean={np.mean(ctrl_vals):+.3f} '
              f'[{", ".join(f"{g}:{v:+.2f}" for g, v in zip(ALL_CONTROLS, ctrl_vals))}]')
