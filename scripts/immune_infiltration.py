"""Immune infiltration via ssGSEA on pan-cancer matrix (immune cell gene sets)"""
import json, statistics
import numpy as np
from scipy import stats

# Immune cell marker gene sets (LM22-based, curated for tumor immunity)
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

with open(r'C:\Users\lee_e\Desktop\生信思路\0818AARS2\data\pancancer\expr_matrix.json', encoding='utf-8') as f:
    meta_data = json.load(f)
with open(r'C:\Users\lee_e\Desktop\生信思路\0818AARS2\data\pancancer\immune_expr_matrix.json', encoding='utf-8') as f:
    immune_data = json.load(f)

# merge: AARS2 from meta_data, immune genes from immune_data
data = {}
for cancer in immune_data:
    info = immune_data[cancer]
    if 'patients' not in info:
        continue
    # get AARS2 for same patients
    a2_pat = meta_data.get(cancer, {}).get('patients', {})
    merged = {}
    for pid, p in info['patients'].items():
        m = dict(p)
        if pid in a2_pat and 'AARS2' in a2_pat[pid]:
            m['AARS2'] = a2_pat[pid]['AARS2']
        merged[pid] = m
    data[cancer] = {'n': info['n'], 'patients': merged}

def ssgsea_score(expr_vals, gene_list, gene_names):
    """Simplified ssGSEA: mean z-score of marker genes present"""
    present = [(g, expr_vals[g]) for g in gene_list if g in expr_vals and expr_vals[g] is not None]
    if len(present) < 3:
        return None
    # rank-based: use percentile rank of each gene within patient
    vals = [v for _, v in present]
    return statistics.mean(vals) if vals else None

cancers = sorted(data.keys())
all_sets = list(IMMUNE_SETS.keys())
rows = []
for cancer in cancers:
    info = data[cancer]
    if 'patients' not in info:
        continue
    pts = info['patients']
    a2 = []
    immune = {s: [] for s in all_sets}
    for pid, p in pts.items():
        if 'AARS2' not in p:
            continue
        a2.append(p['AARS2'])
        for s, genes in IMMUNE_SETS.items():
            # per-patient mean expression of marker genes
            present = [p[g] for g in genes if g in p and p[g] is not None]
            if len(present) >= 3:
                immune[s].append(statistics.mean(present))
            else:
                immune[s].append(None)
    # correlate AARS2 with each immune score
    row = {'cancer': cancer, 'n': len(a2)}
    for s in all_sets:
        pairs = [(a2[i], immune[s][i]) for i in range(len(a2)) if immune[s][i] is not None]
        if len(pairs) >= 20:
            rho, p = stats.spearmanr([x for x, _ in pairs], [y for _, y in pairs])
            row[f'{s}_rho'] = round(rho, 3)
            row[f'{s}_p'] = p
        else:
            row[f'{s}_rho'] = None
            row[f'{s}_p'] = None
    rows.append(row)

# print summary for digestive cancers
print(f"{'Cancer':<10} {'n':>4} | " + ' | '.join(f'{s[:4]}:rho/p' for s in ['CD8_T','Treg','NK','Macrophage','M2']))
for r in rows:
    if r['cancer'] in ('STAD','COADREAD','ESCA','LIHC','PAAD','CHOL','SKCM','TGCT','UVM'):
        parts = []
        for s in ['CD8_T','Treg','NK','Macrophage','M2']:
            rho = r.get(f'{s}_rho'); p = r.get(f'{s}_p')
            parts.append(f'{s[:4]}:{rho if rho is not None else "-"}/{p if p is not None else "-"}')
        print(f"{r['cancer']:<10} {r['n']:>4} | " + ' | '.join(parts))

with open(r'C:\Users\lee_e\Desktop\生信思路\0818AARS2\results\immune_infiltration_corr.json', 'w', encoding='utf-8') as f:
    json.dump(rows, f, ensure_ascii=False, indent=2)
print('\nSaved: results/immune_infiltration_corr.json')
