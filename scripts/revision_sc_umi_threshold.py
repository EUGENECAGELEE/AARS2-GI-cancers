"""P1-5: Correct Fisher OR computation (site fraction vs Normal fraction).
Primary: expr > 0 (>=1 UMI). Robustness: top-25% of positive hepatocytes.
"""
import sys, os, json
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import scanpy as sc
from scipy import stats

base = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
adata = sc.read_h5ad(os.path.join(base, 'data', 'scRNA', 'GSE149614_embedded.h5ad'))
idx = adata.var_names.get_loc('AARS2')
X = adata.X[:, idx]
if hasattr(X, 'toarray'):
    X = X.toarray().ravel()
else:
    X = np.asarray(X).ravel()

obs = adata.obs
hep = (obs['celltype'] == 'Hepatocyte').values

def frac_by_site(pos_mask):
    out = {}
    for site in ['Normal', 'Tumor', 'PVTT']:
        m = hep & (obs['site'].values == site)
        n = m.sum()
        k = (pos_mask & m).sum()
        out[site] = (int(n), int(k), float(k / n) if n else 0.0)
    return out

def fisher_vs_normal(site, frac):
    n_s, k_s, _ = frac[site]
    n_n, k_n, _ = frac['Normal']
    table = [[k_s, n_s - k_s], [k_n, n_n - k_n]]
    orr, p = stats.fisher_exact(np.array(table), alternative='greater')
    return float(orr), float(p), table

results = {}

# ---- Primary: expr > 0 ----
pos = X > 0
frac = frac_by_site(pos)
results['primary_gt0'] = {'frac': frac}
for site in ['Tumor', 'PVTT']:
    orr, p, tbl = fisher_vs_normal(site, frac)
    results['primary_gt0'][f'{site}_vs_Normal'] = {'OR': orr, 'p': p, 'table': tbl}
    print(f'Primary >0: {site} OR={orr:.2f} p={p:.2e}')

# ---- Robustness: top-25% of positive hepatocytes ----
hep_x = X[hep]
pos_vals = hep_x[hep_x > 0]
thr = np.percentile(pos_vals, 75) if len(pos_vals) else 0
pos_top = hep & (X >= thr)
frac2 = frac_by_site(pos_top)
results['top25'] = {'threshold': float(thr), 'frac': frac2}
for site in ['Tumor', 'PVTT']:
    orr, p, tbl = fisher_vs_normal(site, frac2)
    results['top25'][f'{site}_vs_Normal'] = {'OR': orr, 'p': p, 'table': tbl}
    print(f'Top25: {site} OR={orr:.2f} p={p:.2e}')

# ---- Robustness: expr >= 0.42 (threshold of top25) applied globally ----
pos3 = X >= thr
frac3 = frac_by_site(pos3)
results['global_thr'] = {'threshold': float(thr), 'frac': frac3}
for site in ['Tumor', 'PVTT']:
    orr, p, tbl = fisher_vs_normal(site, frac3)
    results['global_thr'][f'{site}_vs_Normal'] = {'OR': orr, 'p': p, 'table': tbl}
    print(f'Global>=thr: {site} OR={orr:.2f} p={p:.2e}')

with open(base + r'\results\revision_sc_threshold.json', 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=1, default=str)
print('\nSaved: results/revision_sc_threshold.json')
