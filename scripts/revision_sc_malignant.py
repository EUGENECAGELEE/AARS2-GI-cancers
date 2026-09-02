"""R4-R7: Malignant vs non-malignant hepatocyte analysis.
GSE149614 has no direct malignant annotation; use sample site (Tumor/PVTT = malignant context,
Normal = non-malignant) + fetal/malignant markers. Approach:
1. Hepatocytes split by site (Normal vs Tumor+PVTT).
2. Correlate AARS2 expression with malignant markers (GPC3, AFP, EPCAM, KRT19, SPP1, CD44)
   within tumor-context hepatocytes, and compare AARS2-high fraction in
   marker-high vs marker-low hepatocytes.
"""
import sys, os, json
import numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
p = os.path.join(base, 'data', 'scRNA', 'GSE149614_embedded.h5ad')
import anndata as ad

MAL_MARKERS = ['GPC3', 'AFP', 'EPCAM', 'KRT19', 'SPP1', 'CD44']
HEP_MARKERS = ['ALB', 'APOA2', 'TTR', 'SLC10A1', 'ASGR1']  # differentiation markers

adata = ad.read_h5ad(p, backed='r')
obs = adata.obs
hep_mask = (obs['celltype'] == 'Hepatocyte').values
print('hepatocytes:', hep_mask.sum())

# site split
sites = obs['site'].values
hep_sites = sites[hep_mask]
print('site counts in hepatocytes:', {s: int((hep_sites == s).sum()) for s in ['Normal', 'Tumor', 'PVTT', 'Lymph']})

# get AARS2 + markers for hepatocytes
gene_idx = {g: list(adata.var_names).index(g) for g in ['AARS2'] + MAL_MARKERS + HEP_MARKERS}
n_hep = int(hep_mask.sum())
mat = np.zeros((n_hep, len(gene_idx)), dtype=np.float32)
rows = np.where(hep_mask)[0]
cols = list(gene_idx.values())
# backed h5py: read row-wise (single vector index allowed), select cols manually
X = adata.X
for bi in range(0, n_hep, 3000):
    sl = rows[bi:bi+3000]
    sub = X[sl, :]  # full row slice (vector index on rows only)
    if hasattr(sub, 'toarray'):
        sub = sub.toarray()
    sub = np.asarray(sub, dtype=np.float32)
    mat[bi:bi+3000] = sub[:, cols]
adata.file.close()

names = list(gene_idx.keys())
idx = {g: i for i, g in enumerate(names)}
aars2 = mat[:, idx['AARS2']]
print('AARS2 in hepatocytes: >0 fraction = %.4f, mean = %.5f' % ((aars2 > 0).mean(), aars2.mean()))

# 1. detection rate & mean by site
res = {}
for site in ['Normal', 'Tumor', 'PVTT']:
    m = hep_sites == site
    res[site] = {
        'n': int(m.sum()),
        'detection_rate': float((aars2[m] > 0).mean()),
        'mean_expr': float(aars2[m].mean()),
        'median_expr': float(np.median(aars2[m])),
    }
print(json.dumps(res, indent=1))

# 2. Within tumor-context hepatocytes: AARS2 vs malignant markers (Spearman)
tumor_ctx = (hep_sites == 'Tumor') | (hep_sites == 'PVTT')
print('\nAARS2 vs markers in tumor-context hepatocytes (n=%d):' % tumor_ctx.sum())
corr = {}
for g in MAL_MARKERS + HEP_MARKERS:
    gi = idx[g]
    x = aars2[tumor_ctx]
    y = mat[tumor_ctx, gi]
    # use cells with any expression of the marker to avoid dropout-dominated rho
    both = (x > 0) & (y > 0)
    if both.sum() > 30:
        r, pv = stats.spearmanr(x[both], y[both])
        corr[g] = {'rho': float(r), 'p': float(pv), 'n_both': int(both.sum())}
    else:
        r_all, pv_all = stats.spearmanr(x, y)
        corr[g] = {'rho': float(r_all), 'p': float(pv_all), 'n_both': 0, 'note': 'few co-expressing cells, used all'}
print(json.dumps(corr, indent=1))

# 3. AARS2-high fraction in GPC3/AFP-high vs low hepatocytes (tumor context)
out3 = {}
for g in ['GPC3', 'AFP']:
    gi = idx[g]
    y = mat[tumor_ctx, gi]
    hi = y > np.percentile(y[y > 0], 75) if (y > 0).sum() > 10 else y > 0
    lo = ~hi
    a_hi = aars2[tumor_ctx]
    fr_hi = (a_hi[hi] > 0).mean()
    fr_lo = (a_hi[lo] > 0).mean()
    tbl = [[int((a_hi[hi] > 0).sum()), int((a_hi[hi] == 0).sum())],
           [int((a_hi[lo] > 0).sum()), int((a_hi[lo] == 0).sum())]]
    or_, pv = stats.fisher_exact(tbl)
    out3[g] = {'marker_high_AARS2pos_frac': float(fr_hi), 'marker_low_AARS2pos_frac': float(fr_lo),
               'fisher_OR': float(or_), 'p': float(pv)}
print('\nAARS2+ fraction by malignant marker status (tumor context):')
print(json.dumps(out3, indent=1))

out_all = {'detection_by_site': res, 'marker_corr_tumor_ctx': corr, 'malignant_strat': out3}
with open(os.path.join(base, 'results', 'revision_sc_malignant.json'), 'w', encoding='utf-8') as f:
    json.dump(out_all, f, ensure_ascii=False, indent=1)
print('\nsaved revision_sc_malignant.json')
