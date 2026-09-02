"""P1-3 (R6-Q6): FDR-correct partial correlation p-values and update paper text.

All partial correlation p-values (hsp|glyc, hsp|hif1a, lac_sub|glyc, lac_sub|hif1a)
across the 5 GI cancers are pooled for BH-FDR (same convention as revision_core.py).
"""
import json, sys, re
sys.stdout.reconfigure(encoding='utf-8')
from statsmodels.stats.multitest import multipletests

BASE = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
d = json.load(open(BASE + r'\results\revision_pathway_corr_fdr.json', encoding='utf-8'))

GI = ['STAD', 'ESCA', 'LIHC', 'PAAD', 'COADREAD']
METRICS = ['hsp|glyc', 'hsp|hif1a', 'lac_sub|glyc', 'lac_sub|hif1a']

# Collect all p values for the four partial metrics across GI cancers
for metric in METRICS:
    ps = [(c, d[c][f'hsp_hsp|glyc_p'] if metric == 'hsp|glyc' else
              d[c][f'hsp_hsp|hif1a_p'] if metric == 'hsp|hif1a' else
              d[c][f'lac_sub_lac_sub|glyc_p'] if metric == 'lac_sub|glyc' else
              d[c][f'lac_sub_lac_sub|hif1a_p']) for c in GI if c in d]
    pvals = [p for _, p in ps]
    if not pvals:
        continue
    _, qs, _, _ = multipletests(pvals, method='fdr_bh')
    print(f'--- {metric} (n={len(pvals)}) ---')
    for (c, p), q in zip(ps, qs):
        rho_key = f'hsp_hsp|glyc_rho' if metric == 'hsp|glyc' else \
                  f'hsp_hsp|hif1a_rho' if metric == 'hsp|hif1a' else \
                  f'lac_sub_lac_sub|glyc_rho' if metric == 'lac_sub|glyc' else \
                  f'lac_sub_lac_sub|hif1a_rho'
        print(f'  {c}: rho={d[c][rho_key]:.3f} p={p:.2e} q={q:.2e}')
