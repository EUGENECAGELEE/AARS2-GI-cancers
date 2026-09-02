"""R4: AARS2 vs HSF1, CD274 (PD-L1), IDO1 correlations in 5 GI cancers.
Converts extra patientIds to uniquePatientKey (base64) to match expr matrix."""
import sys, json, os, base64
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf-8')
base = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'

expr = json.load(open(os.path.join(base, 'data', 'pancancer', 'expr_matrix.json'), encoding='utf-8'))
extra = json.load(open(os.path.join(base, 'data', 'pancancer', 'hsf1_pdl1_ido1_expr.json'), encoding='utf-8'))

CANCERS = ['STAD', 'ESCA', 'LIHC', 'PAAD', 'COADREAD']
out = {}

for cancer in CANCERS:
    info = extra.get(cancer)
    if not info or not info.get('patients'):
        print(f'{cancer}: no extra data')
        continue
    pats = expr[cancer]['patients']  # uniquePatientKey keys
    ex = info['patients']            # plain patientId keys
    # build plain->key map from expr
    key_to_plain = {}
    for k in pats.keys():
        try:
            decoded = base64.b64decode(k + '==').decode('utf-8')
            key_to_plain[decoded.split(':')[0]] = k
        except Exception:
            continue
    aars2_vals, hsf1_vals, cd274_vals, ido1_vals = [], [], [], []
    for pid, gd in ex.items():
        if pid in key_to_plain:
            key = key_to_plain[pid]
        else:
            match = next((k for p, k in key_to_plain.items() if p.startswith(pid) or pid.startswith(p)), None)
            if match is None:
                continue
            key = match
        a = pats[key].get('AARS2')
        h = gd.get('HSF1')
        c = gd.get('CD274')
        i = gd.get('IDO1')
        if a is None or h is None:
            continue
        aars2_vals.append(a); hsf1_vals.append(h)
        cd274_vals.append(c if c is not None else np.nan)
        ido1_vals.append(i if i is not None else np.nan)
    n = len(aars2_vals)
    print(f'== {cancer}: matched {n} patients')
    if n < 30:
        continue
    a = np.array(aars2_vals)
    res = {}
    for nm, y in [('HSF1', np.array(hsf1_vals)), ('CD274', np.array(cd274_vals, dtype=float)),
                  ('IDO1', np.array(ido1_vals, dtype=float))]:
        m = ~np.isnan(y)
        if m.sum() < 30:
            res[nm] = {'n': int(m.sum()), 'note': 'insufficient'}
            continue
        r, p = stats.spearmanr(a[m], y[m])
        res[nm] = {'n': int(m.sum()), 'rho': float(r), 'p': float(p)}
        print(f'  AARS2 vs {nm}: rho={r:.4f} p={p:.2e} n={m.sum()}')
    out[cancer] = res

# FDR across 5 cancers x 3 genes
all_p = []
for c in CANCERS:
    if c in out:
        for g in ['HSF1', 'CD274', 'IDO1']:
            if g in out[c] and 'p' in out[c][g]:
                all_p.append((c, g, out[c][g]['p']))
if all_p:
    _, qvals, _, _ = multipletests([x[2] for x in all_p], method='fdr_bh')
    for (c, g, p), q in zip(all_p, qvals):
        out[c][g]['q'] = float(q)
        print(f'  FDR {c} {g}: q={q:.4f}')

with open(os.path.join(base, 'results', 'revision_hsf1_pdl1_ido1.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
print('\nsaved revision_hsf1_pdl1_ido1.json')
