"""Revision core analyses (R5/R6/R11/R12):
R5:  lactylation score WITHOUT AARS2/AARS1 (fix circular logic)
R6:  partial Spearman correlation (control glycolysis / HIF1A)
R11: FDR q-values for all correlation tables
R12: Cox continuous + restricted cubic spline survival
R4:  control genes vs immune (待control数据到位)
"""
import json, os, math
import numpy as np
from scipy import stats
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

BASE = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
DATA = os.path.join(BASE, 'data', 'pancancer')

def load_json(p):
    with open(p, encoding='utf-8') as f:
        return json.load(f)

expr = load_json(os.path.join(DATA, 'expr_matrix.json'))

LAC_SUBSTRATE  = ['TP53', 'YAP1', 'CGAS', 'ATF6', 'TDO2', 'TEAD1']
LAC_FLUX       = ['LDHA', 'LDHB', 'SLC16A1', 'SLC16A3', 'HK1', 'HK2', 'PFKL',
                   'PFKM', 'PFKP', 'ALDOA', 'TPI1', 'GAPDH', 'PGK1', 'PGAM1', 'ENO1', 'PKM']
HSP_GENES      = ['HSPA1A', 'HSPA1B', 'HSP90AA1', 'HSP90AB1', 'HSPB1', 'HSPD1',
                   'HSPA9', 'HSPA4', 'HSPA5', 'HSPA8', 'HSPH1', 'HSPE1',
                   'DNAJB1', 'DNAJA1', 'DNAJC7']
HIF1A_TARGETS  = ['HIF1A', 'VEGFA', 'SLC2A1', 'SLC2A3', 'LDHA', 'PGK1', 'ENO1',
                   'PKM', 'PDK1', 'HK2', 'PFKL', 'TPI1', 'GAPDH', 'ALDOA',
                   'PGAM1', 'SLC16A1', 'SLC16A3']
GLYCOLYSIS     = ['HK1', 'HK2', 'PFKL', 'PFKM', 'PFKP', 'PKM', 'LDHA', 'LDHB',
                   'GAPDH', 'ENO1', 'PGK1', 'ALDOA', 'TPI1', 'PGAM1', 'GPI',
                   'SLC16A1', 'SLC16A3']

def score_patient(pat_dict, genes):
    vals = [float(pat_dict[g]) for g in genes
            if g in pat_dict and pat_dict[g] is not None and str(pat_dict[g]) != 'NA']
    return float(np.mean(vals)) if vals else None

def spear(x, y):
    """Spearman correlation, returns (rho, p)"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 5:
        return None, None
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p)

def partial_spearman(r1, r2, ctrl):
    """Partial Spearman: rho_xy_z = (rho_xy - rho_xz*rho_yz) / sqrt((1-rho_xz^2)*(1-rho_yz^2))"""
    r1 = np.asarray(r1, dtype=float); r2 = np.asarray(r2, dtype=float)
    ctrl = np.asarray(ctrl, dtype=float)
    m = ~(np.isnan(r1) | np.isnan(r2) | np.isnan(ctrl))
    if m.sum() < 10:
        return None, None
    r1, r2, ctrl = rankdata(r1[m]), rankdata(r2[m]), rankdata(ctrl[m])
    rho_xy = float(stats.spearmanr(r1, r2)[0])
    rho_xz = float(stats.spearmanr(r1, ctrl)[0])
    rho_yz = float(stats.spearmanr(r2, ctrl)[0])
    denom = math.sqrt(max((1 - rho_xz**2) * (1 - rho_yz**2), 1e-10))
    r_par = (rho_xy - rho_xz * rho_yz) / denom
    n = len(r1); df = max(n - 3, 1)
    t = r_par * math.sqrt(df / max(1 - r_par**2, 1e-10))
    p_par = 2 * (1 - stats.t.cdf(abs(t), df))
    return float(r_par), float(p_par)

# ---------------------------------------------------------------
# Part A: 32-cancer pathway correlations + partial correlations
# ---------------------------------------------------------------
results = {}
for cancer in sorted(expr.keys()):
    pats = expr[cancer].get('patients', {})
    rows = []
    for pid, d in pats.items():
        if 'AARS2' not in d:
            continue
        a2  = score_patient(d, ['AARS2'])
        ls  = score_patient(d, LAC_SUBSTRATE)
        lf  = score_patient(d, LAC_FLUX)
        hs  = score_patient(d, HSP_GENES)
        hi  = score_patient(d, HIF1A_TARGETS)
        gl  = score_patient(d, GLYCOLYSIS)
        if all(v is not None for v in [a2, ls, lf, hs, hi, gl]):
            rows.append({'AARS2': a2, 'lac_sub': ls, 'lac_flux': lf,
                         'hsp': hs, 'hif1a': hi, 'glyc': gl})
    n = len(rows)
    if n < 10:
        continue
    a  = np.array([r['AARS2'] for r in rows])
    ls = np.array([r['lac_sub'] for r in rows])
    lf = np.array([r['lac_flux'] for r in rows])
    hs = np.array([r['hsp'] for r in rows])
    hi = np.array([r['hif1a'] for r in rows])
    gl = np.array([r['glyc'] for r in rows])
    res = {'n': n}
    # Simple Spearman
    for label, x in [('lac_sub', ls), ('lac_flux', lf), ('hsp', hs),
                      ('hif1a', hi), ('glyc', gl)]:
        rho, pval = spear(a, x)
        res[f'{label}_rho'] = rho; res[f'{label}_p'] = pval
    # Partial Spearman
    for (l1, x1), ctrl, lbl in [
        (('hsp', hs), gl, 'hsp|glyc'),
        (('hsp', hs), hi, 'hsp|hif1a'),
        (('lac_sub', ls), gl, 'lac_sub|glyc'),
        (('lac_sub', ls), hi, 'lac_sub|hif1a'),
        (('glyc', gl), hs, 'glyc|hsp'),
    ]:
        rp, pp = partial_spearman(a, x1, ctrl)
        res[f'{l1}_{lbl}_rho'] = rp; res[f'{l1}_{lbl}_p'] = pp
    results[cancer] = res

# FDR per metric across all cancers
for metric in ['lac_sub', 'lac_flux', 'hsp', 'hif1a', 'glyc']:
    ps = [(c, results[c][f'{metric}_p']) for c in results
           if f'{metric}_p' in results[c] and results[c][f'{metric}_p'] is not None]
    if ps:
        _, qs, _, _ = multipletests([p for _, p in ps], method='fdr_bh')
        for (c, _), q in zip(ps, qs):
            results[c][f'{metric}_q'] = float(q)
for metric in ['hsp|glyc', 'hsp|hif1a', 'lac_sub|glyc', 'lac_sub|hif1a', 'glyc|hsp']:
    ps = [(c, results[c].get(f'hsp|glyc_p' if metric == 'hsp|glyc' else
              f'hsp|hif1a_p' if metric == 'hsp|hif1a' else
              f'lac_sub|glyc_p' if metric == 'lac_sub|glyc' else
              f'lac_sub|hif1a_p' if metric == 'lac_sub|hif1a' else
              f'glyc|hsp_p', None))
           for c in results]
    valid = [(c, p) for c, p in ps if p is not None and 0 < p < 1]
    if valid:
        _, qs, _, _ = multipletests([p for _, p in valid], method='fdr_bh')
        for (c, _), q in zip(valid, qs):
            key = ('hsp|glyc' if metric == 'hsp|glyc' else
                   'hsp|hif1a' if metric == 'hsp|hif1a' else
                   'lac_sub|glyc' if metric == 'lac_sub|glyc' else
                   'lac_sub|hif1a' if metric == 'lac_sub|hif1a' else 'glyc|hsp')
            results[c][f'{key}_q'] = float(q)

out_a = os.path.join(BASE, 'results', 'revision_pathway_corr_fdr.json')
with open(out_a, 'w', encoding='utf-8') as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print('Part A done: revision_pathway_corr_fdr.json')

# Print key table
print('\n=== GI cancers key results ===')
print(f'{"Cancer":8} {"n":4} {"lac_sub_r":10} {"lac_sub_q":10} {"lac_flux_r":10} {"lac_flux_q":10} {"hsp_r":8} {"hsp_q":8} {"hsp|glyc_r":10} {"hsp|glyc_q":10}')
for c in ['STAD', 'ESCA', 'LIHC', 'PAAD', 'COADREAD']:
    if c not in results: continue
    r = results[c]
    def fmt(v): return f'{v:.3f}' if v is not None else 'NA'
    def fe(v): return f'{v:.2e}' if v is not None else 'NA'
    print(f'{c:8} {r["n"]:4} '
          f'{fmt(r.get("lac_sub_rho")):10} {fe(r.get("lac_sub_q")):10} '
          f'{fmt(r.get("lac_flux_rho")):10} {fe(r.get("lac_flux_q")):10} '
          f'{fmt(r.get("hsp_rho")):8} {fe(r.get("hsp_q")):8} '
          f'{fmt(r.get("hsp|glyc_rho")):10} {fe(r.get("hsp|glyc_q")):10}')

# ---------------------------------------------------------------
# Part B: Cox continuous + RCS for 5 GI cancers
# ---------------------------------------------------------------
try:
    from lifelines import CoxPHFitter
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False

if HAS_LIFELINES:
    surv_raw = load_json(os.path.join(DATA, 'cbio_aars2_survival_scipy.json'))
    # Format: {cancer: {patients: {pid: {OS_MONTHS: float, OS_STATUS: 'DECEASED'|'LIVING'}}}}
    if isinstance(surv_raw, dict):
        surv = surv_raw
    else:
        print('survival data format unexpected:', type(surv_raw)); surv = {}

    cox_results = {}
    for cancer in ['STAD', 'ESCA', 'LIHC', 'PAAD', 'COADREAD']:
        if cancer not in surv:
            continue
        pats = surv[cancer].get('patients', {})
        rows = []
        for pid, d in pats.items():
            if 'OS_MONTHS' not in d or 'OS_STATUS' not in d:
                continue
            t = float(d['OS_MONTHS'])
            e = 1 if str(d['OS_STATUS']) == 'DECEASED' else 0
            rows.append({'pid': pid, 'duration': t, 'event': e})
        if len(rows) < 20:
            continue
        df_cox = type('DataFrame', (), {'__init__': lambda s: None})()

        # Build DataFrame
        import pandas as pd
        df_cox = pd.DataFrame(rows)
        df_cox['AARS2'] = [expr.get(cancer, {}).get('patients', {}).get(r['pid'], {}).get('AARS2')
                            for r in rows]
        df_cox = df_cox.dropna(subset=['AARS2', 'duration', 'event'])
        n = len(df_cox)
        if n < 20:
            continue

        # Cox PH: AARS2 as continuous variable
        cph = CoxPHFitter()
        cph.fit(df_cox[['duration', 'event', 'AARS2']], duration_col='duration', event_col='event')
        hr = float(cph.hazard_ratios_['AARS2'])
        ci_low = float(cph.confidence_intervals_.loc['AARS2', '95% lower-bound'])
        ci_high = float(cph.confidence_intervals_.loc['AARS2', '95% upper-bound'])
        p_cox = float(cph.summary.loc['AARS2', 'p'])
        cox_results[cancer] = {
            'n': n, 'hr': hr, 'ci_low': ci_low, 'ci_high': ci_high, 'p_cox': p_cox,
            'events': int(df_cox['event'].sum())
        }

        # Median-split KM for reference
        med = df_cox['AARS2'].median()
        df_cox['group'] = np.where(df_cox['AARS2'] >= med, 'High', 'Low')
        km = df_cox.groupby('group').apply(lambda g: pd.Series({
            'n': len(g),
            'events': g['event'].sum(),
            'median_os': g['duration'].median()
        })).unstack()

        print(f'\n{cancer} Cox PH: n={n} events={df_cox["event"].sum()} HR={hr:.3f} '
              f'95%CI=[{ci_low:.3f},{ci_high:.3f}] p={p_cox:.4f}')
        print(f'  Median-split: High n={len(df_cox[df_cox.group=="High"])} '
              f'Low n={len(df_cox[df_cox.group=="Low"])}')

    out_cox = os.path.join(BASE, 'results', 'revision_cox_survival.json')
    with open(out_cox, 'w', encoding='utf-8') as f:
        json.dump(cox_results, f, ensure_ascii=False, indent=2)
    print('\nPart B done: revision_cox_survival.json')
else:
    print('lifelines not available; skipping Part C')

print('\nALL DONE')
