"""R12b: Restricted Cubic Spline (RCS) survival analysis for AARS2 (Reviewer 6, comment 12).
Uses Cox PH with RCS basis for continuous AARS2; tests linearity (non-linear term p).
Fetches fresh patient-level data (duration, event, AARS2).
"""
import urllib.request, json, ssl, os, time
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

def get_json(url, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=90, context=ctx) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception:
            if i == retries - 1: raise
            time.sleep(8)

def post_json(url, payload, retries=4):
    for i in range(retries):
        try:
            body = json.dumps(payload).encode('utf-8')
            headers = {'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json', 'Content-Type': 'application/json'}
            req = urllib.request.Request(url, headers=headers, data=body, method='POST')
            with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception:
            if i == retries - 1: raise
            time.sleep(8)

GI_STUDIES = {
    'STAD': 'stad_tcga_pan_can_atlas_2018',
    'COAD': 'coadread_tcga_pan_can_atlas_2018',
    'ESCA': 'esca_tcga_pan_can_atlas_2018',
    'LIHC': 'lihc_tcga_pan_can_atlas_2018',
    'PAAD': 'paad_tcga_pan_can_atlas_2018',
}

BASE = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'

def fetch_surv_df(cancer, study):
    """Fetch (duration, event, AARS2) per patient; cache to data dir."""
    cache = os.path.join(BASE, 'data', 'pancancer', f'surv_df_{cancer}.csv')
    if os.path.exists(cache):
        return pd.read_csv(cache)
    prof_id = f'{study}_rna_seq_v2_mrna'
    slists = get_json(f'https://www.cbioportal.org/api/studies/{study}/sample-lists')
    sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower() and 'all' in s['sampleListId'].lower()), None)
    if not sl:
        sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower()), slists[0])
    expr_raw = post_json(f'https://www.cbioportal.org/api/molecular-profiles/{prof_id}/molecular-data/fetch',
                         {'entrezGeneIds': [57505], 'sampleListId': sl['sampleListId']})
    pat_expr = {}
    for e in expr_raw:
        val = e.get('value')
        if val is None: continue
        pid = e.get('patientId', '')
        if pid:
            pat_expr.setdefault(pid, []).append(float(val))
    pat_expr = {k: float(np.mean(v)) for k, v in pat_expr.items()}
    clin_p = get_json(f'https://www.cbioportal.org/api/studies/{study}/clinical-data?clinicalDataType=PATIENT')
    os_months, os_status = {}, {}
    for c in clin_p:
        aid, v = c.get('clinicalAttributeId'), c.get('value')
        if not v or v == 'NA': continue
        pid = c.get('patientId', '')
        if aid == 'OS_MONTHS': os_months[pid] = float(v)
        elif aid == 'OS_STATUS': os_status[pid] = v
    rows = []
    for pid, expr_val in pat_expr.items():
        if pid not in os_months or pid not in os_status: continue
        t = os_months[pid]
        e = 1 if 'DECEASED' in str(os_status[pid]).upper() else 0
        rows.append({'duration': t, 'event': e, 'AARS2': expr_val})
    df = pd.DataFrame(rows).dropna()
    df = df[df['duration'] > 0]
    df.to_csv(cache, index=False)
    return df

def rcs_basis(x, knots):
    """Restricted cubic spline basis (Harrell). Returns df with columns rcs1, rcs2..."""
    x = np.asarray(x, dtype=float)
    k = np.asarray(knots, dtype=float)
    nk = len(k)
    # basis: x, (x-k_j)+^3 - (x-k_{nk-1})+^3 * (k_nk - k_j)/(k_nk - k_{nk-1}) + (x-k_nk)+^3*(k_{nk-1}-k_j)/(k_nk-k_{nk-1})
    cols = {'x': x}
    for j in range(1, nk - 1):
        col = np.zeros_like(x)
        for i in range(len(x)):
            xi = x[i]
            term = max(xi - k[j], 0)**3
            term -= max(xi - k[nk-1], 0)**3 * (k[nk-1] - k[j]) / (k[nk-1] - k[nk-2])
            term += max(xi - k[nk-2], 0)**3 * (k[nk-2] - k[j]) / (k[nk-1] - k[nk-2])
            col[i] = term
        cols[f'rcs{j}'] = col
    return pd.DataFrame(cols)

def analyze_rcs(cancer, study, n_knots=4):
    df = fetch_surv_df(cancer, study)
    if len(df) < 40:
        return None
    # knots at quantiles
    quants = np.linspace(0, 1, n_knots)
    knots = np.quantile(df['AARS2'], quants)
    # standardize AARS2 to avoid numerical issues
    mu, sd = df['AARS2'].mean(), df['AARS2'].std()
    df['AARS2_z'] = (df['AARS2'] - mu) / sd
    knots_z = (knots - mu) / sd
    basis = rcs_basis(df['AARS2_z'].values, knots_z)
    df2 = pd.concat([df.reset_index(drop=True), basis], axis=1)
    cols = ['duration', 'event'] + [c for c in basis.columns if c != 'x'] + ['x']
    # Cox with rcs terms
    try:
        cph = CoxPHFitter()
        cph.fit(df2[['duration', 'event', 'x'] + [c for c in basis.columns if c != 'x']],
                duration_col='duration', event_col='event')
        # Non-linearity test: joint test of rcs terms
        summ = cph.summary
        rcs_cols = [c for c in basis.columns if c != 'x']
        # Wald test for all rcs terms = 0
        beta = summ.loc[rcs_cols, 'coef'].values
        cov = cph.variance_matrix_.loc[rcs_cols, rcs_cols].values
        try:
            wald = float(beta @ np.linalg.inv(cov) @ beta)
            p_nl = float(1 - __import__('scipy').stats.chi2.sf(wald, len(rcs_cols)))
        except Exception:
            p_nl = None
        # Overall p for x + rcs
        beta_all = summ.loc[['x'] + rcs_cols, 'coef'].values
        cov_all = cph.variance_matrix_.loc[['x'] + rcs_cols, ['x'] + rcs_cols].values
        try:
            wald_all = float(beta_all @ np.linalg.inv(cov_all) @ beta_all)
            p_all = float(__import__('scipy').stats.chi2.sf(wald_all, len(beta_all)))
        except Exception:
            p_all = None
        # HR at percentiles (10th, 50th, 90th) relative to median
        pcts = df2['x'].quantile([0.1, 0.5, 0.9])
        med_row = df2[['x'] + rcs_cols].iloc[[0]].copy()
        hr_at = {}
        # Use partial hazards: exp(beta @ (basis(x) - basis(median)))
        med_basis = rcs_basis(np.array([pcts[0.5]]), knots_z).iloc[0]
        for pct, val in pcts.items():
            b = rcs_basis(np.array([val]), knots_z).iloc[0]
            diff = b - med_basis
            hr = float(np.exp(beta_all @ diff[['x'] + rcs_cols].values))
            hr_at[f'p{pct}'] = hr
        return {
            'n': len(df), 'events': int(df['event'].sum()),
            'knots': knots.tolist(), 'n_knots': n_knots,
            'p_nonlinear': p_nl, 'p_overall': p_all,
            'hr_p10_vs_median': hr_at.get('p0.1'),
            'hr_p50_vs_median': 1.0,
            'hr_p90_vs_median': hr_at.get('p0.9'),
        }
    except Exception as e:
        print(f'  {cancer} RCS FAIL: {e}')
        return None

out = {}
for cancer, study in GI_STUDIES.items():
    print(f'=== {cancer} ===', flush=True)
    try:
        res = analyze_rcs(cancer, study)
        if res:
            out[cancer] = res
            print(f'  n={res["n"]} events={res["events"]} p_overall={res["p_overall"]} '
                  f'p_nonlinear={res["p_nonlinear"]} HR_p90={res["hr_p90_vs_median"]:.3f}')
    except Exception as e:
        print(f'  FAIL {e}')

with open(os.path.join(BASE, 'results', 'revision_rcs_survival.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
print('\nDONE: revision_rcs_survival.json')
