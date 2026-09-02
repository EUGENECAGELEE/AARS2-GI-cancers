"""R12: Cox PH continuous + tertile AARS2 for 5 GI cancers (lifelines).
Merges patient-level OS (cBioPortal clinical data) with sample-level AARS2 expression.
"""
import urllib.request, json, ssl, os, time
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

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
all_out = {}
for cancer, study in GI_STUDIES.items():
    print(f'=== {cancer} ===', flush=True)
    try:
        # 1. sample-level AARS2 expression
        prof_id = f'{study}_rna_seq_v2_mrna'
        slists = get_json(f'https://www.cbioportal.org/api/studies/{study}/sample-lists')
        sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower() and 'all' in s['sampleListId'].lower()), None)
        if not sl:
            sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower()), slists[0])
        expr_raw = post_json(f'https://www.cbioportal.org/api/molecular-profiles/{prof_id}/molecular-data/fetch',
                             {'entrezGeneIds': [57505], 'sampleListId': sl['sampleListId']})
        # expr entries carry sampleId + patientId directly
        pat_expr = {}
        for e in expr_raw:
            val = e.get('value')
            if val is None:
                continue
            pid = e.get('patientId', '')
            if pid:
                pat_expr.setdefault(pid, []).append(float(val))
        pat_expr = {k: float(np.mean(v)) for k, v in pat_expr.items()}

        # 2. patient-level clinical OS
        clin_p = get_json(f'https://www.cbioportal.org/api/studies/{study}/clinical-data?clinicalDataType=PATIENT')
        os_months, os_status = {}, {}
        for c in clin_p:
            aid, v = c.get('clinicalAttributeId'), c.get('value')
            if not v or v == 'NA':
                continue
            pid = c.get('patientId', '')
            if aid == 'OS_MONTHS':
                os_months[pid] = float(v)
            elif aid == 'OS_STATUS':
                os_status[pid] = v

        # 3. merge
        rows = []
        for pid, expr_val in pat_expr.items():
            if pid not in os_months or pid not in os_status:
                continue
            t = os_months[pid]
            e = 1 if str(os_status[pid]).upper() == 'DECEASED' or 'DECEASED' in str(os_status[pid]).upper() else 0
            rows.append({'patient': pid, 'duration': t, 'event': e, 'AARS2': expr_val})
        df = pd.DataFrame(rows).dropna(subset=['duration', 'event', 'AARS2'])
        df = df[df['duration'] > 0]
        n = len(df)
        if n < 20:
            print(f'  only {n} patients; skip')
            continue
        events = int(df['event'].sum())

        # Cox continuous
        cph = CoxPHFitter()
        cph.fit(df[['duration', 'event', 'AARS2']], duration_col='duration', event_col='event')
        hr = float(cph.hazard_ratios_['AARS2'])
        ci = cph.confidence_intervals_.loc['AARS2']
        p_cox = float(cph.summary.loc['AARS2', 'p'])

        # Cox per-SD (standardized)
        df_sd = df.copy()
        df_sd['AARS2_z'] = (df_sd['AARS2'] - df_sd['AARS2'].mean()) / df_sd['AARS2'].std()
        cph_z = CoxPHFitter()
        cph_z.fit(df_sd[['duration', 'event', 'AARS2_z']], duration_col='duration', event_col='event')
        hr_z = float(cph_z.hazard_ratios_['AARS2_z'])
        ci_z = cph_z.confidence_intervals_.loc['AARS2_z']
        p_z = float(cph_z.summary.loc['AARS2_z', 'p'])

        # Tertile groups
        t1, t2 = df['AARS2'].quantile([1/3, 2/3])
        df['group'] = pd.cut(df['AARS2'], bins=[-np.inf, t1, t2, np.inf], labels=['Low', 'Mid', 'High'])
        df['group_ord'] = df['group'].cat.codes
        # Cox with tertile as ordinal
        cph_t = CoxPHFitter()
        cph_t.fit(df[['duration', 'event', 'group_ord']], duration_col='duration', event_col='event')
        hr_tert = float(cph_t.hazard_ratios_['group_ord'])
        p_tert = float(cph_t.summary.loc['group_ord', 'p'])

        # High vs Low log-rank
        hi = df[df['group'] == 'High']
        lo = df[df['group'] == 'Low']
        lr = logrank_test(hi['duration'], lo['duration'], hi['event'], lo['event'])

        all_out[cancer] = {
            'n': n, 'events': events,
            'cox_hr': hr, 'cox_ci_low': float(ci['95% lower-bound']),
            'cox_ci_high': float(ci['95% upper-bound']), 'cox_p': p_cox,
            'cox_hr_perSD': hr_z, 'cox_ci_low_perSD': float(ci_z['95% lower-bound']),
            'cox_ci_high_perSD': float(ci_z['95% upper-bound']), 'cox_p_perSD': p_z,
            'tertile_hr': hr_tert, 'tertile_p': p_tert,
            'hi_vs_lo_logrank_p': float(lr.p_value),
            'hi_n': int(len(hi)), 'lo_n': int(len(lo)),
            'median_AARS2': float(df['AARS2'].median()),
            'q1': float(df['AARS2'].quantile(0.25)), 'q3': float(df['AARS2'].quantile(0.75)),
        }
        print(f'  n={n} events={events} Cox HR={hr:.3f} (95%CI {ci["95% lower-bound"]:.3f}-{ci["95% upper-bound"]:.3f}) p={p_cox:.4f} | perSD HR={hr_z:.3f} p={p_z:.4f} | tertile p={p_tert:.4f} | HiLo logrank p={lr.p_value:.4f}')
    except Exception as e:
        print(f'  FAIL {e}', flush=True)

with open(os.path.join(BASE, 'results', 'revision_cox_survival.json'), 'w', encoding='utf-8') as f:
    json.dump(all_out, f, ensure_ascii=False, indent=2)
print('\nDONE: revision_cox_survival.json')
