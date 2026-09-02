"""Revision R1: control genes expression fetch (mitochondrial markers + aaRS paralogs + TCA)
Fetch via cBioPortal molecular-data/fetch (same pipeline as pan_cancer_fetch.py).
"""
import urllib.request, json, ssl, os, time

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

def get_json(url, retries=4):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=90, context=ctx) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception:
            if i == retries - 1:
                raise
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
            if i == retries - 1:
                raise
            time.sleep(8)

CONTROL_IDS = {"COX4I1": 1327, "COX5A": 9377, "TOMM20": 9804, "TOMM22": 56993,
               "DARS2": 55157, "EARS2": 124454, "GARS1": 2617, "LARS2": 23395,
               "MARS2": 92935, "TARS2": 80222, "YARS2": 51067, "SARS2": 54938,
               "PDHA1": 5160, "SDHB": 6390, "MDH2": 4191, "CS": 1431,
               "IDH2": 3418, "SUCLA2": 8803, "OGDH": 4967, "DLST": 1743,
               "FH": 2271, "TFAM": 7019, "HSPA9": 3313, "HSPD1": 3329}

STUDIES = [
    'acc_tcga_pan_can_atlas_2018','blca_tcga_pan_can_atlas_2018','brca_tcga_pan_can_atlas_2018',
    'cesc_tcga_pan_can_atlas_2018','chol_tcga_pan_can_atlas_2018','coadread_tcga_pan_can_atlas_2018',
    'dlbc_tcga_pan_can_atlas_2018','esca_tcga_pan_can_atlas_2018','gbm_tcga_pan_can_atlas_2018',
    'hnsc_tcga_pan_can_atlas_2018','kich_tcga_pan_can_atlas_2018','kirc_tcga_pan_can_atlas_2018',
    'kirp_tcga_pan_can_atlas_2018','laml_tcga_pan_can_atlas_2018','lgg_tcga_pan_can_atlas_2018',
    'lihc_tcga_pan_can_atlas_2018','luad_tcga_pan_can_atlas_2018','lusc_tcga_pan_can_atlas_2018',
    'meso_tcga_pan_can_atlas_2018','ov_tcga_pan_can_atlas_2018','paad_tcga_pan_can_atlas_2018',
    'pcpg_tcga_pan_can_atlas_2018','prad_tcga_pan_can_atlas_2018','sarc_tcga_pan_can_atlas_2018',
    'skcm_tcga_pan_can_atlas_2018','stad_tcga_pan_can_atlas_2018','tgct_tcga_pan_can_atlas_2018',
    'thca_tcga_pan_can_atlas_2018','thym_tcga_pan_can_atlas_2018','ucs_tcga_pan_can_atlas_2018',
    'ucec_tcga_pan_can_atlas_2018','uvm_tcga_pan_can_atlas_2018',
]

BASE = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'
gene_ids = list(CONTROL_IDS.values())
out = {}
for study in STUDIES:
    cancer = study.split('_')[0].upper()
    try:
        prof_id = f'{study}_rna_seq_v2_mrna'
        slists = get_json(f'https://www.cbioportal.org/api/studies/{study}/sample-lists')
        sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower() and 'all' in s['sampleListId'].lower()), None)
        if not sl:
            sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower()), slists[0])
        expr_raw = post_json(f'https://www.cbioportal.org/api/molecular-profiles/{prof_id}/molecular-data/fetch',
                             {'entrezGeneIds': gene_ids, 'sampleListId': sl['sampleListId']})
        # aggregate per patient
        pat = {}
        for e in expr_raw:
            pid = e['patientId']
            gid = e['entrezGeneId']
            val = e.get('value')
            if val is None or val == 'NA':
                continue
            sym = {v: k for k, v in CONTROL_IDS.items()}[gid]
            pat.setdefault(pid, {})[sym] = float(val)
        out[cancer] = {'study': study, 'n': len(pat), 'patients': pat}
        print(f'{cancer}: {len(pat)} patients', flush=True)
    except Exception as e:
        print(f'{cancer}: FAIL {e}', flush=True)
    time.sleep(0.3)

os.makedirs(os.path.join(BASE, 'data', 'pancancer'), exist_ok=True)
with open(os.path.join(BASE, 'data', 'pancancer', 'control_expr_matrix.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f)
print('DONE: control_expr_matrix.json')
