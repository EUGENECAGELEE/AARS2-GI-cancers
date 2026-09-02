"""R4: Fetch HSF1 (3297), CD274 (29126), IDO1 (3620) for 5 GI cancers via cBioPortal POST
(correct endpoint: /api/molecular-profiles/{profile}/molecular-data/fetch)."""
import sys, json, time, ssl, urllib.request
sys.stdout.reconfigure(encoding='utf-8')
base = r'C:\Users\lee_e\Desktop\生信思路\0818AARS2'

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

def get_json(url, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=90, context=ctx) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception as e:
            print(f'  get retry {i+1}: {e}')
            time.sleep(4 + 2 * i)
    return None

def post_json(url, payload, retries=5):
    for i in range(retries):
        try:
            body = json.dumps(payload).encode('utf-8')
            headers = {'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json', 'Content-Type': 'application/json'}
            req = urllib.request.Request(url, headers=headers, data=body, method='POST')
            with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
                return json.loads(r.read().decode('utf-8'))
        except Exception as e:
            print(f'  post retry {i+1}: {e}')
            time.sleep(4 + 3 * i)
    return None

STUDIES = {
    'STAD': 'stad_tcga_pan_can_atlas_2018',
    'ESCA': 'esca_tcga_pan_can_atlas_2018',
    'LIHC': 'lihc_tcga_pan_can_atlas_2018',
    'PAAD': 'paad_tcga_pan_can_atlas_2018',
    'COADREAD': 'coadread_tcga_pan_can_atlas_2018',
}
GENES = {'HSF1': 3297, 'CD274': 29126, 'IDO1': 3620}

out = {}
for cancer, study in STUDIES.items():
    print(f'== {cancer} ==')
    slists = get_json(f'https://www.cbioportal.org/api/studies/{study}/sample-lists')
    if not slists:
        print('  no sample lists')
        continue
    sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower() and 'all' in s['sampleListId'].lower()), None)
    if not sl:
        sl = next((s for s in slists if 'rna_seq' in s['sampleListId'].lower()), slists[0])
    prof_id = f'{study}_rna_seq_v2_mrna'
    expr_raw = post_json(f'https://www.cbioportal.org/api/molecular-profiles/{prof_id}/molecular-data/fetch',
                         {'entrezGeneIds': list(GENES.values()), 'sampleListId': sl['sampleListId']})
    if not expr_raw:
        print('  FAILED fetch')
        continue
    print(f'  fetched {len(expr_raw)} records, sampleList={sl["sampleListId"]}')
    pat = {}
    rev = {v: k for k, v in GENES.items()}
    for e in expr_raw:
        pid = e.get('patientId')
        gid = e.get('entrezGeneId')
        val = e.get('value')
        if pid is None or gid is None or val is None or val == 'NA':
            continue
        pat.setdefault(pid, {})[rev.get(gid, str(gid))] = float(val)
    print(f'  patients: {len(pat)}')
    out[cancer] = {'study': study, 'sampleListId': sl['sampleListId'], 'patients': pat}
    time.sleep(1)

with open(r'C:\Users\lee_e\Desktop\生信思路\0818AARS2\data\pancancer\hsf1_pdl1_ido1_expr.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, ensure_ascii=False)
print('saved hsf1_pdl1_ido1_expr.json')
