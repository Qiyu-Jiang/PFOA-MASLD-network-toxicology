# -*- coding: utf-8 -*-
"""run_02 step 1: map UniProt IDs from TargetNet + PharmMapper to gene symbols (UniProt REST)"""
import csv, re, os, requests, time, warnings
warnings.filterwarnings('ignore')

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
RAW = os.path.join(BASE, '02_data_raw')
KEEP = os.path.join(BASE, '03_data_processed')
os.makedirs(KEEP, exist_ok=True)

# collect IDs
ids = set()
tn_rows = list(csv.DictReader(open(os.path.join(RAW, 'TargetNet_targets.csv'), encoding='utf-8-sig')))
pm_rows = list(csv.DictReader(open(os.path.join(RAW, 'PharmMapper_targets.csv'), encoding='utf-8-sig')))
for r in tn_rows:
    v = (r['UniProt_ID'] or '').strip()
    if v and v != 'NONE':
        ids.add(v)
for r in pm_rows:
    v = (r['UniProt'] or '').strip()
    if v and v != 'NONE':
        ids.add(v)
print('unique IDs to map:', len(ids))

S = requests.Session()
S.headers.update({'User-Agent': 'PFOA-NAFLD-project/1.0 (academic use)'})

mapping = {}  # id -> (accession, gene)
def flush_batch(batch):
    q = ' OR '.join(f'id:{x}' if '_' in x and not re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$', x) else f'accession:{x}' for x in batch)
    url = 'https://rest.uniprot.org/uniprotkb/stream'
    params = {'query': f'({q})', 'fields': 'accession,id,gene_primary', 'format': 'tsv'}
    r = S.get(url, params=params, timeout=90)
    r.raise_for_status()
    got = 0
    for line in r.text.splitlines()[1:]:
        parts = line.split('\t')
        if len(parts) >= 2:
            acc, entry = parts[0], parts[1]
            gene = parts[2].split()[0] if len(parts) > 2 and parts[2].strip() else ''
            mapping[acc] = (entry, gene)
            got += 1
    return got

all_ids = sorted(ids)
mapped_count = 0
for i in range(0, len(all_ids), 40):
    batch = all_ids[i:i+40]
    for attempt in range(2):
        try:
            mapped_count += flush_batch(batch)
            break
        except Exception as e:
            print(f'batch {i} attempt {attempt+1} failed:', repr(e)[:100])
            time.sleep(2)
    time.sleep(0.5)
    print(f'progress: {min(i+40, len(all_ids))}/{len(all_ids)}, mapped so far: {mapped_count}')

# entry name (mnemonic) is also a key: build a bidirectional index
by_entry = {v[0]: k for k, v in mapping.items()}
out = os.path.join(KEEP, 'uniprot_mapping.tsv')
with open(out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['UniProt', 'Entry', 'GeneSymbol'])
    for acc, (entry, gene) in sorted(mapping.items()):
        w.writerow([acc, entry, gene])
print('mapping saved:', out, 'records:', len(mapping))

# apply the mapping and produce symbol-level files
def apply_map(value):
    v = value.strip()
    if v in mapping:
        return mapping[v][1]
    if v in by_entry:
        return mapping[by_entry[v]][1]
    return ''

tn_out = os.path.join(KEEP, 'TargetNet_targets_symbols.tsv')
unmapped_tn = 0
with open(tn_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Symbol', 'UniProt_ID', 'Protein', 'Prob'])
    for r in tn_rows:
        s = apply_map(r['UniProt_ID'])
        if not s:
            unmapped_tn += 1
        w.writerow([s, r['UniProt_ID'], r['Protein'], r['Prob']])
print('TargetNet mapped:', len(tn_rows) - unmapped_tn, '/', len(tn_rows), 'unmapped:', unmapped_tn)

pm_out = os.path.join(KEEP, 'PharmMapper_targets_symbols.tsv')
unmapped_pm = 0
with open(pm_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Symbol', 'UniProt', 'Protein_Name', 'Norm_Fit', 'Fit', 'Pharma_Model'])
    for r in pm_rows:
        s = apply_map(r['UniProt'])
        if not s:
            unmapped_pm += 1
        w.writerow([s, r['UniProt'], r['Protein_Name'], r['Norm_Fit'], r['Fit'], r['Pharma_Model']])
print('PharmMapper mapped:', len(pm_rows) - unmapped_pm, '/', len(pm_rows), 'unmapped:', unmapped_pm)

# spot check
for name in ['PPARA', 'PPARG', 'NR1H4', 'FABP4', 'CES1', 'HSD11B1', 'ESR1', 'SRC']:
    hits = [r for r in pm_rows if apply_map(r['UniProt']) == name]
    if hits:
        print('spot', name, '->', hits[0]['UniProt'], hits[0]['Norm_Fit'])
