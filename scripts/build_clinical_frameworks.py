"""Build a frozen, offline clinical sequence reference from public TheraSAbDab CSV."""
from pathlib import Path
import csv, hashlib, json, sys
from collections import Counter
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from app.numbering import _model, CDR_RANGES

SOURCE_URL = 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv'
STAGE_COLUMN = "Highest_Clin_Trial (Feb '25)"

def build(source, destination):
    source = Path(source)
    rows = list(csv.DictReader(source.read_text(encoding='utf-8-sig').splitlines()))
    eligible = [r for r in rows if r[STAGE_COLUMN] in {'Approved','Phase-III'}]
    candidates, excluded = [], []
    for row in eligible:
        for column, family in [('HeavySequence','H'),('LightSequence','L'),('HeavySequence(ifbispec)','H'),('LightSequence(ifbispec)','L')]:
            seq = row[column].strip().upper()
            if seq in {'','NA','NONE'}: continue
            if not set(seq) <= set('ACDEFGHIKLMNPQRSTVWY'):
                excluded.append({'drug':row['Therapeutic'],'column':column,'reason':'Noncanonical or ambiguous sequence characters'})
                continue
            candidates.append({'drug':row['Therapeutic'],'stage':row[STAGE_COLUMN], 'development_status':row['Est. Status'], 'reference_chain':column,'family':family,'sequence':seq})
    sequences = {str(i):s for i,s in enumerate(dict.fromkeys(c['sequence'] for c in candidates))}
    model = _model()
    print(f'Numbering {len(sequences)} unique sequences from {len(eligible)} eligible drugs', flush=True)
    imgt = model.number(sequences)
    kabat = model.to_scheme('kabat')
    prepared = {}
    for key,seq in sequences.items():
        schemes = {}
        for scheme, records in [('imgt',imgt),('kabat',kabat)]:
            record = records.get(key,{})
            if record.get('error') or not record.get('numbering'): break
            ct = record.get('chain_type')
            if ct not in {'H','K','L'}: break
            start = int(record.get('query_start') or 0)
            residues=[]
            pos=start
            for (number,insertion),aa in record['numbering']:
                if aa=='-':continue
                pos+=1
                region=next((f'CDR{i}' for i,(lo,hi) in enumerate(CDR_RANGES[scheme]['H' if ct=='H' else 'L'],1) if lo<=number<=hi),'Framework')
                residues.append({'sequence_position':pos,'number':f'{ct}{number}{insertion.strip()}','aa':aa,'region':region})
            if ''.join(r['aa'] for r in residues) != seq or start != 0: break
            schemes[scheme]={'chain_type':ct,'scheme':scheme,'residues':residues}
        if len(schemes)==2:prepared[seq]=schemes
    references=[]
    for item in candidates:
        numbered=prepared.get(item['sequence'])
        if not numbered or ('H' if numbered['imgt']['chain_type']=='H' else 'L') != item['family']:
            excluded.append({'drug':item['drug'],'column':item['reference_chain'],'reason':'Incomplete, unrecognized, multi-domain or incompatible variable-domain numbering'})
            continue
        references.append({**item,'numbering':numbered})
    result={'version':'2026-10-04-clinical-framework-v1','source_url':SOURCE_URL,'source_snapshot':'Highest_Clin_Trial (Feb 2025)', 'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'stage_column':STAGE_COLUMN,
            'coverage':{'eligible_drugs':len(eligible),'reference_chains':len(references),'approved_drugs':sum(r[STAGE_COLUMN]=='Approved' for r in eligible),'phase_iii_drugs':sum(r[STAGE_COLUMN]=='Phase-III' for r in eligible),'represented_drugs':len({r['drug'] for r in references}),'excluded_chains':len(excluded)},'excluded':excluded,'references':references}
    Path(destination).write_text(json.dumps(result,separators=(',',':')))
    print(json.dumps(result['coverage']),flush=True)

if __name__=='__main__':build(sys.argv[1],sys.argv[2])
