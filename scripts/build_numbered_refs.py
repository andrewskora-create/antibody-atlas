"""Build IMGT-positioned OGRDB V-gene reference for exact framework audit."""
import json
from pathlib import Path
from anarcii import Anarcii

root = Path(__file__).resolve().parents[1]
source = json.loads((root / 'data/human_v_ogrdb.json').read_text())
model = Anarcii(seq_type='antibody', cpu=True, ncpu=2, verbose=False)
refs = source['references']
numbered = model.number({f'ref{i}': item['sequence'] for i, item in enumerate(refs)})
out = []
for i, item in enumerate(refs):
    result = numbered[f'ref{i}']
    if result.get('error') or not result.get('numbering'):
        continue
    positions = {f'{pos}{ins.strip()}': aa for (pos, ins), aa in result['numbering'] if aa != '-' and pos <= 104}
    out.append({'name': item['name'], 'locus': item['locus'], 'positions': positions})
output = {'method': 'ANARCII IMGT numbering, positions 1-104 of functional OGRDB V alleles',
          'source': source['source'], 'source_url': source['source_url'],
          'release_versions': source['release_versions'], 'references': out}
(root / 'data/human_v_ogrdb_numbered.json').write_text(json.dumps(output, separators=(',', ':')))
print(f'Numbered {len(out)}/{len(refs)} references')
