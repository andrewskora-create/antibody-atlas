import unittest
from unittest.mock import patch
from app import germline_cdr as g


def fixtures(monkeypatch, chain='H', extra=None):
    positions = {'27':'N','28':'S','38':'Y'}
    ref = {'name':'IGHV1*01','locus':'IGH','positions':dict(positions)}
    if extra: ref['positions'].update(extra)
    monkeypatch.setattr(g, '_references', lambda: {'source':'test functional library','source_url':'https://example.org','release_versions':{'IGH':1,'IGK':1,'IGL':1},'references':[ref]})
    rows = [{'number':chain+k,'aa':v,'sequence_position':i+1,'region':'CDR1'} for i,(k,v) in enumerate(positions.items())]
    n = {'chain_type':chain,'scheme':'imgt','residues':rows}
    return {'position':1,'motif':'NS'},n


class GermlineCDRTests(unittest.TestCase):
    def setUp(self):
        self.patches=[]
    def tearDown(self):
        for p in self.patches: p.stop()
    def setattr(self, obj, key, value):
        p=patch.object(obj,key,value);p.start();self.patches.append(p)
    def test_full_cdr_exact_and_all_alleles(self):
        site,n=fixtures(self)
        data=g._references();data['references'].append(dict(data['references'][0],name='IGHV2*02'))
        self.setattr(g,'_references',lambda:data)
        out=g.annotate_germline_cdr(site,n,n)
        assert out['alleles']==['IGHV1*01','IGHV2*02']
        assert out['status']=='exact_match' and out['motif_matched']
        assert out['sequence']=='NSY'


    def test_motif_only_does_not_establish_whole_cdr_identity(self):
        site,n=fixtures(self);n['residues'][-1]['aa']='F'
        assert g.annotate_germline_cdr(site,n,n)['status']=='no_exact_match'


    def test_extra_reference_insertion_prevents_match(self):
        site,n=fixtures(self,extra={'38A':'G'})
        assert g.annotate_germline_cdr(site,n,n)['status']=='no_exact_match'


    def test_subtype_is_not_interchangeable(self):
        site,n=fixtures(self,chain='K')
        assert g.annotate_germline_cdr(site,n,n)['status']=='no_exact_match'


    def test_cdr3_cannot_be_inferred_from_v_only_library(self):
        site,n=fixtures(self);n['residues'][0].update(number='H105',region='CDR3')
        out=g.annotate_germline_cdr(site,n,n)
        assert out['status']=='unassessed' and not out['alleles']
        assert 'V-only' in out['summary']


    def test_selected_cdr_canonical_framework_is_unassessed(self):
        site,n=fixtures(self);n['residues'][0].update(number='H26',region='Framework')
        selected={'residues':[dict(r,region='CDR1') for r in n['residues']]}
        assert g.annotate_germline_cdr(site,n,selected)['status']=='unassessed'


    def test_boundary_motif_not_called_germline_motif(self):
        site,n=fixtures(self);site.update(position=3,motif='YA')
        n['residues'].append({'number':'H39','aa':'A','sequence_position':4,'region':'Framework'})
        out=g.annotate_germline_cdr(site,n,n)
        assert out['status']=='exact_match' and not out['motif_matched']
