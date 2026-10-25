import io
import zipfile
from unittest.mock import patch
from fastapi.testclient import TestClient
from app import main as api, batch_web


def test_selected_report_names_preserve_commas_unicode_and_percent():
    names=['Parent, reference','Pair α','50% candidate']
    def record(job,name):
        return {'name':name} if name in names else None
    with patch.object(batch_web,'get',return_value={'rows':[{'name':n} for n in names]}), \
         patch.object(batch_web,'record_for',side_effect=record), \
         patch('app.report.make_report',side_effect=lambda r,title: b'%PDF-'+r['name'].encode()):
        client=TestClient(api.app)
        response=client.get('/api/batch/job/reports',params=[('name',n) for n in names]+[('name',names[0])])
        assert response.status_code==200
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            assert archive.namelist()==[n+'_developability.pdf' for n in names]
            assert 'NOT_INCLUDED.txt' not in archive.namelist()
        legacy=client.get('/api/batch/job/reports',params={'names':'Pair α,50% candidate'})
        assert legacy.status_code==200
        assert len(zipfile.ZipFile(io.BytesIO(legacy.content)).namelist())==2


def test_report_keeps_molecular_basis_and_requested_numbering():
    from app import report
    data={'scheme':'kabat','ph':6.0,'chains':{},'molecule_flags':{
        'basis':'full_length','assembly':{'format':'H2L2'},
        'fv':{'pi':8.32,'net_charge_assay_ph':24.78}}}
    with patch.object(report,'Paragraph',wraps=report.Paragraph) as paragraph:
        pdf=report.make_report(data)
    drawn=[str(call.args[0]) for call in paragraph.call_args_list]
    assert pdf.startswith(b'%PDF')
    for text in ['Assembled H2L2 pI','Assembled H2L2 charge (e)','8.32','24.78','Region / KABAT']:
        assert text in drawn
    assert 'Region / IMGT' not in drawn


def test_empty_priority_report_does_not_overlap_methods():
    from app import report
    drawn=[]
    original=report.Paragraph
    class RecordedParagraph(original):
        def drawOn(self, canvas, x, y, *args, **kwargs):
            drawn.append((self.text, y, y+self.height))
            return super().drawOn(canvas,x,y,*args,**kwargs)
    with patch.object(report,'Paragraph',RecordedParagraph):
        report.make_report({'chains':{}})
    empty=next(b for b in drawn if b[0].startswith('No assessed high or medium sites'))
    method=next(b for b in drawn if b[0].startswith('Method: Sequence motifs'))
    assert method[2] < empty[1]
