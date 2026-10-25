from html.parser import HTMLParser
from types import SimpleNamespace
from urllib.parse import unquote
from app.batch_dashboard import render


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []
    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.hrefs.append(dict(attrs)['href'])


def test_punctuation_stays_inside_relative_download_filename():
    name = 'Molecule "A" #2 & sample'
    state = {'options': {'pdf': True, 'fold': True},
             'molecules': {name: {'status': 'complete', 'row': {'name': name}}}}
    links = Links()
    links.feed(render(state, SimpleNamespace(warnings=[])))
    assert [unquote(href) for href in links.hrefs] == [
        name + '_brief.pdf', name + '_report.pdf', name + '.json', name + '.pdb']
    assert all('#' not in href and '"' not in href for href in links.hrefs)


def test_no_pdf_run_does_not_offer_missing_pdf_files():
    state = {'options': {'pdf': False}, 'molecules': {'A': {'status': 'complete'}}}
    links = Links()
    links.feed(render(state, SimpleNamespace(warnings=[])))
    assert links.hrefs == ['A.json']
