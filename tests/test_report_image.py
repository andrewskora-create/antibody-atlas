import base64
import unittest
from io import BytesIO
from PIL import Image
from app.report_image import report_image

class ReportImageTests(unittest.TestCase):
    def payload(self):
        buf=BytesIO();Image.new('RGB',(900,680),'white').save(buf,format='PNG')
        return {'png':'data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode(),'structure_sha256':'a'*64}
    def test_selected_model_png_accepted(self):
        self.assertTrue(report_image(self.payload(),'a'*64).startswith(b'\x89PNG'))
        self.assertIsNone(report_image(None,None))
    def test_wrong_model_and_invalid_data_rejected(self):
        with self.assertRaises(ValueError):report_image(self.payload(),'b'*64)
        for value in ('data:image/png;base64,invalid','data:image/jpeg;base64,AA=='):
            p=self.payload();p['png']=value
            with self.assertRaises(ValueError):report_image(p,'a'*64)
