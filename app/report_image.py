"""Validate a local browser rendering before embedding it in the PDF."""
import base64
import binascii
from io import BytesIO
from PIL import Image


def report_image(payload, expected_hash):
    if payload is None:
        return None
    if not expected_hash or payload.get('structure_sha256') != expected_hash:
        raise ValueError('Structure image does not match the selected model; refresh and export again')
    encoded=payload.get('png','')
    if not encoded.startswith('data:image/png;base64,') or len(encoded)>6_000_000:
        raise ValueError('Structure image must be a PNG smaller than 4 MB')
    try:
        raw=base64.b64decode(encoded.split(',',1)[1],validate=True)
        with Image.open(BytesIO(raw)) as image:
            if image.format!='PNG' or not (64<=image.width<=2400 and 64<=image.height<=2400):
                raise ValueError('Invalid structure image dimensions')
            image.verify()
    except (binascii.Error,OSError,Image.DecompressionBombError) as exc:
        raise ValueError('Invalid structure PNG') from exc
    return raw
