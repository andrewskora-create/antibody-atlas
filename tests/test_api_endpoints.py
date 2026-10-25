"""Exercise the new endpoints over HTTP, not just the modules behind them.

This exists because two bugs shipped that the module tests could not see.
``/api/structure`` passed ``model_index`` and ``sequences`` to a ``prepare()``
that accepted neither, and ``/api/structure/inspect`` called ``.as_dict()`` on
a report class that did not define it. Both modules were well covered and both
endpoints were broken on the first real request: the wiring between them was
what nobody tested.

So these tests do the least interesting and most valuable thing -- call each
route the way a browser does and check it returns something usable.
"""
from __future__ import annotations

import unittest
from pathlib import Path

from fastapi.testclient import TestClient

from app.main import app

ROOT = Path(__file__).resolve().parents[1]
FV = ROOT / "tests" / "fixtures" / "example_fv.pdb"
#: A real antibody/antigen complex, needed only for the epitope path. It is
#: fetched by the co-crystal validation and is not redistributed, so the
#: antigen tests skip without it rather than failing.
COCRYSTAL = ROOT / "validation" / ".cocrystal_work" / "9pi9.cif"

HEAVY = ("EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEG"
         "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS")
LIGHT = ("DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSG"
         "SGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIK")


def client() -> TestClient:
    # Let server errors surface as 500s so a broken route is a failed
    # assertion rather than an exception escaping the test.
    return TestClient(app, raise_server_exceptions=False)


class StructureInspect(unittest.TestCase):
    def test_inspect_returns_a_serialisable_report(self):
        """The regression: this route called a method that did not exist."""
        if not COCRYSTAL.is_file():
            self.skipTest("needs the co-crystal; run validation/paratope_cocrystal.py")
        with COCRYSTAL.open("rb") as handle:
            response = client().post("/api/structure/inspect",
                                     files={"file": ("9pi9.cif", handle, "chemical/x-cif")})
        self.assertEqual(response.status_code, 200, response.text)
        body = response.json()
        for key in ("chains", "fv_pairs", "antigen_chains", "copies", "has_antigen"):
            self.assertIn(key, body)
        self.assertTrue(body["chains"])
        self.assertEqual(body["copies"], len(body["fv_pairs"]))

    def test_inspect_rejects_a_file_that_is_not_a_structure(self):
        response = client().post("/api/structure/inspect",
                                 files={"file": ("notes.txt", b"not a structure", "text/plain")})
        self.assertEqual(response.status_code, 422, response.text)


class StructureUpload(unittest.TestCase):
    """The regression here was a signature mismatch with the module."""

    def _upload(self, path, heavy, light, **extra):
        with Path(path).open("rb") as handle:
            return client().post(
                "/api/structure",
                files={"file": (Path(path).name, handle, "application/octet-stream")},
                data={"heavy": heavy, "light": light, "ph": "7.4",
                      "scheme": "imgt", **extra})

    def test_antibody_only_structure_is_accepted(self):
        response = self._upload(FV, HEAVY, LIGHT)
        self.assertEqual(response.status_code, 200, response.text[:300])
        job = response.json()
        self.assertEqual(job["status"], "complete")
        self.assertEqual(job["structure_source"], "experimental")
        self.assertIn("ingestion", job)

    def test_experimental_source_travels_with_the_job(self):
        """Nothing downstream may present measured and modelled alike."""
        job = self._upload(FV, HEAVY, LIGHT).json()
        self.assertEqual(job["structure_source"], "experimental")
        fetched = client().get(f"/api/fold/{job['id']}").json()
        self.assertEqual(fetched["structure_source"], "experimental")

    def test_complex_yields_a_measured_epitope(self):
        if not COCRYSTAL.is_file():
            self.skipTest("needs the co-crystal; run validation/paratope_cocrystal.py")
        from app.structure_input import inspect
        report = inspect(COCRYSTAL)
        sequences = {c.chain_id: c.sequence for c in report.chains}
        heavy, light = report.fv_pairs[0]
        response = self._upload(COCRYSTAL, sequences[heavy], sequences[light])
        self.assertEqual(response.status_code, 200, response.text[:300])
        job = response.json()
        epitope = job.get("epitope") or {}
        self.assertEqual(epitope.get("status"), "measured")
        self.assertGreater(epitope.get("contact_count", 0), 0)
        # IMGT insertion order, which plain string sorting gets wrong.
        contacts = epitope["contacts"]
        if "H112A" in contacts and "H112" in contacts:
            self.assertLess(contacts.index("H112A"), contacts.index("H112"))

    def test_a_structure_of_a_different_molecule_is_refused(self):
        """A caveat on the wrong structure is still a wrong answer."""
        response = self._upload(FV, LIGHT, LIGHT)   # light sequence as the heavy
        self.assertIn(response.status_code, (200, 422))
        if response.status_code == 200:
            self.assertFalse((response.json().get("ingestion") or {}).get("usable", True))


class Capabilities(unittest.TestCase):
    def test_reports_what_this_process_actually_serves(self):
        """Static files reload per request; Python does not. The page has to
        be able to tell that the server behind it is older than it is."""
        body = client().get("/api/substitutions/capabilities").json()
        self.assertIn("endpoints", body)
        self.assertIn("features", body)
        self.assertIn("/api/structure", body["endpoints"])
        self.assertTrue(body["features"]["structure_upload"])
        self.assertTrue(body["features"]["repack"])


class Repack(unittest.TestCase):
    def test_repack_runs_against_an_uploaded_structure(self):
        with FV.open("rb") as handle:
            job = client().post("/api/structure",
                                files={"file": ("example_fv.pdb", handle, "chemical/x-pdb")},
                                data={"ph": "7.4", "scheme": "imgt"}).json()
        response = client().post("/api/repack", json={
            "fold_job_id": job["id"], "chain": "H", "position": 35,
            "replacement": "A", "ph": 7.4})
        self.assertEqual(response.status_code, 200, response.text[:300])
        body = response.json()
        self.assertEqual(body["source"], "experimental")
        self.assertIn("hydrogen_bonds", body)
        self.assertIn("no_free_energy", body)

    def test_repack_refuses_a_job_that_does_not_exist(self):
        response = client().post("/api/repack", json={
            "fold_job_id": "0" * 32, "chain": "H", "position": 35,
            "replacement": "A", "ph": 7.4})
        self.assertEqual(response.status_code, 404)


if __name__ == "__main__":
    unittest.main()
