"""Constant-region assembly, and the full-length basis it gives the flags."""
import unittest

from app.constant_regions import (ConstantRegionError, assemble, available,
                                  default_light_for)
from app.molecule_flags import molecule_flags

HEAVY = ("EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTIS"
         "ADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS")
LIGHT = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
         "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")


def chains(light_type="K"):
    out = {"H": {"sequence": HEAVY,
                 "numbering_imgt": {"scheme": "imgt", "chain_type": "H", "residues": []},
                 "metrics": {"aggregation_proxy": {"patch_count": 0}}}}
    out["L"] = {"sequence": LIGHT,
                "numbering_imgt": {"scheme": "imgt", "chain_type": light_type, "residues": []},
                "metrics": {"aggregation_proxy": {"patch_count": 0}}}
    return out


class ConstantRegionDataTests(unittest.TestCase):
    def test_heavy_constants_are_the_secreted_isoform_not_the_membrane_one(self):
        """The regression that matters: UniProt's displayed gamma isoform is
        membrane-bound, and appending it would add a transmembrane helix."""
        for entry in available()["heavy"]:
            with self.subTest(entry["key"]):
                self.assertLess(entry["length"], 400, "looks like the membrane-bound isoform")
                self.assertGreater(entry["length"], 300)

    def test_shipped_sequences_end_at_the_ch3_c_terminal_lysine(self):
        from app.constant_regions import _reference
        for key, region in _reference()["regions"].items():
            if region["kind"] != "heavy":
                continue
            with self.subTest(key):
                self.assertTrue(region["sequence"].endswith("GK"), region["sequence"][-12:])

    def test_no_hydrophobic_membrane_helix_in_any_shipped_heavy_constant(self):
        from app.constant_regions import _reference
        for key, region in _reference()["regions"].items():
            sequence = region["sequence"]
            for start in range(len(sequence) - 19):
                stretch = sequence[start:start + 20]
                self.assertLess(sum(a in "AVLIMFWCY" for a in stretch), 17,
                                f"{key}: hydrophobic 20-mer {stretch}")

    def test_licence_and_attribution_travel_with_the_data(self):
        meta = available()
        self.assertEqual(meta["license"], "CC BY 4.0")
        self.assertIn("UniProt", meta["copyright"])


class AssemblyTests(unittest.TestCase):
    def test_trastuzumab_assembles_to_the_known_chain_lengths(self):
        """Heavy 450 aa with the C-terminal lysine, light 214 aa."""
        built = assemble(chains(), clip_c_terminal_lysine=False)
        self.assertEqual(built["heavy"]["full_length"], 450)
        self.assertEqual(built["light"]["full_length"], 214)
        self.assertEqual(built["total_residues"], 450 * 2 + 214 * 2)
        self.assertEqual(built["format"], "H2L2")

    def test_clipping_the_c_terminal_lysine_removes_one_residue_per_heavy_chain(self):
        clipped = assemble(chains(), clip_c_terminal_lysine=True)
        intact = assemble(chains(), clip_c_terminal_lysine=False)
        self.assertEqual(intact["total_residues"] - clipped["total_residues"], 2)
        self.assertTrue(clipped["c_terminal_lysine_clipped"])

    def test_light_constant_follows_the_called_chain_type(self):
        self.assertEqual(default_light_for("K"), "IGKC")
        self.assertEqual(default_light_for("L"), "IGLC2")
        self.assertEqual(assemble(chains("K"))["light"]["constant"], "IGKC")
        self.assertEqual(assemble(chains("L"))["light"]["constant"], "IGLC2")

    def test_heavy_chain_is_required_and_unknown_regions_are_refused(self):
        with self.assertRaises(ConstantRegionError):
            assemble({"L": chains()["L"]})
        with self.assertRaises(ConstantRegionError):
            assemble(chains(), heavy_isotype="IGHG9")
        with self.assertRaises(ConstantRegionError):
            assemble(chains(), heavy_isotype="IGKC")  # a light region asked for as heavy


class FullLengthFlagTests(unittest.TestCase):
    def test_full_length_changes_the_basis_and_clears_the_fv_only_marker(self):
        fv = molecule_flags(chains())
        full = molecule_flags(chains(), assembly=assemble(chains()))
        self.assertEqual(fv["basis"], "variable_domain")
        self.assertEqual(full["basis"], "full_length")
        self.assertTrue(all(f["fv_only"] for f in fv["flags"]))
        self.assertFalse(any(f["fv_only"] for f in full["flags"]))
        self.assertGreater(full["fv"]["residues"], 1000)

    def test_subclass_pi_ordering_matches_constant_region_chemistry(self):
        """IgG4's constant region is the most acidic of the four, IgG1 the least."""
        pis = {iso: molecule_flags(chains(), assembly=assemble(chains(), heavy_isotype=iso))["fv"]["pi"]
               for iso in ("IGHG1", "IGHG2", "IGHG3", "IGHG4")}
        self.assertLess(pis["IGHG4"], pis["IGHG2"])
        self.assertLess(pis["IGHG2"], pis["IGHG1"])

    def test_charge_bands_are_rescaled_so_they_do_not_fire_on_every_molecule(self):
        full = molecule_flags(chains(), assembly=assemble(chains()))
        charge = next(f for f in full["flags"] if f["id"] == "net_charge")
        self.assertGreater(abs(full["fv"]["net_charge_assay_ph"]), 5,
                           "full-length charge exceeds the Fv-scale amber band")
        self.assertEqual(charge["verdict"], "green",
                         "rescaled bands should not flag a well-behaved approved antibody")

    def test_cdr_hydrophobicity_is_unaffected_by_attaching_constant_regions(self):
        fv = molecule_flags(chains())
        full = molecule_flags(chains(), assembly=assemble(chains()))
        pick = lambda r: next(f for f in r["flags"] if f["id"] == "cdr_hydrophobicity")
        self.assertEqual(pick(fv)["value"], pick(full)["value"])


if __name__ == "__main__":
    unittest.main()
