"""The export guard has to actually detect things.

This exists because the guard had no test and quietly rotted. Its name check
hashed only whole tokens, so the compound forms the working files really use
-- a series stem plus chain and mutation suffixes joined by underscores --
went undetected: 154 of 172 real name forms passed. Its sequence check carried
six fragments covering three molecules, so a CDR-H3 from any other member of
the panel passed too. The published tree was clean the whole time, which is
exactly why nobody noticed.

No restricted sequence or name appears here. The tests install synthetic
secrets into the hash tables and check the MECHANISM, which is what broke.
The real hash tables are checked only for shape and size.
"""
from __future__ import annotations

import hashlib
import re
import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("mox", ROOT / "tools" / "make_oss_export.py")
mox = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mox)


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


class NameDetection(unittest.TestCase):
    """A stem in the table must catch every compound built on it."""

    def setUp(self):
        self.saved = dict(mox.FORBIDDEN_NAME_HASHES)
        mox.FORBIDDEN_NAME_HASHES.clear()
        mox.FORBIDDEN_NAME_HASHES[sha("zyxwv")] = "test stem"

    def tearDown(self):
        mox.FORBIDDEN_NAME_HASHES.clear()
        mox.FORBIDDEN_NAME_HASHES.update(self.saved)

    def test_bare_stem(self):
        self.assertEqual(mox.name_leaks("zyxwv run"), {"test stem"})

    def test_case_insensitive(self):
        self.assertEqual(mox.name_leaks("ZyXwV"), {"test stem"})

    def test_underscore_compound(self):
        """The regression: the whole token is not in the table, the stem is."""
        self.assertEqual(mox.name_leaks("ZYXWV_VH_N55A"), {"test stem"})

    def test_hyphen_compound(self):
        self.assertEqual(mox.name_leaks("zyxwv-treated samples"), {"test stem"})

    def test_deep_compound(self):
        self.assertEqual(mox.name_leaks("ZYXWV_CI3_D52E_HM102L"), {"test stem"})

    def test_stem_as_a_later_part(self):
        self.assertEqual(mox.name_leaks("panel_zyxwv_26"), {"test stem"})

    def test_unrelated_text_is_not_flagged(self):
        self.assertEqual(mox.name_leaks("bevacizumab trastuzumab IGHV3-23"), set())

    def test_substring_of_a_longer_word_is_not_flagged(self):
        """Stem matching must be on token parts, not raw substrings, or every
        word containing the stem would fail the export."""
        self.assertEqual(mox.name_leaks("zyxwvwxyz something"), set())


class SequenceDetection(unittest.TestCase):
    def setUp(self):
        self.saved = dict(mox.FORBIDDEN_SEQUENCE_HASHES)
        mox.FORBIDDEN_SEQUENCE_HASHES.clear()
        self.frag = "WWQQMMKKHHRRYY"
        mox.FORBIDDEN_SEQUENCE_HASHES[(sha(self.frag), len(self.frag))] = "test fragment"

    def tearDown(self):
        mox.FORBIDDEN_SEQUENCE_HASHES.clear()
        mox.FORBIDDEN_SEQUENCE_HASHES.update(self.saved)

    def test_exact_fragment(self):
        self.assertEqual(mox.sequence_leaks(f'S = "{self.frag}"'), {"test fragment"})

    def test_fragment_embedded_in_a_longer_chain(self):
        chain = "EVQLVESGGGLVQPGGSLRLSCAAS" + self.frag + "WGQGTLVTVSS"
        self.assertEqual(mox.sequence_leaks(chain), {"test fragment"})

    def test_fragment_split_by_a_non_residue_is_not_matched(self):
        """A run is amino acids only; a fragment broken by punctuation is a
        different string and must not be reported as present."""
        broken = self.frag[:7] + "-" + self.frag[7:]
        self.assertEqual(mox.sequence_leaks(broken), set())

    def test_unrelated_sequence_is_not_flagged(self):
        self.assertEqual(mox.sequence_leaks("EVQLVESGGGLVQPGGSLRLSCAASGFTFS"), set())


class RealTablesShape(unittest.TestCase):
    """Shape only. The values are hashes precisely so this file can ship."""

    def test_sequence_table_covers_more_than_a_handful(self):
        self.assertGreaterEqual(
            len(mox.FORBIDDEN_SEQUENCE_HASHES), 40,
            "the table shrank; it is meant to cover every restricted molecule, "
            "not the few that were enumerated by hand")

    def test_no_table_entry_is_a_plaintext_sequence(self):
        """An exclusion list containing what it excludes is the disclosure."""
        for digest, length in mox.FORBIDDEN_SEQUENCE_HASHES:
            self.assertRegex(digest, r"^[0-9a-f]{64}$")
            self.assertIsInstance(length, int)
        for digest in mox.FORBIDDEN_NAME_HASHES:
            self.assertRegex(digest, r"^[0-9a-f]{64}$")

    def test_labels_do_not_name_a_molecule(self):
        """Labels are opaque on purpose: naming the molecule would leak by
        description what hashing the fragment withholds.

        Checked by running the name detector over the labels themselves rather
        than by matching literal stems -- spelling a stem out here would put
        the very token in the published tree that this gate exists to keep
        out, and would fail the export.
        """
        for label in list(mox.FORBIDDEN_SEQUENCE_HASHES.values()) + list(
                mox.FORBIDDEN_NAME_HASHES.values()):
            self.assertEqual(mox.name_leaks(label), set(),
                             f"label {label!r} names a restricted molecule")

    def test_short_fragments_are_not_hashed(self):
        """Below ten residues a fragment recurs in unrelated molecules and
        would fail the export on data that legitimately ships."""
        too_short = [n for _, n in mox.FORBIDDEN_SEQUENCE_HASHES if n < 9]
        self.assertEqual(too_short, [])


class DocumentedFilesAreExported(unittest.TestCase):
    """A file the docs tell people to use has to be in the export.

    requirements-dev.txt was written to fix "No module named pytest" on a fresh
    clone, and README and CONTRIBUTING were pointed at it -- but it was never
    added to INCLUDE_FILES, so the published package shipped the instruction
    without the file and reproduced the exact failure the file was created to
    remove. The docs and the allowlist are two hand-kept lists; this makes them
    disagree loudly instead of silently.
    """

    DOCS = ("README.md", "CONTRIBUTING.md")

    def test_requirements_files_named_in_docs_are_included(self):
        named = set()
        for doc in self.DOCS:
            path = ROOT / doc
            if not path.is_file():
                continue
            named |= set(re.findall(r"(requirements[\w.-]*\.txt)", path.read_text()))
        self.assertTrue(named, "no requirements file is named in the docs")
        missing = sorted(n for n in named if n not in mox.INCLUDE_FILES)
        self.assertEqual(
            missing, [],
            f"named in {self.DOCS} but absent from INCLUDE_FILES: {missing}")

    def test_every_include_file_exists_on_disk(self):
        """An allowlist entry with no file copies nothing and said nothing."""
        absent = sorted(n for n in mox.INCLUDE_FILES if not (ROOT / n).is_file())
        self.assertEqual(absent, [],
                         f"in INCLUDE_FILES but not on disk: {absent}")


if __name__ == "__main__":
    unittest.main()


class ReleaseBoundaryTests(unittest.TestCase):
    def test_destination_cannot_erase_or_recurse_into_source(self):
        for target in (mox.ROOT, mox.ROOT.parent, mox.ROOT / 'export'):
            with self.subTest(target=target), self.assertRaises(ValueError):
                mox.validate_destination(target)
        mox.validate_destination(mox.ROOT.parent / 'separate-export')

    def test_script_and_workflow_sources_are_scanned(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for suffix in ('.cjs', '.yml', '.command', '.ps1'):
                (root / ('sample' + suffix)).write_text('/Users/' + 'test-user/private.csv')
            problems = mox.verify(root)
            for suffix in ('.cjs', '.yml', '.command', '.ps1'):
                self.assertTrue(any('sample' + suffix + ': build-machine path' in p for p in problems))


class PermissiveReleaseTests(unittest.TestCase):
    def test_unresolved_optional_reference_is_not_distributed(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp); (out/'data').mkdir()
            target=out/'data/experimental_reference_candidates.json'
            target.write_text('{"rows":[]}')
            mox.strip_propermab_descriptors(out)
            self.assertFalse(target.exists())
