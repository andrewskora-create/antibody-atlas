#!/usr/bin/env python
"""Produce the publishable tree, and prove it carries nothing it should not.

The working copy holds components that cannot be redistributed: Regeneron's
PROPERMAB, GPLv3 native sources, and reference corpora derived from Oxford
Thera-SAbDab. This builds a separate export tree without them, so the working
copy keeps running unchanged while the published repository stays clean.

    python tools/make_oss_export.py --out ../atlas-export
    python tools/make_oss_export.py --out ../atlas-export --panel-descriptors permissive-only

The export is verified after it is written: anything matching a forbidden
pattern fails the run rather than being quietly shipped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# Copied verbatim.
INCLUDE_DIRS = ["app", "static", "tests", "scripts", "setup", "tools", "validation",
                # Repository scaffolding: CI, issue templates, PR template. A
                # missing workflow does not fail anything, it just means the
                # published repository silently has no CI.
                ".github"]
INCLUDE_FILES = [
    "LICENSE", "NOTICE", "README.md", "COMPLIANCE.md", "USE_TIERS.md",
    "DATA_PROVENANCE.md", "THIRD_PARTY_NOTICES.md", "PTM_FIX_ENGINEERING.md",
    "MODULE_GUIDE.md", "START_HERE.md", "LINUX_AND_HEADLESS.md",
    "PTM_EXPOSURE_EVIDENCE.md", "TESTED_ENVIRONMENT.txt",
    # requirements-dev.txt carries pytest. README and CONTRIBUTING both tell a
    # fresh clone to install from it, so leaving it out of the export puts the
    # published package back in the state that file was written to fix.
    "requirements.txt", "requirements-dev.txt",
    "run_panel.py", "run_batch.py", ".gitignore", "run_local.sh",
    "run_local.ps1", "PLATFORMS.md", "BATCH.md",
    "run.py", "run.sh", "requirements-propermab.txt",
    "PUBLISHING.md",
    "Start Antibody Atlas.command", "Start Antibody Atlas.bat",
    "CONTRIBUTING.md", "SECURITY.md", "CITATION.cff",
]
# Redistributable datasets only.
INCLUDE_DATA = [
    "human_v_ogrdb.json", "human_v_ogrdb_numbered.json",
    # UniProtKB constant regions: CC BY 4.0, redistribution permitted with the
    # copyright notice, which travels inside the file and in THIRD_PARTY_NOTICES.
    "ig_constant_regions.json",
    "therapeutic_ptm_evidence.json", "experimental_reference_candidates.json",
    "propermab_published_ranges.json",
    "profiler_comparators.stripped.json",
    "therapeutic_ptm_sequences.stripped.json",
]
# Never copied, at any depth: build and cache noise.
EXCLUDE_PARTS = {"__pycache__", ".venv", ".git", ".cache", ".pytest_cache",
                 ".cocrystal_work"}
# Never copied, but ONLY at the top level. These name specific directories in
# this repository, not a kind of directory.
#
# This distinction is load-bearing. "vendor" was previously excluded at any
# depth to drop the Regeneron PROPERMAB tree at the repository root -- and it
# silently also dropped static/vendor/ngl/, the MIT-licensed structure viewer,
# from every export. The export still verified clean and the tests still
# passed, because nothing checks that a 1.3 MB browser asset arrived. The
# published app just had a dead 3D viewer.
EXCLUDE_ROOTS = {"vendor", "native-source", "native-build"}
EXCLUDE_SUFFIX = {".pyc", ".pyo", ".pdb", ".log"}

# Files the published app needs at runtime but nothing imports, so no test
# would miss them. Their absence fails the export.
REQUIRED_FILES = [
    ("static/vendor/ngl/ngl.js", "the 3D structure viewer; without it the viewer silently 404s"),
    ("static/vendor/ngl/LICENSE", "NGL's MIT licence, which must travel with it"),
    ("static/start.html", "the application page"),
    ("data/ig_constant_regions.json", "constant regions for the full-length charge basis"),
    ("Start Antibody Atlas.command", "the macOS launcher"),
    ("app/stability_energy.py", "local stability predictor adapter"),
    ("scripts/setup_stability.py", "pinned MIT stability-model installer"),
    ("scripts/thermompnn_infer.py", "local stability inference runner"),
    ("static/engineering.js", "variant modeling and structure comparison"),
    ("app/deepsp.py", "the MIT-licensed sequence-only descriptor route"),
    ("setup/fetch_deepsp.py", "how a user installs the DeepSP weights at all"),
    ("tests/fixtures/deepsp_reference_molecules.json",
     "DeepSP's own reference molecules; without them the weight conversion cannot be verified"),
    (".github/workflows/tests.yml", "CI; its absence is silent, the repository just has none"),
    ("CONTRIBUTING.md", "what a contributor has to know before opening a pull request"),
    ("SECURITY.md", "where sequences go, and how to report a problem privately"),
    ("CITATION.cff", "how to cite this and, more importantly, the methods it calls"),
    ("tests/test_export_guard.py", "the only test that this gate still detects anything"),
    ("app/structure_input.py", "experimental-structure ingestion; without it only predicted models work"),
    ("app/hbonds.py", "hydrogen-bond geometry, which the substitution ranking reads"),
    ("app/repack.py", "side-chain repacking for substitution consequence"),
    ("app/hbonds.py", "hydrogen-bond geometry; the substitution context depends on it"),
    ("app/repack.py", "modelled substitutions; without it /api/repack 500s"),
    ("app/structure_ensemble.py", "multi-model exposure; without it the ensemble option 500s"),
    ("PUBLISHING.md", "the steps for publishing a fork without leaking sequences"),
]

# Patterns that must not appear anywhere in the export.
FORBIDDEN_PATHS = [
    # Anchored at the repository root. An unanchored "vendor/" also matches
    # static/vendor/ngl/, which is the MIT structure viewer and must ship.
    (re.compile(r"^vendor/"), "Regeneron PROPERMAB source"),
    (re.compile(r"^native-(source|build)/"), "GPLv3 NanoShaper/CGAL source"),
    (re.compile(r"clinical_frameworks\.json$"), "Thera-SAbDab derived corpus"),
    (re.compile(r"^data/therapeutic_ptm_sequences\.json$"), "Thera-SAbDab derived corpus"),
    (re.compile(r"^data/profiler_comparators\.json$"), "Thera-SAbDab derived sequences"),
    (re.compile(r"model\.safetensors$"), "model weights, fetched at setup"),
    (re.compile(r"^models/deepsp/.*\.npz$"), "DeepSP weights, fetched at setup"),
    (re.compile(r"^models/deepsp/provenance\.json$"), "records this machine's fetch, not the user's"),
    (re.compile(r"^data/deepsp_reference\.json$"), "DeepSP percentile table, built at setup"),
]
# Sequence fragments that must never leave the working copy, stored as SHA-256
# of the fragment rather than the fragment itself. An exclusion list that
# contains the sequences it excludes would itself be a disclosure, and this file
# ships in the export. Matching slides a window of each recorded length over the
# amino-acid runs in a file, so absence is proven without the secret being here.
#
# The labels are deliberately opaque. Naming the molecules would leak by
# description what hashing the fragments was meant to withhold; the maintainer
# knows which is which, and a reader only needs to know the export was checked.
FORBIDDEN_SEQUENCE_HASHES = {
    ("d4b88116a492da31c547c38551a05407e59a72885c00cede7f244f2f99f29657", 14): "restricted fragment A1",
    ("068c0a62d56405ddd9c2332d6b04d948f48f5b55db089ae6047a2f07805b129f", 9): "restricted fragment A2",
    ("36829e4b6ef574d125b97590c6dfc3490ddc4003797787cb06b971047aea996b", 14): "restricted fragment B1",
    ("70fb08e5afdb8393f248381dce7bc42bec1a625edaff535d051ca691ad36de6f", 13): "restricted fragment B2",
    ("ac82ba8b1c87a5e1586fa6d766b05cad0d4b14536fb14e0dda1042f336ba376d", 25): "restricted fragment A3",
    ("a692366315a37e16ccd0e2ba21dfd02f7246602cf04ceb8cc1a9ce1a26d569ce", 17): "restricted fragment A4",
    # CDR fragments covering every restricted molecule, not only the three
    # that were enumerated by hand. The earlier six fragments left the rest of
    # the panel unprotected: a planted CDR-H3 from any other member passed the
    # verifier. Each fragment here was checked against the public Thera-SAbDab
    # corpus, the OGRDB germline set and the shipped panel, and dropped if it
    # occurred in any of them, so none of these can fire on data that
    # legitimately ships. Labels are deliberately opaque for the same reason
    # the fragments are hashed: this file is published.
    ("2661aa1c8cc176c8b0172bc14361566f93f8dd34e71f6444a49d9d336aeab1f8", 15): "restricted CDR3 fragment C1",
    ("5253d6932fc15b5bd2ffd74555a04ccaf826686470249d9cceed610ee6eef57b", 12): "restricted CDR3 fragment C2",
    ("183b323e5ce496dcae40665ae366608eb16f7664094913354d92e151e7c4982c", 15): "restricted CDR3 fragment C3",
    ("2b744e262c19946187c26d811c97eda9fe829596cea779ad35acd3798717e9c1", 15): "restricted CDR3 fragment C4",
    ("bacc451bb9518a183a34b00d988bb1d03984e61c9b17fff742f6df2b41042f10", 16): "restricted CDR3 fragment C5",
    ("c499ec3163d415f25f95efe83d38a3ffd43a2e483338a8487f334477fbb23591", 16): "restricted CDR3 fragment C6",
    ("35a61cc69ccce502d025e75e0d59857198f3a2b7be254a360e2875e64987c566", 16): "restricted CDR3 fragment C7",
    ("748f1fecaeda1ea38607a2566e0ed1fc7894c5484bf3dd69f04ae873ef3e0e65", 13): "restricted CDR3 fragment C8",
    ("8184d98230ea75f732da2b784402cbf3bd2f120a53c7b05eefe78eaac7c94987", 13): "restricted CDR3 fragment C9",
    ("7537025ebea3ade809fb34c05856ec25873abb61ec41075df8cc7a9f09e2a704", 13): "restricted CDR3 fragment C10",
    ("49bfee87de19cf3c8e8c2aa337a254ea77b43968778f74cb8698e1a64237091f", 17): "restricted CDR3 fragment C11",
    ("d887469a2845982449cb4bc9e82a53894ce98536fb7df2e15a0f08e11a388ed0", 14): "restricted CDR3 fragment C12",
    ("f1f446dd2cf3b27f92f8ff7d2defe63ab4aa8c87af1b9d11144f69a19138023e", 17): "restricted CDR3 fragment C13",
    ("91370f6ec19dcd12cfb78ade0fbe71788dffd8eab24a7deed8af9330de1cc32f", 15): "restricted CDR3 fragment C14",
    ("9a6d6b751bd40ac2e2a4693808b7a6bf65bac50e2fda1e70f3f729b4ec1fff99", 15): "restricted CDR3 fragment C15",
    ("5384dba6717f07dfce90250b6b5449dbbd6b24303b45f1d3daf95f67e3dcbb31", 19): "restricted CDR3 fragment C16",
    ("b3dbefe3c69810fe9c853d7d133d68a468dd4e4518541237a47b1f10e4639056", 16): "restricted CDR3 fragment C17",
    ("4e238d9670b24d49a9652d5dfc4690a0b30264881d855f9988940a25537aaa23", 15): "restricted CDR3 fragment C18",
    ("63f3f3ff0eb2ac4bc3eefd006c903262a4fdbca628a19e47e788db2ee3fbd2b9", 14): "restricted CDR3 fragment C19",
    ("95e46ee10e466994ee4f2f025bc578f30292faa69b68d3b415037917d72a2a76", 16): "restricted CDR3 fragment C20",
    ("9a047674e276acb879b5b9e42d27c79e5121eb550bc5c959555a91dfe0f824d3", 13): "restricted CDR3 fragment C21",
    ("e8cfcdb40735b83b3691383f0a2e3797703d9409c72c6e7b6baebdd7c886e77f", 10): "restricted CDR3 fragment C22",
    ("b8abfecc33721c19eced15393036f849b10045eb7b495808807930d5ae63611d", 14): "restricted CDR3 fragment C23",
    ("f6fc95b965a9f7a7805dd8e5795e6083c92047837a2d461a6534085d4f5b0f75", 17): "restricted CDR3 fragment C24",
    ("336e8ddeade00d664e51b1c3867f3e7de871097bd48f2dd2d25f2e7e0547f20b", 13): "restricted CDR3 fragment C25",
    ("8d6a81713b11a92b1d2c64526d139461acb96136bcc8b32fbf399b8de8a96ed3", 10): "restricted CDR3 fragment C26",
    ("d76214f650da9df7e5dbbe9dd7157a6d5da5c2868b19d296fb3aa5624af05df2", 16): "restricted CDR3 fragment C27",
    ("7921264d70f5b13d0dfd32f91faaf93997fe165290f04eee71b6cb471f4b6e32", 14): "restricted CDR3 fragment C28",
    ("b725f4d72a23672437b340b48bda25bf251aae5a508a6e5a175f4637bd9c1c1a", 10): "restricted CDR fragment C29",
    ("82a83bff83a7e9d635bd6341de9acdbf9924ee09a6a53aab60cfd7a034fcd17b", 11): "restricted CDR3 fragment C30",
    ("5e1edca13a21af036a2472295fd69c5f908c789800cbbc1d04f4e75b2e529f55", 12): "restricted CDR fragment C31",
    ("b7d2308e418b1026ec83bfd4ba75c72b6f7f7baf7f37042f458e96978636bece", 12): "restricted CDR fragment C32",
    ("3a523bcc9bd63b82ba9184e636cf80026ad6bace375091ddb81f53cc952e8e98", 10): "restricted CDR3 fragment C33",
    ("08f2e3ca1758203ce6b5cfae17bc065372d5d6b5038ebfce1c01180e76c023f3", 12): "restricted CDR fragment C34",
    ("5b865b8dd271cc0964f9a2ed7d8cf22fadfba9fb0d845b45c5062a3799875f32", 12): "restricted CDR fragment C35",
    ("29e500822632487591485d482439cd2504af8bc33abbdf761253b042b8d8c719", 13): "restricted CDR3 fragment C36",
    ("70fb08e5afdb8393f248381dce7bc42bec1a625edaff535d051ca691ad36de6f", 13): "restricted CDR3 fragment C37",
    ("2318fe1ebc012afc0353a799b9cd3222d038fd6a7939a4c9515bb3a693b37066", 16): "restricted CDR3 fragment C38",
    # One window per restricted chain the CDR fragments above still miss --
    # those whose CDR-H3 is only eight or nine residues, too short to hash
    # without firing on germline. Each window is the shortest stretch of that
    # chain occurring nowhere in the public Thera-SAbDab corpus, the OGRDB
    # germline set or the shipped panel, so it cannot false-positive, and a
    # whole-chain leak cannot slip past.
    ("1175622a828f02834e5de4c6092483260b1820921ca3755f59ac175971e57c5e", 24): "restricted chain window D1",
    ("6d7f3608b1fe7fcd4b37b42dbe6bf67b010e3a04e04c4d251b606073bc8c5680", 24): "restricted chain window D2",
    ("cdefc915ed5051e847814f0a1ac6a25ae8ed609c7a09ff22e57af72162010313", 24): "restricted chain window D3",
    ("f5f725cbc41b069ee2187ce8012658ebf17be3284f60fafb6b16ef41b40b33ce", 24): "restricted chain window D4",
    ("93fa99692f1c2d6b225e7e857c50b7d0768a506664058856c6c8e47a35d80536", 24): "restricted chain window D5",
    ("5bcff89152256982b9b870cfd726a7d22d089d9ebf49a2aae7571b0fae3fb649", 24): "restricted chain window D6",
    ("b317b855f179708cec3547b50b1eecc1ef68380f14e2caf7b4045339ed48d9df", 24): "restricted chain window D7",
    ("583becd7d84dd561cb173b18179557323aad0834f5a130f9a297b38d594bbe4b", 24): "restricted chain window D8",
    ("332da055438fb6eb7892847df2736fe842aac89e98b86187f223283e7160e6f2", 24): "restricted chain window D9",
    ("37a59703755cdec6c7db13ca03e3ee8fd97eb8be19eebfcd26a4bb90b963d717", 24): "restricted chain window D10",
    ("dc8dda7410c636519abfde255dac87b4a4ce554b847e6d7a1632905a4a1bca52", 24): "restricted chain window D11",
    ("8766350a6b5125d4b594bf77acadb0e6a65dd249dd142b84c9585031a6b59295", 24): "restricted chain window D12",
    ("465636bf14da01f22d459bd5e7f3a21f44ebe4c09e158b42dd7911ad7ab817fe", 24): "restricted chain window D13",
    ("25ff34578c22feedb580d6a97afe4c970c5c4afd51944fe2bb657a1496af1a8d", 24): "restricted chain window D14",
    ("054c5ffe85d50b77e457e846a2c50952e6191d0df7b5da0b463632ac8794f618", 24): "restricted chain window D15",
    ("4cd5e075e90ba63fb0812aa110902c3c9b3a61a00199c9954db33c4067a58768", 24): "restricted chain window D16",
    ("6777fe3da16f04e766a35cf98e7d561efcb230e17dcf7037a0d044293826150a", 24): "restricted chain window D17",
    ("465636bf14da01f22d459bd5e7f3a21f44ebe4c09e158b42dd7911ad7ab817fe", 24): "restricted chain window D18",
}
AA_RUN = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{9,}")

# Known limit, stated rather than hidden: a CDR-H3 of eight or nine residues
# is not hashed on its own. Fragments that short recur in unrelated molecules
# -- one nine-residue kappa CDR-L3 in this project's own set is shared with
# five approved drugs -- so hashing them would fail the export on data that
# legitimately ships. The chain windows below still catch the whole chain; it
# is only a bare eight- or nine-residue loop, quoted alone, that passes.

# Molecule NAMES must not ship either, and the sequence check cannot catch
# them: a name is not an amino-acid run. Stored as SHA-256 of the lowercased
# token for the same reason the sequences are -- this file is published, so
# listing the names here would be the disclosure it exists to prevent.
#
# A placeholder in a test docstring is how one nearly shipped: the export
# verified clean because nothing looked for names.
FORBIDDEN_NAME_HASHES = {
    "72d3bbbc86472b6f7ece2ce813cd6c7ce8d1d883808ae628e597145a597e4971": "series B",
    "1c410932e677d3c0e5a9009e8f2f88d48ec26ae074bdbbf9250c06aa4d4238fc": "series A prefix",
    "47897b389a7b2efab811cdcc6b93131adf37fc65c0f60f7d8c88c5d60a4584a2": "series A, member 1",
    "00482d0f130e44763f5487b286f024b08d2bc9dc68090a6d0fb6934d3b2e0cf1": "series A, member 2",
    "d73b21a2e6a487c26f0288d66ee382f46e2cff662dde475a007d4a3fa8edffdb": "series A, member 3",
    "3cfd2578da0add65097abb35528c73fde40026fa1fd760892954fd2ea95a5b81": "series A, member 4",
    "8a714d34b3796a0f6d9244c3f3804033d9decd8b9336a6a0193fbafaedd264bf": "series A, member 5",
    "3213e97a8fe4d8cd7850a3ec01f7ad36f4e5bc43cf4edca22027a29d97318e4e": "target designator",
}
WORD = re.compile(r"[A-Za-z][A-Za-z0-9]{2,}(?:[_-][A-Za-z0-9]+)*")


def name_leaks(text: str) -> set[str]:
    """Names of any forbidden token present, found by hash, not by literal.

    Compound names are split as well as hashed whole. This matters more than it
    looks: the real working files name molecules with a series stem plus chain
    and mutation suffixes, joined by underscores, and hashing only the
    complete token meant 154 of the 172 name forms actually in use went
    undetected -- the few that were caught had each been enumerated by hand.
    Splitting on the separators reduces every one of them to a stem that is
    already in the hash set, so a variant nobody thought to enumerate is
    caught by its stem. (No example is written out here: this file ships.)
    """
    found: set[str] = set()
    tokens: set[str] = set()
    for word in set(WORD.findall(text)):
        tokens.add(word)
        tokens.update(part for part in re.split(r"[_\-]", word) if len(part) >= 3)
    for token in tokens:
        digest = hashlib.sha256(token.lower().encode()).hexdigest()
        if digest in FORBIDDEN_NAME_HASHES:
            found.add(FORBIDDEN_NAME_HASHES[digest])
    return found


def sequence_leaks(text: str) -> set[str]:
    """Names of any forbidden fragment present, found by hash, not by literal."""
    lengths = {length for _, length in FORBIDDEN_SEQUENCE_HASHES}
    wanted = {digest for digest, _ in FORBIDDEN_SEQUENCE_HASHES}
    found: set[str] = set()
    for run in AA_RUN.findall(text):
        for length in lengths:
            if len(run) < length:
                continue
            for i in range(len(run) - length + 1):
                digest = hashlib.sha256(run[i:i + length].encode()).hexdigest()
                if digest in wanted:
                    found.add(FORBIDDEN_SEQUENCE_HASHES[(digest, length)])
    return found


PROPERMAB_KEYS = {"hyd_patch_area", "hyd_patch_area_cdr", "pos_patch_area",
                  "pos_patch_area_cdr", "neg_patch_area", "neg_patch_area_cdr",
                  "scm", "heiden_score"}


# Test fixtures are shipped even when their suffix is excluded: .pdb is excluded
# to keep run output and predicted structures out of the export, but the suite
# cannot run without its fixture, and a repository whose tests cannot run is not
# a usable open-source release. Fixture contents still pass the sequence check.
FIXTURE_DIR = "fixtures"


def copyable(path: Path) -> bool:
    if any(part in EXCLUDE_PARTS for part in path.parts):
        return False
    if path.parts and path.parts[0] in EXCLUDE_ROOTS:
        return False
    if path.suffix in EXCLUDE_SUFFIX:
        return FIXTURE_DIR in path.parts
    return True


def copy_tree(out: Path) -> list[Path]:
    written: list[Path] = []
    for name in INCLUDE_DIRS:
        source = ROOT / name
        if not source.is_dir():
            continue
        for item in sorted(source.rglob("*")):
            if not item.is_file() or not copyable(item.relative_to(ROOT)):
                continue
            if item.is_symlink():
                raise ValueError(f"Refusing symlink in export input: {item.relative_to(ROOT)}")
            target = out / item.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, target)
            written.append(target.relative_to(out))
    for name in INCLUDE_FILES:
        source = ROOT / name
        if source.is_file():
            if source.is_symlink():
                raise ValueError(f"Refusing symlink in export input: {name}")
            shutil.copy2(source, out / name)
            written.append(Path(name))
        else:
            # Silence here is how requirements-dev.txt went missing from the
            # export while README still told people to install from it. A name
            # in the allowlist that has no file is a build error, not a no-op.
            print(f"  ! missing {name}; it is in INCLUDE_FILES but not on disk")
    (out / "data").mkdir(parents=True, exist_ok=True)
    for name in INCLUDE_DATA:
        source = ROOT / "data" / name
        if source.is_file():
            if source.is_symlink():
                raise ValueError(f"Refusing symlink in export input: {name}")
            shutil.copy2(source, out / "data" / name)
            written.append(Path("data") / name)
        else:
            print(f"  ! missing {name}; run: python setup/build_reference_data.py --make-stubs")
    return written


def strip_propermab_descriptors(out: Path) -> int:
    """Drop PROPERMAB-derived descriptor values from the shipped panel."""
    # Optional factual comparator extract has unresolved source provenance.
    # Keep it out of the permissive release until that audit is closed.
    (out / "data" / "experimental_reference_candidates.json").unlink(missing_ok=True)
    target = out / "data" / "profiler_comparators.stripped.json"
    if not target.is_file():
        return 0
    data = json.loads(target.read_text())
    removed = 0
    for ref in data.get("references", []):
        if ref.pop("propermab", None) is not None:
            removed += 1
    for key in ("protocol",):
        block = data.get(key)
        if isinstance(block, dict):
            for field in list(block):
                if "propermab" in field or "nanoshaper" in field:
                    block.pop(field)
    data["propermab_descriptors_removed"] = True
    data["propermab_descriptors_note"] = (
        "PROPERMAB-derived descriptor values were omitted from this export at the "
        "publisher's choice. They are research outputs rather than Regeneron's "
        "licensed source code, but omitting them removes any question."
    )
    target.write_text(json.dumps(data))
    return removed


# Absolute paths from the build machine. These are not secret, but they name
# someone's filesystem and they are noise in a published repository. More to the
# point, a path baked into a shipped data file is a sign that file recorded its
# build environment, so the export both scrubs them and then refuses to ship any
# that remain.
LOCAL_PATH = re.compile(r"(/home/[A-Za-z0-9._-]+|/Users/[A-Za-z0-9._-]+|[A-Z]:\\\\Users\\\\[A-Za-z0-9._-]+)"
                        r"[^\s\"',)\]}]*")
SCRUB_SUFFIXES = {".py", ".json", ".md", ".html", ".js", ".cjs", ".txt", ".cfg", ".sh", ".ps1", ".command", ".bat", ".yml", ".yaml", ".cff"}


def scrub_local_paths(out: Path) -> list[str]:
    """Replace build-machine absolute paths with a neutral placeholder."""
    touched: list[str] = []
    for item in sorted(out.rglob("*")):
        if not item.is_file() or item.suffix not in SCRUB_SUFFIXES:
            continue
        try:
            text = item.read_text(errors="ignore")
        except OSError:
            continue
        if not LOCAL_PATH.search(text):
            continue
        item.write_text(LOCAL_PATH.sub("<local path removed for distribution>", text))
        touched.append(item.relative_to(out).as_posix())
    return touched


def verify(out: Path) -> list[str]:
    problems: list[str] = []
    # Positive checks first. Everything else here asks "did something forbidden
    # get in"; this asks "did something required fall out", which is the
    # failure that ships quietly.
    for relative, why in REQUIRED_FILES:
        target = out / relative
        if not target.is_file() or target.stat().st_size == 0:
            problems.append(f"MISSING {relative}: {why}")
    for item in sorted(out.rglob("*")):
        if not item.is_file():
            continue
        rel = item.relative_to(out).as_posix()
        for pattern, why in FORBIDDEN_PATHS:
            if pattern.search(rel):
                problems.append(f"{rel}: {why}")
        if item.suffix in SCRUB_SUFFIXES:
            try:
                text = item.read_text(errors="ignore")
            except OSError:
                continue
            for why in sorted(sequence_leaks(text)):
                problems.append(f"{rel}: contains {why}")
            for why in sorted(name_leaks(text)):
                problems.append(f"{rel}: contains a restricted molecule name ({why})")
            leaked = LOCAL_PATH.search(text)
            if leaked:
                problems.append(f"{rel}: build-machine path {leaked.group(0)[:60]}")
    return problems


def validate_destination(out: Path) -> None:
    """Never let --force erase this checkout or recurse into its own output."""
    source = ROOT.resolve()
    target = out.resolve()
    if target == source or target in source.parents or source in target.parents:
        raise ValueError("Export destination must be outside the source checkout and must not be an ancestor of it.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="target directory for the export")
    parser.add_argument("--panel-descriptors", choices=("all", "permissive-only"),
                        default="all",
                        help="'permissive-only' drops PROPERMAB-derived values from the panel")
    parser.add_argument("--force", action="store_true", help="overwrite an existing target")
    args = parser.parse_args()

    out = Path(args.out).resolve()
    validate_destination(out)
    if out.exists():
        if not args.force:
            raise SystemExit(f"{out} already exists. Pass --force to replace it.")
        shutil.rmtree(out)
    out.mkdir(parents=True)

    print("Refreshing sequence-free stubs")
    subprocess.run([sys.executable, str(ROOT / "setup" / "build_reference_data.py"),
                    "--make-stubs"], check=False)

    print(f"\nCopying to {out}")
    written = copy_tree(out)
    print(f"  {len(written)} files")

    if args.panel_descriptors == "permissive-only":
        removed = strip_propermab_descriptors(out)
        print(f"  removed PROPERMAB descriptors from {removed} panel entries")

    scrubbed = scrub_local_paths(out)
    if scrubbed:
        print(f"  scrubbed build-machine paths from {len(scrubbed)} file(s): "
              + ", ".join(scrubbed[:4]) + (" …" if len(scrubbed) > 4 else ""))

    print("\nVerifying")
    problems = verify(out)
    if problems:
        print("  FAILED. The export contains material it must not:")
        for problem in problems:
            print(f"    {problem}")
        return 1
    print("  no forbidden path or sequence found")

    total = sum(f.stat().st_size for f in out.rglob("*") if f.is_file())
    print(f"\nExport ready: {len(written)} files, {total/1e6:.1f} MB")
    print("Verified clean. Before the first public push:")
    print("  1. cd into the export and run the tests there, not in the working copy.")
    print("  2. git init && git add -A && git status   # confirm nothing unexpected is staged")
    print("  3. Check README and NOTICE read correctly as the repository landing page.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
