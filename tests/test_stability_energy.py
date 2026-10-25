"""Scientific boundary and residue-mapping tests for local ΔΔG prediction."""
import hashlib
from pathlib import Path

import pytest

from app import stability_energy as energy


def pdb_text(residues=(("12", "SER"), ("12A", "THR"), ("15", "ASN"))):
    lines = []
    serial = 1
    for i, (rid, aa) in enumerate(residues):
        number = int(rid.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
        insertion = rid[-1] if rid[-1].isalpha() else " "
        for atom, dx in (("N", 0.0), ("CA", 1.1), ("C", 2.0), ("O", 2.5)):
            lines.append(f"ATOM  {serial:5d} {atom:^4s} {aa} H{number:4d}{insertion}   {i*3.3+dx:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{90.:6.2f}          {atom[0]:>2s}")
            serial += 1
    return "\n".join(lines) + "\nEND\n"


@pytest.fixture
def parent(tmp_path):
    path = tmp_path / "parent.pdb"
    path.write_text(pdb_text())
    return path


@pytest.fixture(autouse=True)
def clean_cache():
    energy._CACHE.clear()
    yield
    energy._CACHE.clear()


def mock_model(monkeypatch):
    calls = []
    monkeypatch.setattr(energy, "_verified_manifest", lambda: {"checkpoint_sha256": "verified"})
    def run(parsed):
        calls.append(parsed)
        return {"sequence": parsed["sequence"], "predictions": [{aa: (-1.25 if aa == "I" else 0.5) for aa in energy.AA} for _ in parsed["sequence"]]}
    monkeypatch.setattr(energy, "_run_model", run)
    return calls


def test_insertion_and_numbering_gaps_preserve_identity(parent):
    parsed = energy._parse_chain(parent.read_text(), "H")
    assert parsed["sequence"] == "STN"
    assert parsed["residue_ids"] == ["12", "12A", "15"]


def test_prediction_uses_exact_residue_and_caches_chain(parent, monkeypatch):
    calls = mock_model(monkeypatch)
    first = energy.predict_stability(parent, "H", "12A", "T", "I")
    second = energy.predict_stability(parent, "H", "15", "N", "A")
    assert first["status"] == second["status"] == "available"
    assert first["value_kcal_mol"] == -1.25
    assert len(calls) == 1
    assert first["provenance"]["pdb_residue"] == "12A"
    assert first["provenance"]["input_sha256"] == hashlib.sha256(parent.read_bytes()).hexdigest()
    assert "negative predicts stabilization" in first["sign_convention"]


def test_mismatched_parent_never_predicts(parent, monkeypatch):
    calls = mock_model(monkeypatch)
    result = energy.predict_stability(parent, "H", "12A", "S", "A")
    assert result["status"] == "unavailable"
    assert "does not match" in result["reason"]
    assert not calls


def test_changed_structure_invalidates_cache(parent, monkeypatch):
    calls = mock_model(monkeypatch)
    energy.predict_stability(parent, "H", "12", "S", "A")
    parent.write_text(parent.read_text() + "REMARK changed parent input\n")
    energy.predict_stability(parent, "H", "12", "S", "A")
    assert len(calls) == 2


@pytest.mark.parametrize("chain,rid,wt,mut", [("Z","12","S","A"), ("H","12B","S","A"), ("H","12","","A"), ("H","12","S","AA"), ("H","12.0","S","A"), ("HH","12","S","A")])
def test_invalid_requests_fail_closed(parent, chain, rid, wt, mut):
    result = energy.predict_stability(parent, chain, rid, wt, mut)
    assert result["status"] == "unavailable"
    assert result["value_kcal_mol"] is None


def test_missing_backbone_is_not_silently_skipped(parent):
    parent.write_text("\n".join(line for line in parent.read_text().splitlines() if not (line[12:16].strip() == "O" and line[26:27] == "A")))
    result = energy.predict_stability(parent, "H", "12", "S", "A")
    assert result["status"] == "unavailable"
    assert "Incomplete backbone" in result["reason"]


def test_noncontiguous_chain_rejected(parent):
    text = parent.read_text().replace("   3.300", "  30.300")
    parent.write_text(text)
    result = energy.predict_stability(parent, "H", "12", "S", "A")
    assert "Discontinuous backbone" in result["reason"]


def test_missing_engine_does_not_break_repack(parent, monkeypatch):
    monkeypatch.setenv("ATLAS_THERMOMPNN_DIR", str(parent.parent / "missing"))
    assert energy.stability_capability()["available"] is False
    result = energy.predict_stability(parent, "H", "12", "S", "A")
    assert result["status"] == "unavailable"
    assert result["value_kcal_mol"] is None


def test_inference_failure_returns_unavailable(parent, monkeypatch):
    mock_model(monkeypatch)
    monkeypatch.setattr(energy, "_run_model", lambda _: (_ for _ in ()).throw(RuntimeError("failed")))
    assert energy.predict_stability(parent,"H","12","S","A")["status"] == "unavailable"


def test_cysteine_warning(parent, monkeypatch):
    mock_model(monkeypatch)
    result = energy.predict_stability(parent,"H","12","S","C")
    assert any("Cysteine" in warning for warning in result["limitations"])


def test_checkpoint_hash_mismatch_detected(tmp_path, monkeypatch):
    monkeypatch.setenv("ATLAS_THERMOMPNN_DIR", str(tmp_path))
    (tmp_path / "LICENSE").write_text("not upstream")
    with pytest.raises(ValueError, match="Missing or modified"):
        energy._verified_manifest()


def test_capability_rejects_modified_model(monkeypatch):
    def modified():
        raise ValueError("Missing or modified ThermoMPNN file: checkpoint")
    monkeypatch.setattr(energy, "_verified_manifest", modified)
    result = energy.stability_capability()
    assert result["available"] is False
    assert "modified" in result["reason"]
