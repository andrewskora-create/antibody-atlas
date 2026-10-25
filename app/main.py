"""Local antibody analysis API. Bind to 127.0.0.1 only."""
from __future__ import annotations

import importlib.util
import hashlib
import copy
import tempfile
import threading
import uuid
from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError
from pathlib import Path
from typing import Literal
from urllib.parse import quote, urlsplit

from fastapi import FastAPI, File, Form, HTTPException, UploadFile, Query
from fastapi.responses import FileResponse, Response, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from .analysis import clean_sequence, sequence_metrics
from .translation import resolve as resolve_sequence
from .humanness import framework_identity
from .numbering import number_chain
from .reference_data import ReferenceDataMissing
from .residue_table import build_residue_table
from .molecule_flags import molecule_flags, DEFAULT_FORMULATION_PH
from . import deepsp as deepsp_module
from . import structure_ensemble
from . import constant_regions
from .structure_analysis import analyze_structure, extract_trace, summarize_site_priorities
from .report import make_report
from .report_brief import make_brief
from . import propermab
from . import eligibility
from .profiler_interpretation import annotate_antiberty, annotate_propermab
from .ptm_precedents import sequence_precedents
from .clinical_frameworks import clinical_framework_precedents
from .cdr_substitutions import propose as propose_substitutions
from . import paratope as paratope_module
from .attribution import citation as method_citation, all_citations
from .sequence_panel import annotate_chain, conclusion as panel_conclusion
from .source_provenance import source_records
from .humanness import REFERENCE_PATH

ROOT = Path(__file__).resolve().parents[1]
# TemporaryDirectory is mode 0700 and is removed when this server exits normally.
RUNS_HANDLE = tempfile.TemporaryDirectory(prefix="antibody-atlas-")
RUNS = Path(RUNS_HANDLE.name)
app = FastAPI(title="Antibody Atlas", version="0.1.0")
app.mount("/static", StaticFiles(directory=ROOT / "static"), name="static")
jobs: dict[str, dict] = {}
analysis_records: dict[str, dict] = {}
jobs_lock = threading.Lock()
fold_lock = threading.Lock()
antiberty_lock = threading.Lock()
antiberty_jobs: dict[str, dict] = {}
propermab_jobs: dict[str, dict] = {}
propermab_lock = threading.Lock()
ANTIBERTY_WEIGHTS = ROOT / "models" / "antiberty"


def antiberty_weights() -> str | None:
    """Resolve AntiBERTy weights: the bundled copy, else the package's own.

    The checkpoint is MIT-licensed but 104 MB, so it is not committed. The
    antiberty package resolves and downloads its own weights, so a clean clone
    works without a bundled copy; requiring one made the feature look broken.
    """
    if (ANTIBERTY_WEIGHTS / "config.json").is_file() and (ANTIBERTY_WEIGHTS / "model.safetensors").is_file():
        return str(ANTIBERTY_WEIGHTS)
    try:
        from antiberty import get_weights
        return get_weights(download=False)
    except Exception:
        return None


@app.middleware("http")
async def no_cache(request, call_next):
    # Binding to loopback alone does not prevent browser DNS rebinding. Only
    # loopback Host values and same-origin browser requests may reach the API.
    host = request.headers.get("host", "")
    try:
        authority = urlsplit("//" + host)
        port = authority.port
        valid_authority = (bool(authority.hostname) and authority.username is None
                           and authority.password is None and not authority.path
                           and not authority.query and not authority.fragment)
        local_host = authority.hostname in {"localhost", "127.0.0.1", "::1"}
        # Starlette's in-process transport is not reachable over the network.
        test_transport = (authority.hostname == "testserver"
                          and request.client is not None
                          and request.client.host == "testclient")
    except ValueError:
        valid_authority = local_host = test_transport = False
    error = None
    if not valid_authority or not (local_host or test_transport):
        error = (400, "Use the local Atlas address (127.0.0.1 or localhost).")
    elif (request.headers.get("sec-fetch-site", "").lower() == "cross-site"
          and not (request.method in {"GET", "HEAD"} and request.url.path == "/"
                   and request.headers.get("sec-fetch-mode") == "navigate"
                   and request.headers.get("sec-fetch-dest") == "document")):
        error = (403, "Cross-site browser requests are not allowed.")
    elif "origin" in request.headers:
        try:
            origin = urlsplit(request.headers["origin"])
            expected_port = port or (443 if request.url.scheme == "https" else 80)
            origin_port = origin.port or (443 if origin.scheme == "https" else 80)
            same_origin = (origin.scheme == request.url.scheme
                           and origin.hostname == authority.hostname
                           and origin_port == expected_port
                           and origin.username is None and origin.password is None
                           and not origin.path and not origin.query and not origin.fragment)
        except ValueError:
            same_origin = False
        if not same_origin:
            error = (403, "Only same-origin browser requests are allowed.")
    if error:
        response = JSONResponse({"detail": error[1]}, status_code=error[0])
    else:
        response = await call_next(request)
    response.headers["Cache-Control"] = "no-store"
    response.headers["Referrer-Policy"] = "no-referrer"
    return response


class AnalyzeRequest(BaseModel):
    heavy: str = Field(default="", max_length=1200)
    light: str = Field(default="", max_length=1200)
    scheme: Literal["imgt", "kabat"] = "imgt"
    ph: float = Field(default=7.4, ge=0, le=14)
    # Net charge is reported at ``ph``; the colloidal margin is evaluated at the
    # formulation pH, which for mAb drug products is usually acidic.
    formulation_ph: float = Field(default=DEFAULT_FORMULATION_PH, ge=0, le=14)
    # Attach constant regions so charge descriptors are computed on a full-length
    # molecule. None keeps the variable-domain-only basis.
    heavy_isotype: str | None = Field(default=None, max_length=12)
    light_constant: str | None = Field(default=None, max_length=12)
    clip_c_terminal_lysine: bool = True


class FoldRequest(BaseModel):
    heavy: str = Field(default="", max_length=1200)
    light: str = Field(default="", max_length=1200)
    ph: float = Field(default=7.4, ge=0, le=14)
    scheme: Literal["imgt", "kabat"] = "imgt"
    analysis_id: str | None = None
    # Number of ABodyBuilder2 models to fold. 1 keeps the historical behaviour
    # and the historical cost; 3 is ABodyBuilder2's and TAP2's recommendation,
    # and roughly triples fold time. Capped at 4 because that is how many
    # models exist.
    models: int = Field(default=1, ge=1, le=4)


class SubstitutionRequest(BaseModel):
    heavy: str = Field(default="", max_length=1200)
    light: str = Field(default="", max_length=1200)
    ph: float = Field(default=7.4, ge=0, le=14)
    scheme: Literal["imgt", "kabat"] = "imgt"
    # When a completed structure job is given, exposure and paratope terms are
    # computed. Without one the chemistry and sequence-level terms still run.
    fold_job_id: str | None = None
    cdr_only: bool = True
    max_per_site: int = Field(default=4, ge=1, le=19)
    preserve: dict[str, list[int]] | None = None


class EligibilityRequest(BaseModel):
    statement: str = Field(default="", max_length=200)


class StructureReportImage(BaseModel):
    png: str = Field(max_length=6_000_000)
    structure_sha256: str = Field(min_length=64, max_length=64)


class ReportRequest(BaseModel):
    structure_image: StructureReportImage | None = None
    title: str = Field(default="Antibody developability report", max_length=100)
    analysis_id: str
    fold_job_id: str | None = None
    antiberty_job_id: str | None = None
    propermab_job_id: str | None = None
    herceptin_tolerance_percent: float = Field(default=10, ge=0, le=20)
    #: "full" is the five-page report with the evidence and the method
    #: citations; "brief" is the one-page summary -- the verdict tiles, the
    #: key numbers, the flagged sites and what to measure next. The brief
    #: generator has existed since early on and was never reachable, so the
    #: cleanest output this project produces could not actually be obtained.
    style: Literal["full", "brief"] = "full"


def _version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unavailable"


def _sequence_digest(sequences: dict[str, str]) -> str:
    canonical = "\n".join(f"{chain}:{sequences[chain]}" for chain in sorted(sequences))
    return hashlib.sha256(canonical.encode()).hexdigest()


def _matched_structure(job_id: str, sequences: dict[str, str], ph: float,
                       scheme: str) -> dict:
    """Reject stale coordinates before joining them to another sequence."""
    with jobs_lock:
        job = copy.deepcopy(jobs.get(job_id))
        record = analysis_records.get((job or {}).get("analysis_id"))
        record = copy.deepcopy(record) if record else None
    if not job or job.get("status") != "complete":
        raise HTTPException(422, "The referenced structure job is not complete in this server session.")
    digest = job.get("input_sha256") or (record or {}).get("provenance", {}).get("input_sha256")
    job_ph = job.get("ph", (record or {}).get("ph"))
    if digest != _sequence_digest(sequences) or job_ph != ph or job.get("scheme") != scheme:
        raise HTTPException(422, "Structure input, pH and numbering convention must match the supplied sequences. Reanalyze and select their own structure.")
    return job


@app.get("/")
def index():
    """One page. The structure viewer, full tables, downloads and method
    sources are collapsible sections on it rather than a separate workbench,
    so there is never a second view of the same antibody saying something
    different."""
    return FileResponse(ROOT / "static" / "start.html")


@app.get("/api/health")
def health():
    return {
        "status": "ok",
        "numbering_model": importlib.util.find_spec("anarcii") is not None,
        "structure_model": importlib.util.find_spec("ImmuneBuilder") is not None,
        "structure_engine": "ImmuneBuilder ABodyBuilder2",
        "surface_engine": "FreeSASA",
        "antiberty_model": importlib.util.find_spec("antiberty") is not None
        and antiberty_weights() is not None,
        "propermab": {**propermab.availability(), "eligibility": eligibility.status()},
        "mode": "local single-user",
    }


@app.post("/api/analyze")
def analyze(request: AnalyzeRequest):
    if not request.heavy.strip() and not request.light.strip():
        raise HTTPException(422, "Enter a heavy or light variable-domain sequence.")
    result = {"scheme": request.scheme, "ph": request.ph, "chains": {},
              "methods": {
                  "numbering": "ANARCII with matching IMGT or Kabat CDR boundaries",
                  "humanness": "Closest functional human V-gene framework identity (OGRDB)",
                  "biophysical": "Biopython ProtParam and disclosed sequence-window heuristic",
                  "ptm": "Sequence motifs and susceptible residues; no occupancy or reaction-rate prediction",
              }}
    for label, raw in (("H", request.heavy), ("L", request.light)):
        if not raw.strip():
            continue
        try:
            # Translate nucleotide input, and record it: a silent translation
            # is worse than none, because the user cannot check the frame.
            raw_sequence, translation = resolve_sequence(raw)
            seq = clean_sequence(raw_sequence, allow_nucleotide=False)
            numbering = number_chain(seq, request.scheme)
            actual = numbering["chain_type"]
            if label == "H" and actual != "H":
                raise ValueError("The heavy-chain field was classified as a light chain.")
            if label == "L" and actual not in {"K", "L"}:
                raise ValueError("The light-chain field was classified as a heavy chain.")
            imgt = numbering if request.scheme == "imgt" else number_chain(seq, "imgt")
            result["chains"][label] = {
                "sequence": seq,
                "translation": translation,
                "numbering": numbering,
                # The residue table joins every layer onto canonical IMGT, so the
                # canonical numbering is kept even when the user selected Kabat.
                "numbering_imgt": imgt,
                "humanness": framework_identity(seq, imgt),
                "metrics": sequence_metrics(seq, request.ph, imgt),
            }
        except ValueError as exc:
            # Say what was received, not just that it was not recognised.
            from .translation import describe
            detail = f"{label} chain: {exc}"
            try:
                detail += ' ' + describe(raw)
            except Exception:
                pass
            raise HTTPException(422, detail) from exc
        except ImportError as exc:
            raise HTTPException(503, f"Numbering dependency unavailable: {exc}") from exc
    # Precedent lookups depend on corpora this project cannot redistribute, so a
    # fresh clone does not have them. They are genuinely optional: numbering,
    # liabilities, the residue table and the molecule-level flags do not use
    # them. Missing data therefore degrades those two sections rather than
    # failing the screen, and the response carries the command to build them.
    def _optional(fn, label):
        try:
            return fn(result["chains"], request.scheme)
        except ReferenceDataMissing as exc:
            result.setdefault("unavailable", []).append(
                {"feature": label, "reason": str(exc).split("\n")[0],
                 "build_command": next((line.strip() for line in str(exc).splitlines()
                                        if line.strip().startswith("python ")), None)})
            return {"cases": [], "sites": [], "unavailable": True, "reason": str(exc)}

    result["ptm_precedents"] = _optional(sequence_precedents, "Documented PTM precedents")
    result["clinical_framework_precedents"] = _optional(
        clinical_framework_precedents, "Approved and Phase-III comparator search")
    result["sources"] = source_records()
    for chain_id, chain in result["chains"].items():
        for site in chain["metrics"]["ptm_sites"]:
            site["precedent_ids"] = [case["id"] for case in result["ptm_precedents"]["cases"]
                if case["chemistry"] == site["kind"] and any(
                    m["query_chain"] == chain_id and m["query_sequence_position"] == site["start"]
                    for m in case["matches"])]
    # Band each sequence descriptor against the clinical reference panel, so the
    # interface can state a conclusion instead of printing a bare number.
    bands = {}
    for label, chain in result["chains"].items():
        chain_type = (chain.get("numbering") or {}).get("chain_type") or label
        chain["sequence_bands"] = annotate_chain(chain.get("metrics", {}), chain_type)
        bands[label] = chain["sequence_bands"]
    result["sequence_conclusion"] = panel_conclusion(bands)

    # The consolidated per-residue table. Sequence-only here: exposure and
    # paratope layers join later, once a structure job has run. The clinical
    # comparator search is already done by this point, so it is passed in
    # rather than left for the fix proposals to show on their own: whether an
    # approved drug carries the same motif belongs on the row where the
    # engineering call is made.
    try:
        result["residue_table"] = build_residue_table(
            result["chains"], precedents=result.get("clinical_framework_precedents"))
    except ValueError as exc:
        result["residue_table"] = {"error": str(exc), "rows": []}

    # Molecule-level charge and hydrophobicity, which are not per-residue and so
    # are reported separately rather than folded into the table.
    assembly = None
    if request.heavy_isotype:
        try:
            assembly = constant_regions.assemble(
                result["chains"], heavy_isotype=request.heavy_isotype,
                light_constant=request.light_constant,
                clip_c_terminal_lysine=request.clip_c_terminal_lysine)
        except constant_regions.ConstantRegionError as exc:
            raise HTTPException(422, str(exc)) from exc
    try:
        result["molecule_flags"] = molecule_flags(
            result["chains"], formulation_ph=request.formulation_ph, assay_ph=request.ph,
            assembly=assembly)
    except ValueError as exc:
        result["molecule_flags"] = {"error": str(exc), "flags": []}

    # DeepSP: the same physical properties the structure route measures, but
    # predicted from sequence in milliseconds and under an MIT licence rather
    # than Regeneron's academic-only one. It runs here, before any folding, so
    # a user who never starts a structure job still gets an aggregation and
    # charge read. It is a paired-Fv model: heavy-only input is skipped, not
    # approximated.
    if deepsp_module.available():
        try:
            result["deepsp"] = deepsp_module.analyze(
                (result["chains"].get("H") or {}).get("numbering_imgt"),
                (result["chains"].get("L") or {}).get("numbering_imgt"))
        except (ValueError, deepsp_module.DeepSPUnavailable) as exc:
            result["deepsp"] = {"status": "error", "reason": str(exc)}
    else:
        result["deepsp"] = {
            "status": "not_installed",
            "reason": ("DeepSP weights are not installed. They are MIT-licensed and "
                       "fetched rather than committed because they are 17 MB."),
            "build_command": "python setup/fetch_deepsp.py"}

    # Each rendered section carries its own attribution, so the citation sits
    # with the claim rather than in a methods appendix at the end.
    result["citations"] = {key: method_citation(key) for key in
                           ("ptm_priority", "biophysical", "humanness", "surface_patches",
                            "substitutions", "paratope", "framework_comparator",
                            "full_length_assembly", "deepsp", "structure_ensemble")}

    result["methods"]["residue_table"] = (
        "Every annotation joined onto canonical IMGT numbering, one row per residue")
    result["methods"]["deepsp"] = (
        "Spatial aggregation propensity and spatial charge maps predicted from the paired "
        "Fv sequence alone; no structure, no Regeneron licence")
    result["methods"]["molecule_flags"] = (
        "pI, net charge and CDR hydrophobicity mapped to aggregation, non-specific binding, "
        "clearance and viscosity; immunogenicity deliberately excluded")

    sequences = {label: chain["sequence"] for label, chain in result["chains"].items()}
    analysis_id = uuid.uuid4().hex
    result["analysis_id"] = analysis_id
    result["provenance"] = {
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": _sequence_digest(sequences),
        "reference_sha256": hashlib.sha256(REFERENCE_PATH.read_bytes()).hexdigest(),
        "reference_file": REFERENCE_PATH.name,
        "numbering_model_version": _version("anarcii"),
        "analysis_version": app.version,
        "variant_lineage": "Not supplied",
        "assay_provenance": "No measured assays supplied",
    }
    with jobs_lock:
        analysis_records[analysis_id] = copy.deepcopy(result)
    return result


def _fold(job_id: str, sequences: dict[str, str], ph: float, scheme: str,
          models: int = 1) -> None:
    with fold_lock:
        with jobs_lock:
            jobs[job_id].update(status="running", input_sha256=_sequence_digest(sequences),
                                ph=ph, scheme=scheme)
        try:
            paired = "L" in sequences
            kind = "antibody" if paired else "nanobody"
            engine = "ABodyBuilder2" if paired else "NanoBodyBuilder2"
            ensemble_plan = structure_ensemble.plan(models, kind) if models > 1 else None

            if ensemble_plan and len(ensemble_plan["model_ids"]) > 1:
                # Every model is folded and analysed separately; model 1 remains
                # the structure shown and the source of the site list, so the
                # viewer and the PDB download mean what they always meant.
                produced = structure_ensemble.fold_models(
                    sequences, ensemble_plan["model_ids"], RUNS, job_id)
                primary = next((p for p in produced if p["model_id"] == ensemble_plan["model_ids"][0]
                                and p.get("path")), None)
                if primary is None:
                    # Report what actually went wrong. Collapsing every folding
                    # failure into "no PDB file" hides the common cause -- a
                    # missing hmmscan, which legacy ANARCI needs and which
                    # ImmuneBuilder numbers through -- behind a message that
                    # sends the reader looking at the wrong thing entirely.
                    reasons = "; ".join(f"model {p['model_id']}: {p['error']}"
                                        for p in produced if p.get("error")) or "no reason reported"
                    raise RuntimeError(f"No structure model produced a PDB file. {reasons}")
                path = RUNS / f"{job_id}.pdb"
                path.write_bytes(primary["path"].read_bytes())
                method = (f"ImmuneBuilder {engine} ensemble of "
                          f"{len([p for p in produced if p.get('path')])} models "
                          f"({', '.join(str(p['model_id']) for p in produced if p.get('path'))}); "
                          f"exposure is the worst case across them")
            else:
                from ImmuneBuilder import ABodyBuilder2, NanoBodyBuilder2
                predictor = (ABodyBuilder2 if paired else NanoBodyBuilder2)(model_ids=[1])
                path = RUNS / f"{job_id}.pdb"
                predictor.predict(sequences).save(str(path))
                produced = [{"model_id": 1, "path": path}]
                method = f"ImmuneBuilder {engine} (single model)"
                if ensemble_plan:
                    method += "; ensemble requested but only one usable checkpoint"

            if not path.exists() or path.stat().st_size == 0:
                raise RuntimeError("The model did not produce a PDB file.")
            model_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()

            try:
                structure_analysis = analyze_structure(path, sequences, ph, scheme)
                analysis_error = None
                structure_trace = None
            except Exception as exc:
                structure_analysis = None
                analysis_error = f"Structure surface analysis failed: {str(exc)[:400]}"
                try:
                    structure_trace = extract_trace(path)
                except Exception:
                    structure_trace = None

            # Hydrogen-bond geometry. Cheap enough (sub-second) to run on every
            # fold, unlike repacking, and it is what turns "this position has a
            # liability" into "this position is also holding the fold together".
            if structure_analysis is not None:
                try:
                    from .hbonds import (find_hydrogen_bonds, protonate,
                                         summarize_position, packing_density)
                    topology, coords = protonate(path, ph)
                    bonds = find_hydrogen_bonds(topology, coords)
                    exposure = {(r["chain"], r["sequence_position"]):
                                r.get("relative_sidechain_sasa")
                                for r in structure_analysis.get("residues", [])
                                if r.get("sequence_position") is not None}
                    by_residue = {}
                    for residue in structure_analysis.get("residues", []):
                        position = residue.get("sequence_position")
                        if position is None:
                            continue
                        # pdb_number already carries any insertion code
                        # ("112A"); parsing it to an int would merge it with
                        # 112, which is a different residue.
                        key = f'{residue["chain"]}{residue["pdb_number"]}'
                        by_residue[key] = summarize_position(
                            bonds, residue["chain"], residue["pdb_number"],
                            relative_sasa=exposure.get((residue["chain"], position)),
                            source="predicted")
                    structure_analysis["hydrogen_bonds"] = {
                        "status": "ok",
                        "source": "predicted",
                        "total": len(bonds),
                        "by_kind": {k: sum(1 for b in bonds if b.kind == k)
                                    for k in ("backbone-backbone", "backbone-sidechain",
                                              "sidechain-sidechain")},
                        "by_residue": by_residue,
                    }
                except Exception as exc:
                    # A geometry failure must not cost the user the structure
                    # analysis they already have.
                    structure_analysis["hydrogen_bonds"] = {
                        "status": "unavailable",
                        "reason": f"{type(exc).__name__}: {str(exc)[:200]}"}

            # Reconcile exposure across the ensemble. A failure here must not
            # lose the primary result the user already has, so it degrades to a
            # note rather than failing the job.
            ensemble = None
            if structure_analysis is not None and len([p for p in produced if p.get("path")]) > 1:
                try:
                    per_model = []
                    for item in produced:
                        if not item.get("path"):
                            continue
                        if item["model_id"] == produced[0]["model_id"]:
                            per_model.append((item["model_id"], structure_analysis))
                        else:
                            per_model.append((item["model_id"],
                                              analyze_structure(item["path"], sequences, ph, scheme)))
                    ensemble = structure_ensemble.consensus(per_model)
                    ensemble["plan"] = ensemble_plan
                    ensemble["failed_models"] = [
                        {"model_id": p["model_id"], "error": p.get("error")}
                        for p in produced if not p.get("path")]
                    # The key is "ptm_sites", not "sites". Getting this wrong
                    # cost nothing visible: reconcile_sites was handed an empty
                    # list, returned an empty list, and the ensemble's
                    # worst-case exposure reached no PTM site at all. Every
                    # unit test still passed, because they call reconcile_sites
                    # directly with synthetic input. Only an end-to-end run
                    # showed the feature doing nothing.
                    structure_analysis["ptm_sites"] = structure_ensemble.reconcile_sites(
                        structure_analysis.get("ptm_sites", []), ensemble["residues"])
                    structure_analysis.setdefault("overall_ptm_risk", {}).update(
                        summarize_site_priorities(structure_analysis["ptm_sites"],
                            structure_analysis.get("unassessed_ptm_sites", [])))
                    structure_analysis["ensemble"] = {
                        k: v for k, v in ensemble.items() if k != "residues"}
                except Exception as exc:
                    ensemble = {"status": "failed", "reason": f"{type(exc).__name__}: {str(exc)[:300]}",
                                "plan": ensemble_plan}

            with jobs_lock:
                jobs[job_id].update(status="complete", method=method,
                                    model_sha256=model_sha256,
                                    model_version=_version("ImmuneBuilder"),
                                    pdb_url=f"/api/fold/{job_id}/pdb",
                                    structure_analysis=structure_analysis,
                                    structure_trace=structure_trace,
                                    ensemble=ensemble,
                                    analysis_error=analysis_error)
        except Exception as exc:
            with jobs_lock:
                jobs[job_id].update(status="failed", error=str(exc)[:500])


@app.post("/api/structure")
async def upload_structure(
    file: UploadFile = File(...),
    heavy: str = Form(default=""),
    light: str = Form(default=""),
    ph: float = Form(default=7.4, ge=0, le=14),
    scheme: str = Form(default="imgt"),
    copy_index: int = Form(default=0, ge=0),   # not "copy": it shadows a BaseModel attribute
    model_index: int = Form(default=0, ge=0),
    analysis_id: str = Form(default=""),
):
    """Use an experimental structure instead of a predicted one.

    The result is a job in the same shape a fold produces, so every consumer
    that takes a ``fold_job_id`` -- substitutions, PROPERMAB, the paratope
    predictor, the viewer -- works against it unchanged. The difference is
    ``structure_source``, which is "experimental" here and travels with every
    finding downstream.

    An antigen, if the file contains one, is kept. That is what turns distance
    to the epitope from an inference into a measurement.
    """
    from . import structure_input

    raw = await file.read()
    if len(raw) > 80 * 1024 * 1024:
        raise HTTPException(413, "Structure files above 80 MB are not accepted.")
    if scheme not in {"imgt", "kabat"}:
        raise HTTPException(422, "Numbering scheme must be IMGT or Kabat.")

    suffix = Path(file.filename or "upload.pdb").suffix.lower() or ".pdb"
    job_id = uuid.uuid4().hex
    source = RUNS / f"{job_id}.upload{suffix}"
    # The upload is stored under the job id, but the report must name the file
    # the USER chose. Showing them "396a4b72...upload.cif" when they picked
    # "9pi9.cif" is a needless puzzle, and in a session with several uploads
    # it makes the panels impossible to tell apart.
    original_name = Path(file.filename or "uploaded structure").name
    prepared = RUNS / f"{job_id}.pdb"

    sequences = {}
    for label, text in (("H", heavy), ("L", light)):
        if text.strip():
            try:
                sequences[label] = clean_sequence(resolve_sequence(text)[0],
                                                  allow_nucleotide=False)
            except ValueError as exc:
                raise HTTPException(422, f"{label} chain: {exc}") from exc

    # A structure cannot be attached to a different sequence record, even if
    # the caller supplies an otherwise valid analysis id.
    if analysis_id:
        with jobs_lock:
            stored = copy.deepcopy(analysis_records.get(analysis_id))
        if not stored:
            raise HTTPException(422, "Analysis record expired; analyze again before attaching a structure.")
        expected = {label: chain["sequence"] for label, chain in stored["chains"].items()}
        if (sequences and sequences != expected) or stored["ph"] != ph or stored["scheme"] != scheme:
            raise HTTPException(422, "Structure input, pH and numbering convention must match the analyzed record.")
        sequences = expected

    source.write_bytes(raw)

    try:
        report = structure_input.prepare(
            source, prepared, copy=copy_index, model_index=model_index,
            sequences=sequences or None)
    except structure_input.StructureInputError as exc:
        source.unlink(missing_ok=True)
        raise HTTPException(422, str(exc)) from exc
    except Exception as exc:
        source.unlink(missing_ok=True)
        raise HTTPException(422, f"Could not read the structure: {exc}") from exc
    finally:
        # The upload itself is not kept: the prepared Fv is what the rest of
        # the app uses, and holding a second copy of someone's coordinates
        # serves no purpose.
        source.unlink(missing_ok=True)

    # Numbering comes from the structure when the user supplied no sequence,
    # so a file alone is enough to get started.
    if not sequences:
        inspected = structure_input.inspect(prepared)
        for chain in inspected.chains:
            if chain.role == "heavy":
                sequences.setdefault("H", chain.sequence)
            elif chain.role == "light":
                sequences.setdefault("L", chain.sequence)

    structure_analysis = None
    analysis_error = None
    try:
        structure_analysis = analyze_structure(prepared, sequences, ph, scheme)
    except Exception as exc:
        analysis_error = f"Structure surface analysis failed: {str(exc)[:400]}"

    epitope = None
    if report.get("antigen_chains"):
        try:
            epitope = structure_input.epitope_distances(
                prepared, report["antigen_chains"])
        except Exception as exc:
            epitope = {"status": "failed", "reason": str(exc)[:300]}

    with jobs_lock:
        report["source_file"] = original_name
        jobs[job_id] = {
            "id": job_id, "status": "complete",
            "analysis_id": analysis_id or None, "scheme": scheme,
            "ph": ph, "input_sha256": _sequence_digest(sequences),
            "structure_source": "experimental",
            "method": (f"Experimental structure supplied by the user "
                       f"({report['source_format']}, copy {report['copy_used'] + 1} "
                       f"of {report['copies_available']})"),
            "model_sha256": hashlib.sha256(prepared.read_bytes()).hexdigest(),
            "pdb_url": f"/api/fold/{job_id}/pdb",
            "structure_analysis": structure_analysis,
            "analysis_error": analysis_error,
            "ingestion": report,
            "epitope": epitope,
            "sequences_from_structure": sequences,
            "warning": ("Experimental coordinates, so the structure "
                        "prediction's error is gone. Resolution, refinement "
                        "and unresolved density remain, and a crystal or grid "
                        "is not the molecule in formulation."),
        }
    return jobs[job_id]


@app.post("/api/structure/inspect")
async def inspect_structure(file: UploadFile = File(...)):
    """What is in this file, without committing to it.

    Lets the interface show which chains were found, how many copies there
    are and whether an antigen is present, so the user can choose a copy
    before anything is analysed.
    """
    from . import structure_input

    raw = await file.read()
    if len(raw) > 80 * 1024 * 1024:
        raise HTTPException(413, "Structure files above 80 MB are not accepted.")
    suffix = Path(file.filename or "upload.pdb").suffix.lower() or ".pdb"
    temporary = RUNS / f"inspect-{uuid.uuid4().hex}{suffix}"
    temporary.write_bytes(raw)
    try:
        return structure_input.inspect(temporary).as_dict()
    except structure_input.StructureInputError as exc:
        raise HTTPException(422, str(exc)) from exc
    finally:
        temporary.unlink(missing_ok=True)


@app.post("/api/fold")
def fold(request: FoldRequest):
    if importlib.util.find_spec("ImmuneBuilder") is None:
        raise HTTPException(503, "Install the optional ImmuneBuilder dependency to predict local structures.")
    if not request.heavy.strip():
        raise HTTPException(422, "A heavy variable-domain sequence is required for structure prediction.")
    try:
        sequences = {"H": clean_sequence(request.heavy)}
        if request.light.strip():
            sequences["L"] = clean_sequence(request.light)
        # The folding engine expects Fv sequences, not full-length Ig chains.
        if any(not 90 <= len(seq) <= 150 for seq in sequences.values()):
            raise ValueError("Structure prediction expects 90–150 residues per variable domain.")
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    if request.analysis_id:
        with jobs_lock:
            record = analysis_records.get(request.analysis_id)
        if not record or record["provenance"]["input_sha256"] != _sequence_digest(sequences) or record["scheme"] != request.scheme or record["ph"] != request.ph:
            raise HTTPException(422, "Fold input, pH and CDR convention must match the analyzed record.")
    job_id = uuid.uuid4().hex
    with jobs_lock:
        jobs[job_id] = {"id": job_id, "status": "queued", "method": None,
                        "analysis_id": request.analysis_id, "scheme": request.scheme,
                        "ph": request.ph, "input_sha256": _sequence_digest(sequences),
                        "models_requested": request.models,
                        # Set here so every consumer can label its findings.
                        # Structure ingestion will set "experimental"; today
                        # every structure this app has is one it predicted.
                        "structure_source": "predicted",
                        "warning": (
                            "Predicted coordinates are a hypothesis, especially for flexible "
                            "CDR-H3. Agreement between ensemble models is not accuracy: they "
                            "share training data and can agree and all be wrong. Validate "
                            "experimentally."
                            if request.models > 1 else
                            "Single-model predicted coordinates are a hypothesis, especially for "
                            "flexible CDR-H3; no ensemble confidence estimate is available. "
                            "Validate experimentally.")}
    threading.Thread(target=_fold,
                     args=(job_id, sequences, request.ph, request.scheme, request.models),
                     daemon=True).start()
    return jobs[job_id]


@app.get("/api/fold/{job_id}")
def fold_status(job_id: str):
    with jobs_lock:
        job = jobs.get(job_id)
        if job is None:
            raise HTTPException(404, "Structure job not found in this server session.")
        return dict(job)


@app.get("/api/fold/{job_id}/pdb")
def fold_pdb(job_id: str):
    with jobs_lock:
        job = jobs.get(job_id)
    if not job or job["status"] != "complete":
        raise HTTPException(404, "Structure file is not available.")
    return FileResponse(RUNS / f"{job_id}.pdb", media_type="chemical/x-pdb",
                        filename=f"antibody-{job_id[:8]}.pdb")


def _run_antiberty(job_id: str, sequences: dict[str, str], scheme: str = "imgt") -> None:
    with antiberty_lock:
        with jobs_lock:
            antiberty_jobs[job_id]["status"] = "running"
        try:
            from antiberty import AntiBERTyRunner
            runner = AntiBERTyRunner(device="cpu", checkpoint_path=antiberty_weights())
            from .antiberty_analysis import analyze_antiberty
            result = analyze_antiberty(runner, sequences, scheme=scheme, batch_size=16)
            chains = {label: {"numbering": number_chain(seq, scheme)} for label, seq in sequences.items()}
            result = annotate_antiberty(result, chains=chains)
            with jobs_lock:
                antiberty_jobs[job_id].update(status="complete", results=result,
                                             model_version=_version("antiberty"))
        except Exception as exc:
            with jobs_lock:
                antiberty_jobs[job_id].update(status="failed", error=str(exc)[:300])


@app.post("/api/antiberty")
def antiberty_start(request: FoldRequest):
    if antiberty_weights() is None:
        raise HTTPException(503, "Local AntiBERTy model weights are not installed.")
    if importlib.util.find_spec("antiberty") is None:
        raise HTTPException(503, "The antiberty Python package is not installed.")
    try:
        sequences = {label: clean_sequence(raw) for label, raw in (("H", request.heavy), ("L", request.light))
                     if raw.strip()}
        if not sequences or any(not 90 <= len(seq) <= 150 for seq in sequences.values()):
            raise ValueError("AntiBERTy expects 90-150 amino acids per variable domain.")
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    if request.analysis_id:
        with jobs_lock:
            record = analysis_records.get(request.analysis_id)
        if not record or record["provenance"]["input_sha256"] != _sequence_digest(sequences):
            raise HTTPException(422, "AntiBERTy input must match the analyzed record.")
    scheme = record["scheme"] if request.analysis_id else request.scheme
    job_id = uuid.uuid4().hex
    with jobs_lock:
        antiberty_jobs[job_id] = {"id": job_id, "status": "queued", "analysis_id": request.analysis_id}
    threading.Thread(target=_run_antiberty, args=(job_id, sequences, scheme), daemon=True).start()
    return antiberty_jobs[job_id]


@app.get("/api/antiberty/{job_id}")
def antiberty_status(job_id: str):
    with jobs_lock:
        job = antiberty_jobs.get(job_id)
        if job is None:
            raise HTTPException(404, "AntiBERTy job not found in this server session.")
        return dict(job)


class PropermabRequest(BaseModel):
    analysis_id: str
    fold_job_id: str | None = None


def _run_propermab(job_id, record, pdb_path):
    with propermab_lock:
        with jobs_lock:
            propermab_jobs[job_id]["status"] = "running"
        try:
            repaired_path=RUNS / f'{job_id}-terminal-repaired.pdb'
            result = annotate_propermab(propermab.run(record, pdb_path, repair_output_path=repaired_path))
            if result.get('terminal_repair',{}).get('status')=='complete' and repaired_path.exists():
                result['terminal_repair']['pdb_url']=f'/api/propermab/{job_id}/repaired-pdb'
            with jobs_lock:
                propermab_jobs[job_id].update(status="complete", results=result)
        except Exception as exc:
            with jobs_lock:
                propermab_jobs[job_id].update(status="failed", error=str(exc)[:300])


@app.get("/api/propermab/eligibility")
def propermab_eligibility():
    """Whether PROPERMAB may be used here, and what affirming means."""
    return eligibility.status()


@app.post("/api/propermab/eligibility")
def propermab_affirm(request: EligibilityRequest):
    if not eligibility.installed():
        raise HTTPException(409, "PROPERMAB is not installed, so there is nothing to affirm.")
    try:
        eligibility.record(request.statement, source="application")
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    return eligibility.status()


@app.delete("/api/propermab/eligibility")
def propermab_revoke():
    """Withdraw the affirmation. PROPERMAB stops being used immediately."""
    eligibility.revoke()
    return eligibility.status()


@app.post("/api/propermab")
def propermab_start(request: PropermabRequest):
    # Academic-only licence, enforced here rather than only documented.
    blocked = eligibility.gate()
    if blocked:
        raise HTTPException(403, blocked)
    if not propermab.availability()["available"]:
        raise HTTPException(503, "Local PROPERMAB dependencies are unavailable.")
    with jobs_lock:
        record = copy.deepcopy(analysis_records.get(request.analysis_id))
        fold_job = copy.deepcopy(jobs.get(request.fold_job_id)) if request.fold_job_id else None
    if not record:
        raise HTTPException(404, "Analyze again before running PROPERMAB.")
    if set(record['chains']) != {'H','L'} or not 4 <= record['ph'] <= 10:
        raise HTTPException(422, "PROPERMAB requires paired VH/VL and pH 4–10.")
    if request.fold_job_id and (not fold_job or fold_job.get('status')!='complete' or fold_job.get('analysis_id')!=request.analysis_id):
        raise HTTPException(422, "The structure job does not belong to this analysis record.")
    job_id = uuid.uuid4().hex
    with jobs_lock:
        propermab_jobs[job_id] = {'id':job_id,'status':'queued','analysis_id':request.analysis_id,
                                 'fold_job_id':request.fold_job_id}
    pdb_path = RUNS/f'{request.fold_job_id}.pdb' if fold_job else None
    threading.Thread(target=_run_propermab,args=(job_id,record,pdb_path),daemon=True).start()
    return dict(propermab_jobs[job_id])


@app.get("/api/propermab/{job_id}")
def propermab_status(job_id: str):
    with jobs_lock:
        job = propermab_jobs.get(job_id)
        if job is None:
            raise HTTPException(404, "PROPERMAB job not found in this server session.")
        return copy.deepcopy(job)


@app.get('/api/propermab/{job_id}/repaired-pdb')
def repaired_pdb(job_id: str):
    with jobs_lock:
        job=propermab_jobs.get(job_id)
    if not job or job.get('status')!='complete' or job.get('results',{}).get('terminal_repair',{}).get('status')!='complete':
        raise HTTPException(404,'No completed repaired model for this job.')
    path=RUNS / f'{job_id}-terminal-repaired.pdb'
    if not path.exists():raise HTTPException(404,'Repaired model expired.')
    return FileResponse(path,media_type='chemical/x-pdb',filename='terminal-repaired-hypothesis.pdb')


def assemble_report(request: ReportRequest):
    with jobs_lock:
        stored = analysis_records.get(request.analysis_id)
        fold_job = jobs.get(request.fold_job_id) if request.fold_job_id else None
        anti_job = antiberty_jobs.get(request.antiberty_job_id) if request.antiberty_job_id else None
    if not stored:
        raise HTTPException(404, "Analysis record expired; analyze again before exporting.")
    report = copy.deepcopy(stored)
    if request.fold_job_id:
        if not fold_job or fold_job.get("status") != "complete" or fold_job.get("analysis_id") != request.analysis_id:
            raise HTTPException(422, "The structure job does not belong to this analysis record.")
        report["structure_analysis"] = fold_job.get("structure_analysis")
        report["provenance"].update(structure_model=fold_job.get("method"),
                                    structure_model_version=fold_job.get("model_version"),
                                    structure_pdb_sha256=fold_job.get("model_sha256"))
    if request.antiberty_job_id:
        if not anti_job or anti_job.get("status") != "complete" or anti_job.get("analysis_id") != request.analysis_id:
            raise HTTPException(422, "The AntiBERTy job does not belong to this analysis record.")
        for label, result in anti_job["results"].items():
            report["chains"][label]["antiberty"] = result
        report["provenance"]["antiberty_model_version"] = anti_job.get("model_version")
    if request.propermab_job_id:
        with jobs_lock:
            pm_job = copy.deepcopy(propermab_jobs.get(request.propermab_job_id))
        if not pm_job or pm_job.get("status") != "complete" or pm_job.get("analysis_id") != request.analysis_id:
            raise HTTPException(422, "The PROPERMAB job does not belong to this analysis record.")
        if request.fold_job_id and pm_job.get("fold_job_id") and request.fold_job_id != pm_job["fold_job_id"]:
            raise HTTPException(422, "PROPERMAB was calculated from a different structure job; rerun it for the selected structure.")
        report["propermab"] = pm_job["results"]
        from .model_selection import select_propermab_structure
        try:
            select_propermab_structure(report,pm_job['results'])
        except ValueError as exc:
            raise HTTPException(422,str(exc)) from exc
        report["provenance"]["propermab_source_commit"] = pm_job["results"]["source_commit"]
    report["assessment_policy"] = {"herceptin_tolerance_percent": request.herceptin_tolerance_percent}
    report["provenance"]["report_timestamp_utc"] = datetime.now(timezone.utc).isoformat()
    return report


@app.get("/api/methods")
def methods():
    """Attribution for every assessment this application produces."""
    return {"assessments": all_citations(),
            "note": "Citation is not evidence that a method ran, and no cited paper "
                    "validates this application's combined heuristics."}


@app.get("/api/substitutions/capabilities")
def substitution_capabilities():
    """What the PTM-fix engineering step can do in this installation."""
    from .reference_data import missing_datasets
    from .stability_energy import stability_capability
    missing = missing_datasets()
    # What this RUNNING PROCESS can do, not what the files on disk say.
    #
    # Static files are re-read on every request; Python modules are loaded
    # once at startup. Update the code without restarting and the browser
    # gets the new page against the old server, so the page offers a control
    # that posts to an endpoint this process has never heard of. The failure
    # surfaces as an opaque 404 in the middle of an analysis.
    #
    # Routes are read from the live app rather than hardcoded, so this stays
    # honest without anyone remembering to bump a version number.
    served = {getattr(route, "path", None) for route in app.routes}
    return {
        "paratope": paratope_module.available(),
        "stability_energy": stability_capability(),
        "clinical_precedent_available": "clinical_frameworks" not in missing,
        "missing_reference_data": missing,
        "setup_hint": ("python setup/build_reference_data.py --all" if missing else None),
        "endpoints": sorted(p for p in served if isinstance(p, str) and p.startswith("/api/")),
        "features": {
            "structure_upload": "/api/structure" in served,
            "repack": "/api/repack" in served,
            "stability_prediction": True,
            "top_five_models": True,
            "deepsp": deepsp_module.available(),
            "structure_ensemble": bool(structure_ensemble.usable_models()[0]),
        },
    }


@app.post("/api/substitutions")
def substitutions(request: SubstitutionRequest):
    """Propose and rank PTM fixes, with charge, aggregation and paratope context."""
    try:
        sequences = {label: clean_sequence(raw)
                     for label, raw in (("H", request.heavy), ("L", request.light))
                     if raw.strip()}
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    if not sequences:
        raise HTTPException(422, "Enter a heavy or light variable-domain sequence.")

    chains: dict[str, dict] = {}
    for label, sequence in sequences.items():
        try:
            canonical = number_chain(sequence, "imgt")
            selected = canonical if request.scheme == "imgt" else number_chain(sequence, request.scheme)
        except ValueError as exc:
            raise HTTPException(422, f"Chain {label}: {exc}") from exc
        chains[label] = {"sequence": sequence, "numbering": selected,
                         "numbering_imgt": canonical,
                         "metrics": sequence_metrics(sequence, request.ph, canonical),
                         "humanness": framework_identity(sequence, canonical)}

    structure_analysis = None
    paratope = None
    structure_note = "No structure was supplied, so exposure and paratope terms were not computed."
    if request.fold_job_id:
        job = _matched_structure(request.fold_job_id, sequences, request.ph, request.scheme)
        structure_analysis = job.get("structure_analysis")
        pdb_path = RUNS / f"{request.fold_job_id}.pdb"
        if pdb_path.is_file():
            paratope = paratope_module.predict(
                pdb_path, heavy_chain_id="H",
                light_chain_id="L" if "L" in sequences else None,
                imgt_to_position=paratope_module.imgt_position_map(chains))
            structure_note = (("Exposure calculated from supplied experimental coordinates. "
                               if job.get("structure_source") == "experimental"
                               else "Exposure calculated from a predicted structure. ")
                              + f"Paratope prediction {paratope.get('status')}.")
        else:
            structure_note = "Structure job is complete but its coordinates are no longer on disk."

    preserve = ({label: set(positions) for label, positions in request.preserve.items()}
                if request.preserve else None)
    try:
        result = propose_substitutions(
            chains, cdr_only=request.cdr_only, max_per_site=request.max_per_site,
            structure=structure_analysis, paratope=paratope, ph=request.ph,
            preserve=preserve)
    except Exception as exc:
        from .reference_data import ReferenceDataMissing
        if isinstance(exc, ReferenceDataMissing):
            raise HTTPException(503, str(exc)) from exc
        raise HTTPException(500, f"Substitution proposal failed: {str(exc)[:300]}") from exc

    # Framework liabilities are searched against approved and Phase-III drugs
    # automatically, so the interface can show the comparator alongside the site.
    precedents = None
    try:
        from .clinical_frameworks import clinical_framework_precedents
        precedents = clinical_framework_precedents(chains, request.scheme)
        index = {(p.get('query_chain'), p.get('position')): p for p in precedents.get('sites', [])}
        for label, chain_result in result.get("chains", {}).items():
            for site in chain_result.get("sites", []):
                entry = index.get((label, site.get("position")))
                matches = (entry or {}).get("matches", [])
                approved = {m["drug"] for m in matches if m.get("stage") == "Approved"}
                phase_iii = {m["drug"] for m in matches if m.get("stage") == "Phase-III"}
                site["comparators"] = {
                    "searched": entry is not None,
                    "approved": len(approved), "phase_iii": len(phase_iii),
                    "drugs": [f"{m['drug']} ({m['stage']})" for m in matches[:6]],
                }
        result["comparator_coverage"] = precedents.get("coverage")
    except Exception as exc:  # reference data absent, or numbering unavailable
        result["comparator_coverage"] = {"error": str(exc)[:200]}

    result["citations"] = {key: method_citation(key) for key in
                           ("ptm_priority", "substitutions", "paratope", "biophysical",
                            "structure", "framework_comparator")}
    result["structure_note"] = structure_note
    result["paratope_detail"] = paratope

    # Rebuild the per-residue table here, where the structure and the paratope
    # are both in hand. The copy made during /api/analyze predates the fold, so
    # its exposure and contact columns are necessarily empty; serving that one
    # after a successful fold would show the table's two most useful columns as
    # blank on a run that actually has the data.
    try:
        result["residue_table"] = build_residue_table(
            chains, structure=structure_analysis, paratope=paratope,
            precedents=precedents)
    except ValueError as exc:
        result["residue_table"] = {"error": str(exc), "rows": []}
    result["scheme"] = request.scheme
    result["ph"] = request.ph
    return result


#: Which category tiles rest on coordinates, and which would read the same
#: without any structure at all. A single page-level "this is predicted"
#: banner flattens that: germline identity is sequence-based, while CDR
#: liability triage uses structural exposure as well as motif chemistry.
CATEGORY_BASIS = {
    "CDR liabilities": "structure",
    "Aggregation": "structure",
    "Charge surface": "structure",
    "Concentration": "structure",
    "Germline FR": "sequence",
}


@app.post("/api/assessment")
def assessment(request: ReportRequest):
    """The screen's conclusions, including the category tiles.

    The tiles come from ``report_brief.signals`` -- the same function that
    draws them on the one-page PDF -- rather than being recomputed here. If
    the page banded them itself the screen and the brief would drift apart,
    and a reader comparing the two would have no way to tell which was
    right. No threshold is invented at this layer.
    """
    from .assessment import build_assessment
    from .report_brief import signals

    report = assemble_report(request)
    built = build_assessment(report)
    source = "experimental" if (
        (jobs.get(request.fold_job_id) or {}).get("structure_source") == "experimental"
    ) else "predicted"
    try:
        lights, framework, compared, above = signals(report, built)
    except Exception as exc:                     # a missing layer, not a bug
        built["categories"] = {"available": False, "reason": str(exc)[:200]}
        return built
    built["categories"] = {
        "available": True,
        "tiles": [{
            "label": label,
            "state": state,
            "detail": detail,
            # "sequence" means a structure would not change this call, so it
            # is neither measured nor modelled and must not be badged as if
            # the question of coordinates applied to it.
            "basis": ("sequence" if CATEGORY_BASIS.get(label) == "sequence"
                      else source if state != "unknown" else "none"),
        } for label, state, detail in lights],
        "worst": ("high" if any(s == "high" for _, s, _ in lights)
                  else "review" if any(s == "review" for _, s, _ in lights)
                  else "unknown" if any(s == "unknown" for _, s, _ in lights)
                  else "good"),
        "framework_sites": len(framework),
        "panel": {"compared": compared, "above": above},
        "state_labels": {"good": "OK", "review": "REVIEW",
                         "high": "ATTENTION", "unknown": "NO DATA"},
        "unknown_note": ("NO DATA means the evidence for that category was not "
                         "produced, usually because no structure was supplied. "
                         "It is not a pass."),
    }
    return built


@app.post("/api/report")
def export_report(request: ReportRequest):
    report = assemble_report(request)
    try:
        from .report_image import report_image
        image = report_image(request.structure_image.model_dump() if request.structure_image else None, report.get("provenance",{}).get("structure_pdb_sha256"))
        pdf = (make_brief(report, request.title, structure_image=image)
               if request.style == "brief"
               else make_report(report, request.title, structure_image=image))
    except (KeyError, TypeError, ValueError) as exc:
        raise HTTPException(422, f"Report data is incomplete: {exc}") from exc
    return Response(content=pdf, media_type="application/pdf",
                    headers={"Content-Disposition":
                             f'attachment; filename="antibody-developability-'
                             f'{"brief" if request.style == "brief" else "report"}.pdf"'})

# Loaded after core entry points so the orchestration can reuse the same local engines.
from .variant_jobs import router as variant_router
app.include_router(variant_router)


class PositionRequest(BaseModel):
    heavy: str = Field(default="", max_length=1200)
    light: str = Field(default="", max_length=1200)
    chain: Literal["H", "L"] = "H"
    position: int = Field(ge=1, le=600)
    scheme: Literal["imgt", "kabat"] = "imgt"
    ph: float = Field(default=7.4, ge=0, le=14)
    fold_job_id: str | None = None
    max_candidates: int = Field(default=8, ge=1, le=19)
    #: How many of the top-ranked candidates to model rather than infer.
    #: Each costs roughly three seconds, so the whole set of nineteen would
    #: take about a minute; the default models the few a user would actually
    #: consider and leaves the rest to the sequence-level ranking.
    repack_top: int = Field(default=5, ge=0, le=19)


repack_previews: dict[str, dict] = {}
repack_results: dict[tuple, dict] = {}
repack_lock = threading.Lock()


class RepackRequest(BaseModel):
    """One substitution, modelled rather than inferred.

    Deliberately not part of /api/substitutions: a repack costs seconds where
    the rest of the ranking costs milliseconds, so running it over every
    candidate at every site would turn a screen into a batch job. It is asked
    for per substitution, once the ranking has narrowed the field.
    """
    fold_job_id: str
    chain: Literal["H", "L"] = "H"
    position: int = Field(ge=1, le=600)
    replacement: str = Field(min_length=1, max_length=1)
    ph: float = Field(default=7.4, ge=0, le=14)


@app.post("/api/repack")
def repack_substitution(request: RepackRequest):
    """Build the mutant side chain, relax it, and report the geometry change."""
    # Serialise the CPU-heavy calculation and reuse identical local previews.
    with repack_lock:
        return _repack_substitution(request)


def _repack_substitution(request: RepackRequest):
    from .repack import RepackError, repack as run_repack
    from .stability_energy import predict_stability

    with jobs_lock:
        job = jobs.get(request.fold_job_id)
    if not job or job.get("status") != "complete":
        raise HTTPException(404, "No completed structure for that fold job.")
    pdb_path = RUNS / f"{request.fold_job_id}.pdb"
    if not pdb_path.is_file():
        raise HTTPException(404, "The structure for that job is no longer available.")
    # Whether these coordinates were measured or predicted changes what the
    # answer is worth, so it travels with the result rather than being
    # inferred by whoever reads it.
    source = job.get("structure_source", "predicted")
    structure = job.get("structure_analysis") or {}
    if job.get("ph", structure.get("ph")) != request.ph:
        raise HTTPException(422, "Use the same pH as the analyzed parent structure.")
    matches = [r for r in structure.get("residues", [])
               if r.get("chain") == request.chain and r.get("sequence_position") == request.position
               and r.get("pdb_number")]
    if len(matches) != 1:
        raise HTTPException(422, "The selected sequence position has no unique coordinate mapping.")
    pdb_residue = str(matches[0]["pdb_number"])
    cache_key = (request.fold_job_id, hashlib.sha256(pdb_path.read_bytes()).hexdigest(),
                 request.chain, request.position, pdb_residue, request.replacement.upper(), request.ph, source)
    cached = repack_results.get(cache_key)
    if cached and all(Path(p).is_file() for p in cached['paths']):
        result = copy.deepcopy(cached['result'])
        # An engine may have been installed since a geometry-only preview.
        if result.get('free_energy', {}).get('status') != 'available':
            result['free_energy'] = predict_stability(pdb_path, request.chain, pdb_residue,
                                                      result.get('from', matches[0].get('aa', '')),
                                                      request.replacement.upper())
            cached['result'] = copy.deepcopy(result)
        return result
    preview_id = uuid.uuid4().hex
    paths = tuple(RUNS / f"repack-{preview_id}-{state}.pdb" for state in ("baseline", "candidate"))
    try:
        result = run_repack(pdb_path, request.chain, pdb_residue,
                            request.replacement, ph=request.ph, source=source, output_paths=paths)
        if not all(path.is_file() for path in paths):
            raise RuntimeError("Local preview coordinates were not produced.")
    except Exception as exc:
        for path in paths:
            path.unlink(missing_ok=True)
        if isinstance(exc, RepackError):
            raise HTTPException(422, str(exc)) from exc
        raise HTTPException(500, f"Repacking failed: {type(exc).__name__}: {exc}") from exc
    with jobs_lock:
        repack_previews[preview_id] = {"baseline": paths[0], "candidate": paths[1]}
    parent_url, mutant_url = (f"/api/repack/{preview_id}/{state}/pdb" for state in ("baseline", "candidate"))
    result.update(parent_pdb_url=parent_url, mutant_pdb_url=mutant_url,
                  sequence_position=request.position, ph=request.ph,
                  structures={"baseline": {"pdb_url": parent_url}, "candidate": {"pdb_url": mutant_url},
                      "selection": {"chain": request.chain, "residue": pdb_residue},
                      "source": source, "method": result["method"]},
                  free_energy=predict_stability(pdb_path, request.chain, pdb_residue,
                                                result.get('from', matches[0].get('aa', '')),
                                                request.replacement.upper()))
    result["citation"] = method_citation("repack")
    repack_results[cache_key] = {'result': copy.deepcopy(result), 'paths': paths}
    if len(repack_results) > 128:
        repack_results.pop(next(iter(repack_results)))
    return result


@app.get("/api/repack/{preview_id}/{state}/pdb")
def repack_pdb(preview_id: str, state: Literal["baseline", "candidate"]):
    """Serve only coordinates registered by this session's local preview job."""
    with jobs_lock:
        preview = repack_previews.get(preview_id)
    path = preview.get(state) if preview else None
    if path is None or not path.is_file():
        raise HTTPException(404, "Local preview coordinates are no longer available.")
    return FileResponse(path, media_type="chemical/x-pdb", filename=f"local-{state}-{preview_id[:8]}.pdb")


@app.post("/api/substitutions/position")
def substitutions_at_position(request: PositionRequest):
    """Rank substitutions at a position the user picked.

    The automatic proposals only cover detected liabilities. This answers the
    other question -- "what could go here instead" -- anywhere in the chain.
    """
    from .cdr_substitutions import propose_at_position

    try:
        sequences = {label: clean_sequence(resolve_sequence(raw)[0])
                     for label, raw in (("H", request.heavy), ("L", request.light))
                     if raw.strip()}
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    if request.chain not in sequences:
        raise HTTPException(422, f"No {request.chain} chain was supplied.")

    chains: dict[str, dict] = {}
    for label, sequence in sequences.items():
        try:
            canonical = number_chain(sequence, "imgt")
        except ValueError as exc:
            raise HTTPException(422, f"Chain {label}: {exc}") from exc
        selected = canonical if request.scheme == "imgt" else number_chain(sequence, request.scheme)
        chains[label] = {"sequence": sequence, "numbering": selected,
                         "numbering_imgt": canonical}

    structure_analysis = paratope = None
    if request.fold_job_id:
        job = _matched_structure(request.fold_job_id, sequences, request.ph, request.scheme)
        if job:
            structure_analysis = job.get("structure_analysis")
            pdb_path = RUNS / f"{request.fold_job_id}.pdb"
            if pdb_path.is_file():
                paratope = paratope_module.predict(
                    pdb_path, heavy_chain_id="H",
                    light_chain_id="L" if "L" in sequences else None,
                    imgt_to_position=paratope_module.imgt_position_map(chains))

    try:
        result = propose_at_position(
            chains, request.chain, request.position,
            max_candidates=max(request.max_candidates, request.repack_top), structure=structure_analysis,
            paratope=paratope, ph=request.ph)
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    result["citation"] = method_citation("substitutions")

    # Model the top few rather than inferring them. The sequence-level ranking
    # decides WHICH candidates are worth the compute; repacking then says what
    # each would actually cost in hydrogen bonds, which is the part subtraction
    # gets wrong on conservative swaps.
    if request.repack_top and request.fold_job_id and structure_analysis is not None:
        pdb_path = RUNS / f"{request.fold_job_id}.pdb"
        if pdb_path.is_file():
            with jobs_lock:
                job = jobs.get(request.fold_job_id) or {}
            source = job.get("structure_source", "predicted")
            modelled = []
            for candidate in (result.get("candidates") or [])[:request.repack_top]:
                replacement = candidate.get("replacement_aa")
                if not replacement:
                    continue
                try:
                    detail = repack_substitution(RepackRequest(
                        fold_job_id=request.fold_job_id, chain=request.chain,
                        position=request.position, replacement=replacement, ph=request.ph))
                except Exception as exc:
                    candidate["repack"] = {"status": "failed",
                                           "reason": f"{type(exc).__name__}: {str(exc)[:160]}"}
                    continue
                candidate["free_energy"] = detail.get("free_energy")
                hb = detail["hydrogen_bonds"]
                candidate["repack"] = {
                    "status": "ok",
                    "source": detail["source"],
                    "hydrogen_bonds_lost": hb["lost_at_site"],
                    "hydrogen_bonds_gained": hb["gained_at_site"],
                    "net": hb["net_at_site"],
                    "clashes_introduced": detail["clashes"]["introduced"],
                    "protocol_noise": detail["protocol_noise"][
                        "bonds_changed_away_from_the_site"],
                    "detail": detail,
                }
                modelled.append(replacement)
            result["repacked"] = {
                "modelled": modelled,
                "requested": request.repack_top,
                "source": source,
                "note": ("The top candidates selected using residue similarity and clinical frequency were built and relaxed rather than "
                         "scored by subtraction. Candidates below them carry the "
                         "sequence-level ranking only, which cannot see a compensating "
                         "bond forming and so under-rates conservative swaps."),
                "caveat": ("Hydrogen bonds and clashes are geometric estimates from "
                           + ("locally modeled coordinates based on an experimental input; the mutant is not an experimental observation." if source == "experimental"
                              else "a PREDICTED backbone, where side-chain placement is "
                                   "the least reliable part. Treat it as triage, not "
                                   "evidence.") + " Folding ΔΔG is a separate learned prediction from the parent chain, not from local relaxation energies."),
            }
    return result


class BatchRequest(BaseModel):
    panel: str = Field(max_length=400_000)
    scheme: Literal["imgt", "kabat"] = "imgt"
    ph: float = Field(default=7.4, ge=0, le=14)
    formulation_ph: float = Field(default=6.0, ge=0, le=14)
    isotype: str | None = Field(default=None, max_length=12)
    precedents: bool = True


@app.post("/api/batch")
def batch_start(request: BatchRequest):
    from . import batch_web
    try:
        return batch_web.start(request.panel, {
            "scheme": request.scheme, "ph": request.ph,
            "formulation_ph": request.formulation_ph,
            "isotype": request.isotype, "precedents": request.precedents})
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc


@app.get("/api/batch/{job_id}")
def batch_status(job_id: str):
    from . import batch_web
    job = batch_web.get(job_id)
    if not job:
        raise HTTPException(404, "That batch is not in this server session.")
    return job


@app.get("/api/batch/{job_id}/sequences/{name}")
def batch_sequences(job_id: str, name: str):
    """The input sequences for one molecule, so the page can promote a row
    into the full single-molecule screen without the user re-pasting."""
    from . import batch_web
    sequences = batch_web.sequences_for(job_id, name)
    if sequences is None:
        raise HTTPException(404, f"No construct {name!r} in that batch.")
    return {"name": name, "chains": sequences}


@app.get("/api/batch/{job_id}/report/{name}")
def batch_report(job_id: str, name: str):
    """The full developability report for one screened molecule."""
    from fastapi.responses import Response
    from . import batch_web
    from .report import make_report
    record = batch_web.record_for(job_id, name)
    if record is None:
        raise HTTPException(404, f"No completed result for {name!r} in that batch.")
    try:
        pdf = make_report(record, title=f"Developability report - {name}")
    except Exception as exc:
        raise HTTPException(500, f"Report generation failed: {str(exc)[:200]}") from exc
    return Response(pdf, media_type="application/pdf",
                    headers={"Content-Disposition":
                             'attachment; filename="developability-report.pdf"; '
                             f"filename*=UTF-8''{quote(name + '_developability.pdf', safe='')}"})


@app.get("/api/batch/{job_id}/reports")
def batch_reports(job_id: str, names: str = "", name: list[str] | None = Query(default=None)):
    """Reports for several molecules at once, as a zip."""
    import io
    import zipfile
    from fastapi.responses import Response
    from . import batch_web
    from .report import make_report

    job = batch_web.get(job_id)
    if not job:
        raise HTTPException(404, "That batch is not in this server session.")
    wanted = (name if name is not None else
              [n.strip() for n in names.split(",") if n.strip()] or
              [r["name"] for r in job.get("rows", [])])
    wanted = list(dict.fromkeys(wanted))

    buffer = io.BytesIO()
    written, missing = [], []
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name in wanted:
            record = batch_web.record_for(job_id, name)
            if record is None:
                missing.append(name)
                continue
            try:
                archive.writestr(f"{name}_developability.pdf",
                                 make_report(record, title=f"Developability report - {name}"))
                written.append(name)
            except Exception as exc:                 # one bad report must not lose the rest
                missing.append(f"{name} ({type(exc).__name__})")
        if missing:
            archive.writestr("NOT_INCLUDED.txt",
                             "No report was produced for:\n  " + "\n  ".join(missing) + "\n")
    if not written:
        raise HTTPException(404, "No completed results to report on.")
    return Response(buffer.getvalue(), media_type="application/zip",
                    headers={"Content-Disposition":
                             f'attachment; filename="batch_{job_id[:8]}_reports.zip"'})


@app.get("/api/batch/{job_id}/csv")
def batch_csv(job_id: str):
    from fastapi.responses import Response
    from . import batch_web
    if not batch_web.get(job_id):
        raise HTTPException(404, "That batch is not in this server session.")
    return Response(batch_web.to_csv(job_id), media_type="text/csv",
                    headers={"Content-Disposition":
                             f'attachment; filename="batch_{job_id[:8]}.csv"'})


@app.get("/api/constant-regions")
def constant_region_options():
    """What can be attached to the variable domains, for the isotype selector."""
    try:
        return constant_regions.available()
    except constant_regions.ConstantRegionError as exc:
        raise HTTPException(503, str(exc)) from exc
