"""Describe actual profiler scales against a small, protocol-matched local panel.

Comparisons are descriptive, never therapeutic pass/fail thresholds or population
percentiles. No network request or sequence comparison service is used.
"""
from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
from bisect import bisect_right
from statistics import median

COMPARATOR_PATH = Path(__file__).resolve().parents[1] / "data" / "profiler_comparators.json"
PANEL_BOUNDARY = ("This therapeutic reference panel was not selected for favorable experimental "
                  "behavior; it is every molecule meeting its stated criteria in the source "
                  "snapshot. Its range is the observed minimum to maximum under the same "
                  "calculation protocol, not an acceptance interval, and a percentile within it "
                  "is a position among approved antibodies, not a verdict on acceptability. "
                  "Within range means similar on this descriptor, not low risk, safe, or suitable "
                  "for high concentration; approval or Phase III status does not establish those "
                  "properties. Approved antibodies are themselves a survivorship-filtered "
                  "population.")


def _load_panel():
    try:
        panel = json.loads(COMPARATOR_PATH.read_text())
        return panel if isinstance(panel.get("references"), list) else {}
    except (OSError, ValueError, TypeError):
        return {}


def _panel_note(panel):
    names = [str(r.get('name', r.get('id', 'Unnamed'))) for r in panel.get('references', [])]
    text = panel.get('description', 'Local comparator panel unavailable.')
    return text + (' Panel members: ' + ', '.join(names) + '.' if names else '')


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


#: Below this many references, a percentile is noise dressed as precision: with
#: n references the finest distinction expressible is 100/n percentile points,
#: and the tails are decided by one or two molecules. The panel's own range is
#: reported either way; the percentile is withheld until it carries information.
MIN_REFERENCES_FOR_PERCENTILE = 60


def _percentile_of(value, sorted_values):
    """Fraction of references at or below `value`, as a percentage.

    The plain empirical definition, deliberately not interpolated. Interpolating
    between two adjacent approved antibodies would invent precision the panel
    does not have.
    """
    if not sorted_values:
        return None
    at_or_below = bisect_right(sorted_values, value)
    return round(100.0 * at_or_below / len(sorted_values), 1)


def _benchmark(key, label, value, unit, references, direction, scale, note=""):
    references = [r for r in references if _finite(r.get("value"))]
    values = [r["value"] for r in references]
    low, high = (min(values), max(values)) if values else (None, None)
    relation = "unavailable" if not values or not _finite(value) else "below_range" if value < low else "above_range" if value > high else "within_range"
    benchmark = {"key": key, "label": label, "query_value": value, "unit": unit,
                 "reference_n": len(values), "min": low, "median": median(values) if values else None,
                 "max": high, "query_relation": relation, "reference_values": references,
                 "direction": direction, "scale_description": scale, "comparison_note": note}
    # Percentile position within the panel. This is strictly more informative
    # than within/above/below range -- "inside the range" says nothing about
    # whether a molecule sits with the bulk or against the wall -- but it is
    # only meaningful once the panel is both large enough and complete within
    # its own definition. Before the panel became the full approved set it was
    # a curated subset of 36, where a percentile would have described the
    # curation rather than the population.
    if len(values) >= MIN_REFERENCES_FOR_PERCENTILE and _finite(value):
        ordered = sorted(values)
        benchmark["percentile"] = _percentile_of(value, ordered)
        benchmark["quartiles"] = {
            "p5": _quantile(ordered, 0.05), "p25": _quantile(ordered, 0.25),
            "p75": _quantile(ordered, 0.75), "p95": _quantile(ordered, 0.95)}
        benchmark["percentile_basis"] = (
            f"Position among {len(values)} approved human or humanised Fv pairs scored "
            f"under this protocol. A percentile is a position, not a verdict.")
    else:
        benchmark["percentile"] = None
        if values:
            benchmark["percentile_basis"] = (
                f"Withheld: {len(values)} references is below the {MIN_REFERENCES_FOR_PERCENTILE} "
                f"needed for a percentile to carry more information than the observed range.")
    return benchmark


def _quantile(ordered, fraction):
    """Linear-interpolated quantile of an already-sorted list."""
    if not ordered:
        return None
    if len(ordered) == 1:
        return ordered[0]
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def _anti_signature(result):
    provenance = result.get("provenance", {})
    artifacts = provenance.get("checkpoint_artifacts", [])
    hashes = {a.get("filename"): a.get("sha256") for a in artifacts}
    if not provenance.get("version") or not hashes.get("model.safetensors") or not hashes.get("config.json"):
        return None
    return provenance["version"], hashes["model.safetensors"], hashes["config.json"]


def annotate_antiberty(results: dict, chains: dict | None = None) -> dict:
    """Return enriched per-chain copies; pass analysis['chains'] for H/K/L type.

    AntiBERTy Heavy/Light classification is not used to infer kappa/lambda.
    """
    output = deepcopy(results)
    panel = _load_panel()
    chains = chains or {}
    for label, result in output.items():
        chain_type = chains.get(label, {}).get("numbering", {}).get("chain_type")
        signature = _anti_signature(result)
        references = []
        for reference in panel.get("references", []):
            if chain_type not in {"H", "K", "L"}:
                continue
            anti = reference.get("antiberty", {})
            if anti.get("status") != "complete":
                continue
            for ref_label, item in anti.get("chains", {}).items():
                if reference.get("chain_types", {}).get(ref_label) == chain_type and signature is not None and _anti_signature(item) == signature:
                    value = item.get("pseudo_log_likelihood")
                    if _finite(value):
                        references.append({"name": reference.get("name", reference.get("id", "Reference")), "value": value})
        type_name = {"H": "heavy", "K": "kappa", "L": "lambda"}.get(chain_type, "unknown subclass")
        note = f"Same {type_name} chain type, AntiBERTy package version and checkpoint/config hashes required."
        if not references:
            note += " No compatible scored controls available; no benchmark range is inferred."
        pll = result.get("pseudo_log_likelihood")
        pll_scale = "Mean masked observed-residue log-probability (natural log); at most 0. Higher/less negative means greater compatibility with this model, not greater clinical safety."
        benchmark = _benchmark("pseudo_log_likelihood", "Mean pseudo log-likelihood", pll, "nats/residue", references,
                               "Higher means greater sequence-model compatibility", pll_scale, note)
        perplexity = result.get("pseudo_perplexity")
        if perplexity is None and _finite(pll) and -700 < pll <= 0:
            perplexity = math.exp(-pll)
        perplexity_refs = [{"name": r["name"], "value": math.exp(-r["value"])} for r in references if -700 < r["value"] <= 0]
        perplexity_benchmark = _benchmark("pseudo_perplexity", "Pseudo-perplexity", perplexity, "", perplexity_refs,
                                         "Lower means greater sequence-model compatibility",
                                         "exp(-mean pseudo log-likelihood), minimum 1; a model score, not a PTM or immunogenicity percentage.", note)
        relation = benchmark["query_relation"]
        if relation == "unavailable":
            conclusion = "AntiBERTy describes sequence-model compatibility; no matching local reference range is available for this chain."
        else:
            where = {"below_range": "lower than all", "above_range": "higher than all", "within_range": "within the range of the"}[relation]
            conclusion = f"Mean pseudo log-likelihood is {where} {benchmark['reference_n']} compatible {type_name} reference chain scores. This is model context, not a developability verdict."
        regions = [r for r in result.get("region_summaries", []) if _finite(r.get("mean_log_likelihood"))]
        findings = []
        if regions:
            lowest = min(regions, key=lambda r: r["mean_log_likelihood"])
            findings.append(f"{lowest['region']} has the lowest mean sequence-model compatibility in this chain ({lowest['mean_log_likelihood']:.3f} nats/residue). This does not establish a chemical liability or recommend mutation.")
        result["region_conclusion"] = findings[0] if findings else "Region summaries unavailable."
        if result.get("species"):
            findings.append(f"Top classifier labels: {result.get('species')}, {result.get('chain_class', 'unknown chain')}, {result.get('graft_class', 'unknown graft class')}. Their softmax scores are model classification outputs, not measured humanness.")
        result["conclusion"] = conclusion
        result["summary"] = {"headline": "Sequence-model context", "conclusion": conclusion,
                              "findings": findings, "limitations": [PANEL_BOUNDARY, "Reference overlap with AntiBERTy training data is unknown, so these controls are not an independent model-validation set.", "PLL is not a validated predictor of affinity, PTM rate, immunogenicity or clinical performance. Lowest-decile residue markers are within-chain rankings, not therapeutic-risk cutoffs."],
                              "benchmarks": [benchmark, perplexity_benchmark], "panel_note": _panel_note(panel), "chain_type": chain_type}
    return output


def _scale(key):
    if "patch_area" in key:
        return "Larger means more calculated patch area", "Area in Å² from this surface/APBS protocol; higher means a larger computed patch, not a validated aggregation, PK or viscosity threshold."
    if key in {"hyd_asa", "hph_asa", "aromatic_asa"}:
        return "Larger means more accessible area", "Accessible area in Å² of the stated atom/residue class. This is total accessible area, not necessarily one contiguous patch or a calibrated liability."
    if key.startswith("ann_index"):
        return "Below 1 indicates closer spacing than the randomized reference", "Observed/expected nearest-neighbor distance ratio; below 1 means closer spacing, above 1 more dispersed. The randomized null and finite sampling affect the estimate."
    if key.startswith("ripley_k"):
        return "Above 1 indicates more nearby pairs than the randomized reference", "Observed/expected local pair-count ratio at the specified distance; 1 is the randomized-null ratio. This is not a significance probability."
    if key == "fv_csp":
        return "Sign indicates same-sign versus opposite-sign chain charges", "VH charge multiplied by VL charge (e²); negative indicates opposite-sign net chain charges, positive same-sign. Magnitude is not a calibrated self-association measure."
    if key in {"vh_charge", "vl_charge", "fv_charge", "net_charge", "exposed_net_charge", "net_charge_cdr", "fv_chml", "exposed_fv_chml"}:
        return "Signed charge; no universally favorable direction", "Charge in elementary-charge units under the stated calculation/protonation conditions. Positive and negative indicate sign; chain-difference metrics are VH minus VL. Neither sign alone predicts PK or viscosity."
    if key == "dipole_moment":
        return "Larger means a larger computed dipole magnitude", "Debye; the adapter uses coordinates centered at the atomic centroid. Charge assignment and conformation affect this descriptor; there is no clinical cutoff."
    if key in {"cdr_length", "aromatic_cdr", "exposed_aromatic", "n_charged_res_fv"}:
        return "Larger means a greater residue count", "Integer residue count under the stated region/exposure definition. Neither a higher nor lower count is universally more developable."
    if key == "heiden_score":
        return "Protocol-specific surface hydrophobic potential", "Heiden hydrophobic potential descriptor from the upstream implementation. Compare only within the same surface preparation and computation protocol; no universal therapeutic cutoff."
    if key == "scm":
        return "Protocol-specific spatial charge descriptor", "Spatial charge map score from the upstream SCM implementation; not a viscosity measurement, probability or universally calibrated risk scale."
    return "No validated favorable direction assigned", "Upstream descriptor under the stated calculation protocol; no clinical acceptance range is inferred."


def _pm_compatible(query, reference):
    if not _finite(query.get("ph")) or not _finite(reference.get("ph")) or abs(query["ph"] - 7.4) > 1e-6 or abs(reference["ph"] - 7.4) > 1e-6:
        return False
    for key in ("source_commit", "adapter_version"):
        if not query.get(key) or query.get(key) != reference.get(key):
            return False
    for key in ("model_checkpoint", "protocol_version", "software_versions", "structure_model"):
        if query.get(key) != reference.get(key):
            return False
    return True


def _protocol_compatible(query_protocol, panel_protocol, feature_key):
    sequence_metric = feature_key in {"n_charged_res_fv", "vh_charge", "vl_charge", "fv_charge", "fv_csp"}
    software = ("anarcii", "numpy", "scipy") if sequence_metric else ("ImmuneBuilder", "anarcii", "freesasa", "openmm", "numpy", "scipy", "numba")
    if "patch_area" in feature_key or feature_key == "heiden_score":
        software += ("pdb2pqr", "apbs-binary")
        query_native=query_protocol.get('nanoshaper_override_sha256')
        reference_native=panel_protocol.get('nanoshaper_override_sha256')
        if (query_native or reference_native) and (not query_native or query_native=='unavailable' or query_native!=reference_native):
            return False
    for name in software:
        query_version = query_protocol.get("software_versions", {}).get(name)
        expected_version = panel_protocol.get("software_versions", {}).get(name)
        if not query_version or not expected_version or query_version != expected_version:
            return False
    if not sequence_metric:
        if query_protocol.get('partial_structure') or panel_protocol.get('partial_structure'):
            return False
        for key in ("structure_model", "preparation"):
            if not query_protocol.get(key) or query_protocol.get(key) != panel_protocol.get(key):
                return False
        def weights(protocol):
            return {item.get("filename"): item.get("sha256") for item in protocol.get("structure_checkpoint_artifacts", [])}
        query_weights, expected_weights = weights(query_protocol), weights(panel_protocol)
        if not query_weights or not expected_weights or query_weights != expected_weights or any(not value for value in query_weights.values()):
            return False
    return True


def annotate_propermab(result: dict, protocol: dict | None = None) -> dict:
    """Return an enriched copy; comparisons require matched pH 7.4 protocol."""
    output = deepcopy(result)
    panel = _load_panel()
    query_protocol = protocol or result.get("comparison_protocol", {})
    panel_protocol = panel.get("protocol", {})
    benchmarks, findings = [], []
    for feature in output.get("features", []):
        key, value = feature.get("key", ""), feature.get("value")
        references = []
        for reference in panel.get("references", []):
            run = reference.get("propermab", {})
            ref = run.get("results", {})
            if run.get("status") != "complete" or not _pm_compatible(output, ref) or not _protocol_compatible(query_protocol, ref.get("comparison_protocol") or panel_protocol, key):
                continue
            if ("patch_area" in key or key == "heiden_score") and (not output.get("patch_pipeline") or output.get("patch_pipeline") != ref.get("patch_pipeline") or output.get("patch_conditions") != ref.get("patch_conditions")):
                continue
            matched = next((item for item in ref.get("features", []) if item.get("key") == key and item.get("unit", "") == feature.get("unit", "")), None)
            if matched and _finite(matched.get("value")):
                references.append({"name": reference.get("name", reference.get("id", "Reference")), "value": matched["value"]})
        direction, scale = _scale(key)
        note = "Requires matching source commit, adapter, known protocol/checkpoint metadata and pH 7.4; surface patches additionally require the same native pipeline and conditions."
        if not references:
            note += " Compatible computed reference values are unavailable or protocol metadata is incomplete; no range is inferred."
        benchmark = _benchmark(key, feature.get("label", key), value, feature.get("unit", ""), references, direction, scale, note)
        anchor=next((r for r in references if r['name'].lower()=='trastuzumab'),None)
        if anchor and _finite(value):
            baseline=anchor['value']
            benchmark['reference_anchor']={'name':'Herceptin (trastuzumab)','value':baseline,
                'delta':value-baseline,'relative_difference_percent':100*(value-baseline)/abs(baseline) if baseline else None,
                'note':'Named therapeutic benchmark; numerical proximity is not a validated acceptance criterion.'}
        feature["benchmark"] = benchmark
        benchmarks.append(benchmark)
        if benchmark["query_relation"] == "above_range" and ("patch_area" in key or key in {"hyd_asa", "aromatic_asa"}):
            findings.append(f"{feature.get('label', key)} is larger than all {benchmark['reference_n']} compatible therapeutic references. Prioritize the corresponding surface/self-interaction assay assessment; this comparison does not establish a liability.")
    compared = sum(b["query_relation"] != "unavailable" for b in benchmarks)
    outside = sum(b["query_relation"] in {"above_range", "below_range"} for b in benchmarks)
    differences = [b for b in benchmarks if b["query_relation"] in {"above_range", "below_range"}]
    if differences:
        labels = [f"{b['label']} ({'above' if b['query_relation'] == 'above_range' else 'below'})" for b in differences[:8]]
        extra = f"; plus {len(differences)-8} more in the comparison table" if len(differences) > 8 else ""
        findings.insert(0, "Outside the observed reference span: " + "; ".join(labels) + extra + ". Direction alone does not establish a developability concern.")
    if compared:
        conclusion = f"{compared} descriptors have compatible local therapeutic reference comparisons; {outside} fall outside their observed reference ranges. These are descriptive differences, not validated failure calls."
    else:
        conclusion = "PROPERMAB descriptors are available with intrinsic scale explanations; no compatible local therapeutic reference range is available."
    if not output.get("features"):
        conclusion = "No PROPERMAB descriptors were available to interpret."
    if output.get("ph") != 7.4:
        findings.append("Reference comparisons are withheld because the query pH differs from the panel's pH 7.4 preparation protocol.")
    output["conclusion"] = conclusion
    output["summary"] = {"headline": "Developability descriptor context", "conclusion": conclusion,
                         "findings": findings, "limitations": [PANEL_BOUNDARY, "ANN and Ripley descriptors use 1,000 random samples in the upstream implementation; their values can vary between runs. Small differences from reference ranges should not be treated as robust outliers.", "These descriptors do not directly predict PK, concentration limits, viscosity, stability or CMC readiness. Formulation and measured assay data remain necessary."],
                         "benchmarks": benchmarks, "panel_note": _panel_note(panel), "published_context": _published_context(), "experimental_context": _experimental_context(panel, output)}
    return output


def _published_context():
    try:
        return json.loads((COMPARATOR_PATH.parent / "propermab_published_ranges.json").read_text())
    except (OSError, ValueError):
        return None


def _experimental_context(panel, query):
    try:
        data = json.loads((COMPARATOR_PATH.parent / "experimental_reference_candidates.json").read_text())
    except (OSError, ValueError):
        return None
    examples=[]
    for candidate in data['rows']:
        if not candidate.get('computed_example_eligible'):
            continue
        reference=next((r for r in panel.get('references',[]) if r.get('id')==candidate['name'].lower()),None)
        if not reference or reference.get('sequence_sha256')!=candidate.get('sequence_sha256'):
            continue
        run=reference.get('propermab',{});ref=run.get('results',{})
        if (run.get('status')!='complete' or not _pm_compatible(query,ref)
            or not query.get('patch_pipeline') or query.get('patch_pipeline')!=ref.get('patch_pipeline')
            or query.get('patch_conditions')!=ref.get('patch_conditions')):
            continue
        values={f['key']:f['value'] for f in ref.get('features',[]) if f['key'] in
                ('hyd_patch_area','hyd_patch_area_cdr','pos_patch_area','neg_patch_area')
                and _finite(f.get('value')) and _protocol_compatible(query.get('comparison_protocol',{}),ref.get('comparison_protocol',{}),f['key'])}
        if len(values)==4:examples.append({'name':candidate['name'],'values':values})
    data['computed_examples']=examples
    return data
