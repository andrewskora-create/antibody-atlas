"""Local candidate comparison orchestration. Each structure belongs to its own input."""
from __future__ import annotations

import copy
import threading
import uuid
from typing import Literal

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field, model_validator

from .analysis import AA, clean_sequence

router = APIRouter()
variant_jobs: dict[str, dict] = {}
variant_lock = threading.Lock()
batch_lock = threading.Lock()


class Mutation(BaseModel):
    chain: Literal['H', 'L']
    position: int = Field(ge=1, le=600, strict=True)
    from_aa: str = Field(min_length=1, max_length=1)
    to_aa: str = Field(min_length=1, max_length=1)

    @model_validator(mode='after')
    def validate_residues(self):
        if self.from_aa not in AA or self.to_aa not in AA or self.from_aa == self.to_aa:
            raise ValueError('Use different uppercase standard amino acids for a substitution.')
        return self


class Variant(BaseModel):
    name: str = Field(min_length=1, max_length=80)
    mutations: list[Mutation] = Field(default_factory=list, max_length=30)
    heavy: str | None = Field(default=None, max_length=1200)
    light: str | None = Field(default=None, max_length=1200)

    @model_validator(mode='after')
    def validate_mode(self):
        if not self.name.strip():
            raise ValueError('Candidate name cannot be blank.')
        has_sequence = self.heavy is not None or self.light is not None
        if bool(self.mutations) == has_sequence:
            raise ValueError('Supply either mutations or candidate sequences, not both.')
        return self


class VariantRequest(BaseModel):
    analysis_id: str
    variants: list[Variant] = Field(min_length=1, max_length=8)
    isotype: Literal["fv", "IGHG1", "IGHG2", "IGHG3", "IGHG4"] = "fv"
    light_constant: Literal["IGKC", "IGLC2"] | None = None
    clip_c_terminal_lysine: bool = True
    run_structure: bool = False
    run_antiberty: bool = False
    run_propermab: bool = False


def candidate_sequences(baseline: dict, variant: Variant) -> tuple[dict, int]:
    original = {label: item['sequence'] for label, item in baseline['chains'].items()}
    if variant.mutations:
        edited = {label: list(seq) for label, seq in original.items()}
        seen = set()
        for mutation in variant.mutations:
            key = (mutation.chain, mutation.position)
            if key in seen:
                raise ValueError('Duplicate mutation position; provide each chain position once.')
            seen.add(key)
            if mutation.chain not in original or mutation.position > len(original[mutation.chain]):
                raise ValueError('Mutation position or chain is outside the analyzed input.')
            if original[mutation.chain][mutation.position - 1] != mutation.from_aa:
                raise ValueError('Mutation source residue does not match the baseline at the input position.')
            edited[mutation.chain][mutation.position - 1] = mutation.to_aa
        sequences = {label: ''.join(seq) for label, seq in edited.items()}
    else:
        sequences = {label: clean_sequence(raw) for label, raw in (('H', variant.heavy), ('L', variant.light))
                     if raw and raw.strip()}
    if set(sequences) != set(original):
        raise ValueError('Candidate sequences must include the same chains as the baseline.')
    if any(len(sequences[label]) != len(seq) for label, seq in original.items()):
        raise ValueError('This release compares substitutions only; insertions/deletions are not supported.')
    count = sum(a != b for label, seq in original.items() for a, b in zip(seq, sequences[label]))
    if count == 0:
        raise ValueError('Candidate is identical to the baseline.')
    return sequences, count


def _set(job_id, **changes):
    with variant_lock:
        variant_jobs[job_id].update(changes)


def _models(record, request, progress):
    from . import main as api
    seqs = {key: value['sequence'] for key, value in record['chains'].items()}
    states = {key: {'status': 'not_requested'} for key in ('structure', 'antiberty', 'propermab')}
    pdb_path = None
    if request.run_structure:
        progress('Predicting its own structure')
        jid = uuid.uuid4().hex
        with api.jobs_lock:
            api.jobs[jid] = {'id': jid, 'status': 'queued', 'analysis_id': record['analysis_id'],
                             'scheme': record['scheme'], 'ph': record['ph']}
        if 'H' not in seqs or any(not 90 <= len(seq) <= 150 for seq in seqs.values()):
            states['structure'] = {'status': 'unavailable', 'reason': 'Structure requires H and 90–150 residues per chain.'}
        else:
            api._fold(jid, seqs, record['ph'], record['scheme'])
            with api.jobs_lock:
                job = copy.deepcopy(api.jobs[jid])
            structure = job.get('structure_analysis')
            if job['status'] == 'complete' and structure:
                record['structure_analysis'] = structure
                record['provenance'].update(structure_model=job['method'],
                    structure_model_version=job['model_version'], structure_pdb_sha256=job['model_sha256'],
                    structure_input_sha256=record['provenance']['input_sha256'])
                pdb_path = api.RUNS / f'{jid}.pdb'
                states['structure'] = {'status': 'complete', 'job_id': jid, 'pdb_url': job['pdb_url']}
            else:
                states['structure'] = {'status': 'unavailable', 'reason': job.get('analysis_error') or job.get('error') or 'Structure unavailable.'}
    if request.run_antiberty:
        progress('Scoring both available chains with AntiBERTy')
        jid = uuid.uuid4().hex
        with api.jobs_lock:
            api.antiberty_jobs[jid] = {'id': jid, 'status': 'queued', 'analysis_id': record['analysis_id']}
        api._run_antiberty(jid, seqs, record['scheme'])
        with api.jobs_lock:
            job = copy.deepcopy(api.antiberty_jobs[jid])
        if job['status'] == 'complete':
            for label, result in job['results'].items():
                record['chains'][label]['antiberty'] = result
            record['provenance']['antiberty_model_version'] = job['model_version']
            states['antiberty'] = {'status': 'complete', 'job_id': jid}
        else:
            states['antiberty'] = {'status': 'unavailable', 'reason': job.get('error', 'AntiBERTy unavailable.')}
    if request.run_propermab:
        if not pdb_path:
            states['propermab'] = {'status': 'unavailable', 'reason': 'Own structure failed; structural PROPERMAB comparison withheld. No baseline structure reused.'}
        else:
            progress('Computing PROPERMAB on its own structure')
            jid = uuid.uuid4().hex
            with api.jobs_lock:
                api.propermab_jobs[jid] = {'id': jid, 'status': 'queued', 'analysis_id': record['analysis_id']}
            api._run_propermab(jid, record, pdb_path)
            with api.jobs_lock:
                job = copy.deepcopy(api.propermab_jobs[jid])
            if job['status'] == 'complete':
                record['propermab'] = job['results']
                from .model_selection import select_propermab_structure
                selected_repair = select_propermab_structure(record,job['results'])
                if selected_repair:
                    # Surface descriptors and the displayed coordinates must be the same model.
                    states['structure']['pdb_url'] = job['results']['terminal_repair'].get('pdb_url')
                    states['structure']['preparation'] = 'terminal_repaired'
                    if not states['structure']['pdb_url']:
                        states['structure'].update(status='unavailable', reason='Selected repaired coordinates are unavailable.')
                record['provenance'].update(propermab_input_sha256=record['provenance']['input_sha256'],
                    propermab_structure_pdb_sha256=record['provenance']['structure_pdb_sha256'])
                states['propermab'] = {'status': 'complete', 'job_id': jid}
            else:
                states['propermab'] = {'status': 'unavailable', 'reason': job.get('error', 'PROPERMAB unavailable.')}
    record['model_status'] = states
    return states


def _run(job_id, baseline, prepared, request):
    from . import main as api
    from .variant_comparison import compare_records
    with batch_lock:
        try:
            _set(job_id, status='running', progress='Preparing baseline')
            base_states = _models(baseline, request, lambda stage: _set(job_id, progress=f'Baseline: {stage}'))
            _set(job_id, baseline_models=base_states, baseline_structure_analysis=baseline.get("structure_analysis"))
            results = []
            for index, (variant, sequences, count) in enumerate(prepared, 1):
                prefix = f'Candidate {index}/{len(prepared)}'
                _set(job_id, progress=f'{prefix}: analyzing sequence')
                item = {'name': variant.name, 'mutation_count': count}
                try:
                    candidate = api.analyze(api.AnalyzeRequest(heavy=sequences.get('H', ''),
                        light=sequences.get('L', ''), scheme=baseline['scheme'], ph=baseline['ph']))
                    candidate['provenance']['variant_lineage'] = {'baseline_analysis_id': baseline['analysis_id'],
                        'baseline_input_sha256': baseline['provenance']['input_sha256'], 'mutation_count': count}
                    item['model_status'] = _models(candidate, request, lambda stage: _set(job_id, progress=f'{prefix}: {stage}'))
                    item.update(status='complete', analysis_id=candidate['analysis_id'], structure_analysis=candidate.get('structure_analysis'),
                                comparison=compare_records(baseline, candidate, isotype=request.isotype,
                                    light_constant=request.light_constant, clip_c_terminal_lysine=request.clip_c_terminal_lysine))
                    # Server records retain lineage. Profiler jobs are separate, as in the main app.
                    with api.jobs_lock:
                        api.analysis_records[candidate['analysis_id']]['provenance']['variant_lineage'] = candidate['provenance']['variant_lineage']
                except Exception as exc:
                    detail = exc.detail if isinstance(exc, HTTPException) else str(exc)
                    item.update(status='failed', error=str(detail)[:400])
                results.append(item)
                _set(job_id, results=copy.deepcopy(results))
            _set(job_id, status='complete', progress=f'Compared {sum(r["status"] == "complete" for r in results)}/{len(results)} candidates; unavailable results remain explicit.')
        except Exception as exc:
            _set(job_id, status='failed', error=str(exc)[:400], progress='Comparison stopped; no overall ranking assigned.')


@router.post('/api/variants')
def start_variants(request: VariantRequest):
    from . import main as api
    with api.jobs_lock:
        baseline = copy.deepcopy(api.analysis_records.get(request.analysis_id))
    if not baseline:
        raise HTTPException(404, 'Analyze the baseline again in this server session.')
    if len({v.name.strip() for v in request.variants}) != len(request.variants):
        raise HTTPException(422, 'Use distinct candidate names.')
    if request.run_propermab:
        from .eligibility import gate as propermab_gate
        blocked = propermab_gate()
        if blocked:
            raise HTTPException(403, blocked)
        if set(baseline['chains']) != {'H', 'L'} or not 4 <= baseline['ph'] <= 10:
            raise HTTPException(422, 'PROPERMAB requires paired VH/VL and pH 4–10.')
        request = request.model_copy(update={'run_structure': True})
    try:
        prepared = [(variant, *candidate_sequences(baseline, variant)) for variant in request.variants]
    except ValueError as exc:
        raise HTTPException(422, str(exc)) from exc
    # Preserve supplied order within mutation count: no score-derived ranking.
    prepared.sort(key=lambda item: item[2])
    jid = uuid.uuid4().hex
    state = {'id': jid, 'status': 'queued', 'progress': 'Queued for local comparison', 'results': [],
        'baseline_analysis_id': request.analysis_id, 'isotype': request.isotype,
        'note': 'Single substitutions run before combinations. Predicted deltas are not proof of improvement; no composite score or clinical pass/fail is assigned.'}
    with variant_lock:
        variant_jobs[jid] = state
    threading.Thread(target=_run, args=(jid, baseline, prepared, request), daemon=True).start()
    return copy.deepcopy(state)


@router.get('/api/variants/{job_id}')
def get_variants(job_id: str):
    with variant_lock:
        result = variant_jobs.get(job_id)
        if not result:
            raise HTTPException(404, 'Comparison job not found in this server session.')
        return copy.deepcopy(result)
