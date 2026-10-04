#!/usr/bin/env python3
"""Expand declared downstream sensitivities into auditable candidate workflows."""
import argparse
import hashlib
import itertools
import json
import re
import subprocess
import sys
from pathlib import Path

from plan_downstream import plan


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    temporary = Path(path).with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def source_sha256():
    root = Path(__file__).resolve().parent.parent
    h = hashlib.sha256()
    for path in sorted((root / 'scripts').glob('*.py')) + [root / 'pixi.lock']:
        h.update(path.name.encode())
        h.update(sha(path).encode())
    return h.hexdigest()


def expand(spec):
    sets = spec.get('artifact_sets')
    if not isinstance(sets, list) or not sets:
        raise ValueError('matrix requires nonempty artifact_sets')
    axes = spec.get('axes', {})
    defaults = {'mapq': [20], 'remove_duplicates': [False], 'fraction': [2/3],
                'universe': ['support_core'], 'blacklist_mode': ['subtract'], 'min_width': [50]}
    unknown = set(axes) - set(defaults)
    if unknown:
        raise ValueError('Unknown strategy axes: ' + ', '.join(sorted(unknown)))
    for key, values in axes.items():
        if not isinstance(values, list) or not values:
            raise ValueError('Axis must be a nonempty list: ' + key)
        defaults[key] = values
    candidates = []
    keys = list(defaults)
    for artifacts in sets:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', artifacts.get('id', '')):
            raise ValueError('Unsafe artifact set id')
        if not artifacts.get('manifest'):
            raise ValueError('Each artifact set requires manifest')
        for values in itertools.product(*(defaults[k] for k in keys)):
            params = dict(zip(keys, values))
            if type(params['remove_duplicates']) is not bool:
                raise ValueError('remove_duplicates axis requires booleans')
            if params['universe'] not in ('support_core', 'reproducible_union'):
                raise ValueError('Unsupported universe')
            if params['blacklist_mode'] not in ('subtract', 'drop'):
                raise ValueError('Unsupported blacklist_mode')
            if not isinstance(params['mapq'], int) or not 0 <= params['mapq'] <= 255:
                raise ValueError('Invalid MAPQ')
            if not isinstance(params['fraction'], (int, float)) or not 0 < params['fraction'] <= 1:
                raise ValueError('Invalid support fraction')
            if not isinstance(params['min_width'], int) or params['min_width'] < 1:
                raise ValueError('Invalid minimum width')
            name = artifacts['id'] + '_q%s_%s_f%s_%s_%s_w%s' % (
                params['mapq'], 'dedup' if params['remove_duplicates'] else 'keepdup',
                format(params['fraction'], '.4g').replace('.', 'p'), params['universe'],
                params['blacklist_mode'], params['min_width'])
            candidates.append({'candidate_id': name, 'artifact_set': artifacts['id'],
                               'artifacts_manifest': artifacts['manifest'], **params})
    names = [x['candidate_id'] for x in candidates]
    if len(names) != len(set(names)):
        raise ValueError('Candidate IDs collide; simplify or distinguish axis values')
    baseline = spec.get('baseline_candidate')
    if baseline is not None and baseline not in names:
        raise ValueError('baseline_candidate is not present in expanded matrix')
    return candidates


def build(config, matrix_path, out, max_candidates=48, execute=False, resume=False):
    config = Path(config).resolve()
    out = Path(out).resolve()
    if not config.is_file():
        raise ValueError('Missing project config: ' + str(config))
    matrix_path = matrix_path.resolve()
    spec = json.loads(matrix_path.read_text())
    candidates = expand(spec)
    if len(candidates) > max_candidates:
        raise ValueError('Expanded matrix has %d candidates; limit is %d' % (len(candidates), max_candidates))
    manifest_paths = {}
    for candidate in candidates:
        manifest = (matrix_path.parent / candidate['artifacts_manifest']).resolve()
        if not manifest.is_file():
            raise ValueError('Missing artifact manifest: ' + str(manifest))
        manifest_paths[candidate['candidate_id']] = manifest
    if resume:
        registry = out / 'strategy_matrix.json'
        if not out.is_dir() or not registry.is_file():
            raise ValueError('Resume needs an existing strategy_matrix.json')
        data = json.loads(registry.read_text())
        if data.get('config_sha256') != sha(config) or data.get('matrix_sha256') != sha(matrix_path):
            raise ValueError('Config or matrix changed since this run; use a new output directory')
        if data.get('template_source_sha256') != source_sha256():
            raise ValueError('Template code or environment lock changed since this run; use a new output directory')
        records = data['candidates']
        if [x['candidate_id'] for x in records] != [x['candidate_id'] for x in candidates]:
            raise ValueError('Expanded candidate list changed; use a new output directory')
    else:
        if out.exists():
            raise ValueError('Output exists; choose a new matrix directory or use --resume')
        out.mkdir(parents=True)
        records = []
        for candidate in candidates:
            manifest = manifest_paths[candidate['candidate_id']]
            candidate_out = out / 'candidates' / candidate['candidate_id']
            record = {**candidate, 'artifacts_manifest': str(manifest), 'artifacts_sha256': sha(manifest),
                      'plan_dir': str(candidate_out / 'plan'), 'run_dir': str(candidate_out / 'execution'),
                      'state': 'PLANNED'}
            records.append(record)
        data = {'schema_version': 1, 'state': 'REVIEW_REQUIRED', 'config': str(config.resolve()),
                'config_sha256': sha(config), 'matrix': str(matrix_path), 'matrix_sha256': sha(matrix_path),
                'template_source_sha256': source_sha256(),
                'baseline_candidate': spec.get('baseline_candidate'), 'candidate_count': len(records),
                'candidates': records,
                'interpretation': 'Candidates are sensitivity analyses. No scientific method is selected automatically.'}
        save(out / 'strategy_matrix.json', data)

    for index, candidate in enumerate(candidates):
        record = records[index]
        if record.get('state') == 'COMPUTATIONAL_PASS':
            if resume:
                candidate_out = out / 'candidates' / candidate['candidate_id']
                command = [sys.executable, str(Path(__file__).with_name('workflow.py')),
                           '--workflow', str(candidate_out / 'plan/workflow.json'),
                           '--out', record['run_dir'], '--resume']
                result = subprocess.run(command)
                if result.returncode:
                    record.update(state='OUTPUT_VALIDATION_FAILED', returncode=result.returncode)
                    data.update(execution_requested=execute, candidates=records)
                    save(out / 'strategy_matrix.json', data)
            continue
        manifest = Path(record['artifacts_manifest'])
        if not manifest.is_file() or sha(manifest) != record['artifacts_sha256']:
            record.update(state='INPUT_CHANGED', error='Artifact manifest missing or changed')
            save(out / 'strategy_matrix.json', data)
            continue
        candidate_out = out / 'candidates' / candidate['candidate_id']
        args = (config, manifest, candidate_out / 'plan', candidate['mapq'], candidate['fraction'],
                candidate['universe'], candidate['blacklist_mode'], candidate['min_width'], candidate['remove_duplicates'])
        if not (candidate_out / 'plan/workflow.json').is_file():
            if record.get('state') == 'PLAN_FAILED':
                continue
            try:
                plan(*args)
            except Exception as error:
                record.update(state='PLAN_FAILED', error=str(error))
                candidate_out.mkdir(parents=True, exist_ok=True)
                (candidate_out / 'plan_error.txt').write_text(str(error) + '\n')
        if record.get('state') == 'PLANNED' and (candidate_out / 'plan/workflow.json').is_file():
            record['workflow_sha256'] = sha(candidate_out / 'plan/workflow.json')
        if execute and record.get('state') == 'PLANNED':
            if sha(candidate_out / 'plan/workflow.json') != record['workflow_sha256']:
                record.update(state='PLAN_CHANGED', error='Workflow plan changed since registry creation')
            else:
                command = [sys.executable, str(Path(__file__).with_name('workflow.py')),
                           '--workflow', str(candidate_out / 'plan/workflow.json'), '--out', record['run_dir']]
                if Path(record['run_dir']).is_dir():
                    command.append('--resume')
                result = subprocess.run(command)
                record['returncode'] = result.returncode
                record['state'] = 'COMPUTATIONAL_PASS' if result.returncode == 0 else 'FAILED'
        data.update(execution_requested=execute, candidates=records)
        save(out / 'strategy_matrix.json', data)
    data.update(execution_requested=execute, candidates=records)
    save(out / 'strategy_matrix.json', data)
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', required=True, type=Path)
    p.add_argument('--matrix', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--max-candidates', type=int, default=48)
    p.add_argument('--execute', action='store_true', help='Run every independent candidate workflow sequentially')
    p.add_argument('--resume', action='store_true', help='Resume an unchanged matrix and candidate workflows')
    a = p.parse_args()
    if a.max_candidates < 1:
        p.error('--max-candidates must be positive')
    result = build(a.config, a.matrix, a.out, a.max_candidates, a.execute, a.resume)
    print(json.dumps({'state': result['state'], 'candidate_count': result['candidate_count'],
                      'execution_requested': result['execution_requested']}))


if __name__ == '__main__':
    main()
