#!/usr/bin/env python3
"""Explicit DAG runner; resume only content-identical completed steps."""
import argparse
import datetime
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import shutil
from decisions import sha256


def save(path, data):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def fingerprint(paths, cwd):
    return {str((cwd / name).resolve()): sha256(cwd / name) for name in paths}


def validate(spec):
    import re
    seen = set()
    if not spec.get('steps'):
        raise ValueError('Workflow has no steps')
    for step in spec['steps']:
        sid = step['id']
        if not isinstance(sid, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', sid) or sid in seen:
            raise ValueError('Unsafe or duplicate step ID')
        if any(dep not in seen for dep in step.get('depends_on', [])):
            raise ValueError('Dependencies must precede their consumers; no unknown IDs or cycles')
        command = step.get('command')
        if not isinstance(command, list) or not command or not all(isinstance(x, str) and x for x in command):
            raise ValueError('Step command must be a nonempty argv array')
        for key in ('inputs', 'outputs'):
            if not isinstance(step.get(key, []), list) or not all(isinstance(x, str) and x for x in step.get(key, [])):
                raise ValueError(key + ' must be an array of file paths')
        if type(step.get('continue_on_error', False)) is not bool:
            raise ValueError('continue_on_error must be boolean')
        seen.add(sid)


def run(spec, out, resume):
    cwd = Path(spec.get('working_directory', Path.cwd())).resolve()
    if not cwd.is_dir():
        raise ValueError('Working directory absent')
    snapshot = out / 'workflow.json'
    if resume:
        if not snapshot.is_file() or json.loads(snapshot.read_text()) != spec:
            raise ValueError('Workflow changed: create a new run and new output paths')
    else:
        save(snapshot, spec)
    states, summary, stopped = {}, [], False
    for step in spec['steps']:
        sid = step['id']
        directory = out / sid
        directory.mkdir(exist_ok=True)
        previous = json.loads((directory / 'status.json').read_text()) if (directory / 'status.json').exists() else {}
        record = {'id': sid, 'state': 'BLOCKED'}
        unmet = [dep for dep in step.get('depends_on', []) if states[dep]['state'] != 'COMPUTATIONAL_PASS']
        if stopped or unmet:
            record.update(reason='upstream failure', unmet=unmet)
        else:
            try:
                inputs = fingerprint(step.get('inputs', []), cwd)
                executable = shutil.which(step['command'][0])
                if not executable:
                    raise ValueError('Executable not found: '+step['command'][0])
                payload = {'step': step, 'cwd': str(cwd), 'inputs': inputs, 'executable_sha256': sha256(Path(executable)),
                           'dependencies': {dep: states[dep].get('signature') for dep in step.get('depends_on', [])}}
                signature = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
                record.update(signature=signature, inputs=inputs)
                if previous.get('state') == 'COMPUTATIONAL_PASS':
                    if not step.get('outputs') or previous.get('signature') != signature or previous.get('outputs') != fingerprint(step['outputs'], cwd):
                        raise ValueError('STALE: inputs/output changed or outputs undeclared; use new run/output paths')
                    record = dict(previous, resumed=True)
                else:
                    if any((cwd / name).exists() for name in step.get('outputs', [])):
                        raise ValueError('Output exists from incomplete/foreign run; use new output paths')
                    attempt = previous.get('attempt', 0) + 1
                    record.update(state='RUNNING', attempt=attempt, started=datetime.datetime.now(datetime.timezone.utc).isoformat())
                    save(directory / 'status.json', record)
                    save(directory / 'command.json', step['command'])
                    with (directory / ('execution.%d.log' % attempt)).open('w') as log:
                        result = subprocess.run(step['command'], cwd=cwd, stdout=log, stderr=subprocess.STDOUT)
                    record['returncode'] = result.returncode
                    if result.returncode:
                        raise ValueError('Command failed with exit code %d' % result.returncode)
                    missing = [name for name in step.get('outputs', []) if not (cwd / name).is_file()]
                    if missing:
                        record['missing_outputs'] = missing
                        raise ValueError('Required outputs absent')
                    record['outputs'] = fingerprint(step.get('outputs', []), cwd)
                    if inputs != fingerprint(step.get('inputs', []), cwd):
                        raise ValueError('Input changed during execution')
                    record.update(state='COMPUTATIONAL_PASS', finished=datetime.datetime.now(datetime.timezone.utc).isoformat())
            except Exception as error:
                record.update(state='FAILED', message=str(error))
            if record['state'] == 'FAILED' and not step.get('continue_on_error', False):
                stopped = True
        states[sid] = record
        summary.append(record)
        save(directory / 'status.json', record)
        save(out / 'workflow_status.json', {'state': 'RUNNING', 'steps': summary})
    passed = all(r['state'] == 'COMPUTATIONAL_PASS' for r in summary)
    save(out / 'workflow_status.json', {'state': 'COMPUTATIONAL_PASS' if passed else 'REVIEW_REQUIRED', 'steps': summary})
    return 0 if passed else 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--workflow', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    spec = json.loads(a.workflow.read_text())
    validate(spec)
    if a.resume and not a.out.is_dir():
        p.error('Resume run does not exist')
    if not a.resume:
        a.out.mkdir(parents=True, exist_ok=False)
    with (a.out / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return run(spec, a.out, a.resume)


if __name__ == '__main__':
    raise SystemExit(main())
