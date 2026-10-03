#!/usr/bin/env python3
"""Create a clean project from allowlisted template files, excluding data and credentials."""
import argparse
import json
from pathlib import Path
import re
import shutil
from decisions import sha256


def initialize(source, destination, project_id):
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]*', project_id):
        raise ValueError('Invalid project ID')
    if destination.exists():
        raise ValueError('Destination exists; refusing to merge/overwrite')
    files = [source / name for name in ('README.md', 'PROJECT_STATUS.md', 'CHANGELOG.md', '.gitignore', 'pixi.toml', 'pixi.lock')]
    for directory, suffixes in [('scripts', {'.py', '.R'}), ('tests', {'.py'}), ('docs', {'.md', '.json', '.tsv'}), ('examples', {'.json', '.tsv', '.md'}), ('config', {'.json', '.tsv', '.config'}), ('.github/workflows', {'.yml', '.yaml'})]:
        files.extend(p for p in (source / directory).rglob('*') if p.is_file() and p.suffix in suffixes and '__pycache__' not in p.parts)
    records = []
    for path in files:
        if path.is_symlink():
            raise ValueError('Template source must not contain symlinks: ' + str(path))
        records.append({'path': str(path.relative_to(source)), 'sha256': sha256(path)})
    destination.mkdir(parents=True)
    for path in files:
        target = destination / path.relative_to(source)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    config = destination / 'config/project.json'
    cfg = json.loads(config.read_text())
    cfg['project_id'] = project_id
    cfg['qc'].update(upstream_accepted=False, peaks_accepted=False, official_test_run_id='', production_review='')
    cfg['analysis']['replicates_confirmed'] = False
    cfg['spikein'].update(equal_amount_confirmed=False, calibration_accepted=False, added_at='unknown', calibration_scope='')
    config.write_text(json.dumps(cfg, indent=2) + '\n')
    for name in ('data/raw', 'data/external', 'results', 'logs', 'work', 'shared_cache', 'reviews'):
        (destination / name).mkdir(parents=True, exist_ok=True)
    (destination / 'template_origin.json').write_text(json.dumps({'template_version': cfg['template_version'], 'source_files': records}, indent=2) + '\n')
    return destination


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--project-id', required=True)
    a = p.parse_args()
    print(initialize(Path(__file__).resolve().parents[1], a.out.resolve(), a.project_id))
