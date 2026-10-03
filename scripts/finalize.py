#!/usr/bin/env python3
"""Copy an inspectable portable delivery package, rewriting selected IGV resources."""
import argparse
import html
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tempfile
from urllib.parse import quote
import xml.etree.ElementTree as ET
from decisions import sha256


def safe_relative(value):
    path = PurePosixPath(value)
    if path.is_absolute() or '..' in path.parts or not path.parts or '\\' in value:
        raise ValueError('Unsafe package destination: ' + value)
    if path.parts[0] in ('manifest.json', 'index.html'):
        raise ValueError('Reserved package filename')
    return Path(path)


def package(manifest, out):
    if out.exists():
        raise ValueError('Output exists')
    spec = json.loads(manifest.read_text())
    if not spec.get('artifacts'):
        raise ValueError('No delivery artifacts')
    records, copies, destinations, sources = [], [], set(), {}
    for i, item in enumerate(spec['artifacts'], 1):
        source = (manifest.parent / item['path']).resolve()
        destination = safe_relative(item.get('destination', 'artifacts/%03d_%s' % (i, source.name)))
        if type(item.get('required', True)) is not bool:
            raise ValueError('required must be boolean')
        if not source.exists():
            if item.get('required', True):
                raise ValueError('Missing/broken artifact ' + str(source))
            records.append(dict(item, status='MISSING'))
            continue
        # Never recursively capture a package's own output or unrelated cache trees.
        if source.is_dir() and (out.resolve() == source or source in out.resolve().parents):
            raise ValueError('Package output cannot be inside a selected source directory')
        selected = sorted(source.rglob('*')) if source.is_dir() else [source]
        if any(path.is_symlink() for path in selected):
            raise ValueError('Directory artifacts contain symlinks; select resolved files explicitly')
        for path in selected:
            if path.is_dir():
                continue
            if not path.is_file():
                raise ValueError('Non-regular artifact')
            target = destination / path.relative_to(source) if source.is_dir() else destination
            if str(target) in destinations or path.resolve() in sources:
                raise ValueError('Duplicate source or destination in delivery manifest')
            destinations.add(str(target))
            sources[path.resolve()] = target
            copies.append((path, target, item.get('label', path.name)))
    for dest in destinations:
        if any(str(parent) in destinations for parent in Path(dest).parents):
            raise ValueError('Destination file/directory collision')
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.cuttag-package-', dir=out.parent) as temporary:
        staging = Path(temporary) / 'package'
        staging.mkdir()
        for source, destination, label in copies:
            target = staging / destination
            target.parent.mkdir(parents=True, exist_ok=True)
            source_hash = sha256(source)
            shutil.copy2(source, target)
            if sha256(target) != source_hash or sha256(source) != source_hash:
                raise ValueError('Artifact changed during copy')
            rewritten = False
            if source.name.endswith('.igv.xml'):
                tree = ET.parse(target)
                for element in tree.iter():
                    for key in ('path', 'genome'):
                        value = element.get(key)
                        if not value:
                            continue
                        if key == 'genome' and '/' not in value and not (source.parent / value).is_file():
                            continue  # named IGV genome identifier
                        if '://' in value:
                            raise ValueError('Remote IGV resource must be supplied as a local artifact: ' + value)
                        original = (source.parent / value).resolve()
                        if original not in sources:
                            raise ValueError('IGV resource missing from package specification: ' + value)
                        element.set(key, os.path.relpath(staging / sources[original], target.parent))
                        rewritten = True
                tree.write(target, encoding='UTF-8', xml_declaration=True)
            records.append({'path': destination.as_posix(), 'label': label, 'status': 'PRESENT',
                            'bytes': target.stat().st_size, 'sha256': sha256(target),
                            'source': str(source), 'source_sha256': source_hash, 'igv_paths_rewritten': rewritten})
        data = {'schema_version': 1, 'project': spec['project'], 'limitations': spec.get('limitations', []),
                'artifacts': records, 'state': 'REVIEW_REQUIRED'}
        (staging / 'manifest.json').write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
        body = '<h1>' + html.escape(spec['project']) + '</h1><p>Scientific release review required. Local files are copied; verify manifest hashes after transfer.</p>'
        body += '<ul>' + ''.join('<li>' + html.escape(x) + '</li>' for x in data['limitations']) + '</ul>'
        body += '<table><tr><th>Artifact</th><th>Status</th><th>SHA256</th></tr>'
        for record in records:
            label = html.escape(record.get('label', record['path']))
            if record['status'] == 'PRESENT':
                label = '<a href="' + quote(record['path']) + '">' + label + '</a>'
            body += '<tr><td>' + label + '</td><td>' + record['status'] + '</td><td>' + record.get('sha256', '') + '</td></tr>'
        (staging / 'index.html').write_text('<!doctype html><meta charset="utf-8"><title>CUTTag delivery</title><style>body{font:16px sans-serif;margin:40px}td,th{padding:8px;border:1px solid #ddd;word-break:break-all}table{border-collapse:collapse}</style>' + body + '</table>')
        staging.rename(out)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    package(a.manifest.resolve(), a.out.resolve())
