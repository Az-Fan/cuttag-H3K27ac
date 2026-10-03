#!/usr/bin/env python3
"""Verify package hashes after transfer without accessing original source paths."""
import argparse
import json
from pathlib import Path
from decisions import sha256
from finalize import safe_relative


def verify(package):
    manifest = json.loads((package / 'manifest.json').read_text())
    seen = set()
    for record in manifest['artifacts']:
        if record['status'] != 'PRESENT':
            continue
        relative = safe_relative(record['path'])
        path = package / relative
        if relative in seen or not path.is_file() or path.is_symlink() or package.resolve() not in path.resolve().parents:
            raise ValueError('Invalid package path: ' + str(relative))
        seen.add(relative)
        if path.stat().st_size != record['bytes'] or sha256(path) != record['sha256']:
            raise ValueError('Package content mismatch: ' + str(relative))
    return {'state': 'COMPUTATIONAL_PASS', 'files_verified': len(seen), 'scientific_review': manifest['state']}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--package', required=True, type=Path)
    a = p.parse_args()
    print(json.dumps(verify(a.package), indent=2))
