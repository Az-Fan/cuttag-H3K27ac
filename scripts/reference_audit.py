#!/usr/bin/env python3
"""Hash references and verify sequence dictionaries and BED/GTF coordinate bounds."""
import argparse
import json
from pathlib import Path
from decisions import sha256
from cuttag import load, resolve


def fasta_sizes(path):
    sizes, name = {}, None
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in sizes:
                    raise ValueError('Duplicate FASTA sequence: ' + name)
                sizes[name] = 0
            elif line.strip():
                if name is None:
                    raise ValueError('FASTA sequence before header')
                sequence = line.strip()
                if not set(sequence.upper()).issubset(set('ACGTRYSWKMBDHVN')):
                    raise ValueError('Invalid nucleotide FASTA sequence')
                sizes[name] += len(sequence)
    if not sizes or any(n == 0 for n in sizes.values()):
        raise ValueError('Empty reference sequence')
    return sizes


def audit(root, cfg):
    paths = {key: resolve(root, cfg['reference'][key]).resolve() for key in ('fasta', 'gtf', 'blacklist')}
    before = {key: sha256(path) for key, path in paths.items()}
    sizes = fasta_sizes(paths['fasta'])
    fai = Path(str(paths['fasta']) + '.fai')
    if fai.exists():
        rows = [line.split('\t') for line in fai.read_text().splitlines()]
        indexed = {r[0]: int(r[1]) for r in rows}
        if len(indexed) != len(rows) or indexed != sizes:
            raise ValueError('FASTA index sequence lengths/names differ from FASTA')
        paths['fai'] = fai
        before['fai'] = sha256(fai)
    counts = {}
    for key in ('gtf', 'blacklist'):
        n = 0
        with paths[key].open() as f:
            for line in f:
                if not line.strip() or line.startswith(('#', 'track ', 'browser ')):
                    continue
                fields = line.rstrip('\n').split('\t')
                if key == 'gtf':
                    if len(fields) != 9:
                        raise ValueError('GTF must have nine tab-separated columns')
                    chrom, start, end = fields[0], int(fields[3])-1, int(fields[4])
                else:
                    chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if chrom not in sizes or start < 0 or end <= start or end > sizes[chrom]:
                    raise ValueError(key + ' interval outside FASTA dictionary: ' + line[:100])
                n += 1
        if key == 'gtf' and not n:
            raise ValueError('Empty GTF')
        counts[key] = n
    mito = cfg['reference'].get('mito_name')
    if mito and mito not in sizes:
        raise ValueError('Configured mitochondrial sequence absent from FASTA')
    if cfg['spikein']['enabled']:
        paths['spikein_fasta'] = resolve(root, cfg['spikein']['fasta']).resolve()
        before['spikein_fasta'] = sha256(paths['spikein_fasta'])
        spike_sizes = fasta_sizes(paths['spikein_fasta'])
        if paths['spikein_fasta'] == paths['fasta'] or before['spikein_fasta'] == before['fasta']:
            raise ValueError('Spike-in FASTA is identical to target FASTA')
    else:
        spike_sizes = {}
    if before != {key: sha256(path) for key, path in paths.items()}:
        raise ValueError('Reference changed during audit')
    return {'state': 'COMPUTATIONAL_PASS', 'genome_label': cfg['reference']['genome'],
            'files': {key: {'path': str(path), 'sha256': before[key], 'bytes': path.stat().st_size} for key, path in paths.items()},
            'target_sizes': sizes, 'spikein_sizes': spike_sizes, 'interval_counts': counts,
            'limitations': ['Coordinate compatibility does not prove assembly identity or annotation release; confirm reference source/release in scientific review.']}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    root, cfg, _ = load(a.config)
    result = audit(root, cfg)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
