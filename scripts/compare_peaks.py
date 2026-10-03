#!/usr/bin/env python3
"""Compare peak strategies by base coverage; diagnostic agreement is not scientific acceptance."""
import argparse
import csv
import itertools
import json
from pathlib import Path
from consensus import intervals, merge
from decisions import sha256


def overlap_bp(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        ca, sa, ea = a[i]
        cb, sb, eb = b[j]
        if ca == cb:
            total += max(0, min(ea, eb)-max(sa, sb))
        if ca < cb or (ca == cb and ea <= eb):
            i += 1
        else:
            j += 1
    return total


def compare(manifest, out):
    if out.exists():
        raise ValueError('Output exists')
    with manifest.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    grouped, summary = {}, []
    for row in rows:
        key = (row['sample_id'], row['strategy'])
        if key in grouped:
            raise ValueError('Duplicate sample/strategy')
        path = (manifest.parent/row['peaks_bed']).resolve()
        original = intervals(path)
        merged = merge(original)
        grouped[key] = merged
        summary.append(dict(row, peak_count=len(original), union_bp=sum(e-s for _, s, e in merged), sha256=sha256(path)))
    if not rows:
        raise ValueError('Empty manifest')
    pairs = []
    for sample in sorted({key[0] for key in grouped}):
        strategies = sorted(key[1] for key in grouped if key[0] == sample)
        for left, right in itertools.combinations(strategies, 2):
            a, b = grouped[(sample, left)], grouped[(sample, right)]
            na, nb = sum(e-s for _, s, e in a), sum(e-s for _, s, e in b)
            overlap = overlap_bp(a, b)
            pairs.append({'sample_id': sample, 'left': left, 'right': right, 'intersection_bp': overlap,
                          'union_bp': na+nb-overlap, 'jaccard_bp': overlap/(na+nb-overlap) if na+nb-overlap else None,
                          'left_bp_recovered': overlap/na if na else None, 'right_bp_recovered': overlap/nb if nb else None})
    out.mkdir(parents=True)
    (out/'comparison.json').write_text(json.dumps({'state': 'REVIEW_REQUIRED', 'strategies': summary, 'pairwise': pairs,
        'recommendation_rule': 'Assess background depth/structure, common-universe FRiP, biological replicate agreement and loci alongside coverage agreement. Do not maximize peak count or accept scale-large/no-control automatically.'}, indent=2)+'\n')
    return pairs


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    compare(a.manifest.resolve(), a.out.resolve())
