#!/usr/bin/env python3
"""Collect fragment-level QC without inventing unavailable duplicate/spike-in metrics."""
import argparse
import bisect
import collections
import csv
import json
from pathlib import Path
from consensus import intervals, merge
from decisions import sha256


def weighted_median(histogram):
    n = sum(histogram.values())
    positions = ((n - 1) // 2, n // 2)
    selected, cumulative = [], 0
    for length, count in sorted(histogram.items()):
        selected.extend(length for pos in positions if cumulative <= pos < cumulative + count)
        cumulative += count
    return sum(selected) / 2


def collect(manifest, out):
    if out.exists():
        raise ValueError('Output exists')
    with manifest.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    ids = [r['sample_id'] for r in rows]
    if not rows or any(not sid for sid in ids) or len(set(ids)) != len(ids):
        raise ValueError('One unique sample_id per biological sample required')
    metrics, provenance, histograms = [], [], []
    for row in rows:
        fragment = (manifest.parent / row['fragments_bed']).resolve()
        peak = (manifest.parent / row['peaks_bed']).resolve()
        before = {'fragments_sha256': sha256(fragment), 'peaks_sha256': sha256(peak)}
        peaks = intervals(peak)
        index = collections.defaultdict(lambda: ([], []))
        for chrom, start, end in merge(peaks):
            index[chrom][0].append(start)
            index[chrom][1].append(end)
        histogram = collections.Counter()
        count = in_peaks = 0
        with fragment.open() as f:
            for line in f:
                fields = line.split()
                if len(fields) < 3:
                    raise ValueError('Invalid fragment BED row')
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                if start < 0 or end <= start:
                    raise ValueError('Invalid fragment coordinates')
                count += 1
                histogram[end - start] += 1
                starts, ends = index[chrom]
                candidate = bisect.bisect_right(ends, start)
                if candidate < len(starts) and starts[candidate] < end:
                    in_peaks += 1  # each fragment counted once, including duplicate rows
        if not count:
            raise ValueError('No retained fragments for ' + row['sample_id'])
        if before != {'fragments_sha256': sha256(fragment), 'peaks_sha256': sha256(peak)}:
            raise ValueError('QC input changed while reading')
        metrics.append({'sample_id': row['sample_id'], 'target_fragments': count,
                        'fragments_in_peaks': in_peaks, 'frip': in_peaks / count,
                        'peak_count': len(peaks), 'peak_union_bp': sum(e-s for _, s, e in merge(peaks)),
                        'median_fragment_length': weighted_median(histogram),
                        'duplication_rate': 'NA', 'spikein_fraction': 'NA'})
        histograms.extend((row['sample_id'], length, n) for length, n in sorted(histogram.items()))
        provenance.append({'sample_id': row['sample_id'], 'fragments': str(fragment), 'peaks': str(peak), **before})
    out.mkdir(parents=True)
    with (out / 'metrics.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(metrics[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(metrics)
    with (out / 'length_histograms.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_id', 'length', 'fragments'])
        writer.writerows(histograms)
    (out / 'provenance.json').write_text(json.dumps({
        'state': 'REVIEW_REQUIRED', 'sources': provenance,
        'frip_definition': 'Retained paired fragments overlapping >=1 peak / all retained paired fragments. BED half-open. Duplicate rows retained.',
        'comparison_requirement': 'Use the same peak universe and fragment filtering policy for between-sample comparison.',
        'unavailable': {'duplication_rate': 'Cannot infer PCR duplicates from fragment coordinates.',
                        'spikein_fraction': 'Requires independently audited spike-in mapping and an explicit denominator.'}
    }, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    collect(a.manifest.resolve(), a.out)
