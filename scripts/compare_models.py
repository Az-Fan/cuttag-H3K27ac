#!/usr/bin/env python3
"""Compare reviewed differential models by reciprocal-overlap matched regions."""
import argparse
import bisect
import csv
import itertools
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def read_results(path):
    with Path(path).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    required = {'peak_id', 'chrom', 'start', 'end', 'padj', 'log2FoldChange', 'direction'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('DE result needs peak coordinates, padj, log2FoldChange and direction')
    out = []
    seen = set()
    for row in rows:
        chrom, start, end = row['chrom'], int(row['start']), int(row['end'])
        if not chrom or start < 0 or end <= start or row['peak_id'] in seen:
            raise ValueError('Invalid or duplicate peak in ' + str(path))
        seen.add(row['peak_id'])
        try:
            lfc = float(row['log2FoldChange'])
        except (TypeError, ValueError):
            lfc = math.nan
        try:
            padj = float(row['padj'])
        except (TypeError, ValueError):
            padj = math.nan
        out.append({'peak_id': row['peak_id'], 'chrom': chrom, 'start': start, 'end': end,
                    'lfc': lfc, 'padj': padj, 'direction': row['direction']})
    out.sort(key=lambda z: (z['chrom'], z['start'], z['end'], z['peak_id']))
    return out


def validate_review(path, normalization):
    review_path = Path(path)
    review = json.loads(review_path.read_text())
    if not all(isinstance(review.get(k), str) and review[k].strip() for k in ('reviewer', 'reviewed_at', 'reason')):
        raise ValueError('Each model requires a completed scientific review record')
    if review.get('model', {}).get('normalization') != normalization:
        raise ValueError('Review normalization differs from comparison manifest')
    required = ['upstream_accepted', 'peaks_accepted', 'replicates_confirmed', 'simple_design_accepted']
    if normalization == 'spikein':
        required += ['spikein_counting_accepted', 'spikein_calibration_accepted']
    if any(review.get('decisions', {}).get(key) is not True for key in required):
        raise ValueError('Model review has not accepted all required decisions')
    model = review.get('model', {})
    signature = tuple(model.get(key) for key in ('numerator', 'denominator', 'design', 'alpha', 'abs_log2fc', 'min_total_count'))
    signature += (review.get('inputs', {}).get('metadata_sha256'),)
    if any(value is None for value in signature):
        raise ValueError('Review is missing sample identity, contrast or filtering thresholds')
    return signature


def match_regions(left, right, reciprocal_overlap):
    right_by_chrom = defaultdict(list)
    for j, row in enumerate(right):
        right_by_chrom[row['chrom']].append((row['start'], row['end'], j))
    starts = {chrom: [x[0] for x in items] for chrom, items in right_by_chrom.items()}
    possible = []
    for i, a in enumerate(left):
        items = right_by_chrom.get(a['chrom'], [])
        if not items:
            continue
        stop = bisect.bisect_left(starts[a['chrom']], a['end'])
        for bstart, bend, j in items[:stop]:
            if bend <= a['start']:
                continue
            overlap = min(a['end'], bend) - max(a['start'], bstart)
            if overlap <= 0 or overlap / (a['end'] - a['start']) < reciprocal_overlap:
                continue
            if overlap / (bend - bstart) < reciprocal_overlap:
                continue
            possible.append((overlap, min(overlap/(a['end']-a['start']), overlap/(bend-bstart)), i, j))
    used_left, used_right, matches = set(), set(), []
    for overlap, reciprocal, i, j in sorted(possible, reverse=True):
        if i in used_left or j in used_right:
            continue
        used_left.add(i); used_right.add(j); matches.append((left[i], right[j]))
    return matches, len(used_left), len(used_right)


def pearson(xs, ys):
    if len(xs) < 2:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    dx, dy = [x-mx for x in xs], [y-my for y in ys]
    den = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy)) / den if den else None


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2
        for k in range(i, j):
            result[order[k]] = rank
        i = j
    return result


def compare(manifest, out, reciprocal_overlap=0.5):
    if out.exists():
        raise ValueError('Output exists')
    with Path(manifest).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if len(rows) < 2 or not {'candidate_id', 'normalization', 'result_tsv', 'review_json'}.issubset(rows[0]):
        raise ValueError('Manifest columns: candidate_id normalization result_tsv review_json')
    if any(row['normalization'] not in ('conventional', 'spikein') for row in rows):
        raise ValueError('normalization must be conventional or spikein')
    ids = [r['candidate_id'] for r in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('Candidate IDs must be unique')
    loaded = []
    signatures = []
    for row in rows:
        result_path = Path(row['result_tsv']); result_path = result_path if result_path.is_absolute() else Path(manifest).resolve().parent / result_path
        review_path = Path(row['review_json']); review_path = review_path if review_path.is_absolute() else Path(manifest).resolve().parent / review_path
        signatures.append(validate_review(review_path, row['normalization']))
        loaded.append((row, read_results(result_path)))
    if len(set(signatures)) != 1:
        raise ValueError('Candidate models must use the same contrast, design and filtering thresholds')
    candidate_summaries = []
    for row, result_rows in loaded:
        candidate_summaries.append({'candidate_id': row['candidate_id'], 'normalization': row['normalization'],
            'tested_regions': len(result_rows), 'gain_regions': sum(x['direction'] == 'gain' for x in result_rows),
            'loss_regions': sum(x['direction'] == 'loss' for x in result_rows),
            'not_significant_regions': sum(x['direction'] == 'not_significant' for x in result_rows),
            'missing_or_unclassified_regions': sum(x['direction'] not in ('gain', 'loss', 'not_significant') for x in result_rows)})
    pairs = []
    for (left_meta, left), (right_meta, right) in itertools.combinations(loaded, 2):
        matches, matched_left, matched_right = match_regions(left, right, reciprocal_overlap)
        finite = [(a, b) for a, b in matches if math.isfinite(a['lfc']) and math.isfinite(b['lfc'])]
        xs, ys = [a['lfc'] for a, _ in finite], [b['lfc'] for _, b in finite]
        sig = [(a, b) for a, b in finite if a['direction'] in ('gain', 'loss') and b['direction'] in ('gain', 'loss')]
        any_sig = [(a, b) for a, b in finite if (a['direction'] in ('gain', 'loss')) != (b['direction'] in ('gain', 'loss'))]
        same_sig = sum(a['direction'] == b['direction'] for a, b in sig)
        pairs.append({'left': left_meta['candidate_id'], 'right': right_meta['candidate_id'],
                      'left_normalization': left_meta['normalization'], 'right_normalization': right_meta['normalization'],
                      'left_tested_regions': len(left), 'right_tested_regions': len(right), 'matched_regions': len(matches),
                      'left_match_fraction': matched_left / len(left) if left else None,
                      'right_match_fraction': matched_right / len(right) if right else None,
                      'lfc_pearson': pearson(xs, ys), 'lfc_spearman': pearson(ranks(xs), ranks(ys)) if len(xs) >= 2 else None,
                      'matched_log2fc_sign_concordance': sum((a['lfc'] >= 0) == (b['lfc'] >= 0) for a, b in finite) / len(finite) if finite else None,
                      'both_significant_regions': len(sig),
                      'significance_status_disagreements': len(any_sig),
                      'gain_loss_direction_concordance': same_sig / len(sig) if sig else None,
                      'gain_loss_direction_flips': len(sig) - same_sig})
    out.mkdir(parents=True)
    (out / 'comparison.json').write_text(json.dumps({'state': 'REVIEW_REQUIRED', 'reciprocal_overlap_threshold': reciprocal_overlap,
        'candidate_summaries': candidate_summaries, 'pairwise': pairs,
        'recommendation': 'No normalization is selected automatically. Review calibration evidence and matched-region effect stability.',
        'limitations': ['One-to-one reciprocal-overlap matching is coordinate based; it is not a validated enhancer-to-gene mapping.',
                        'Different tested universes can leave many unmatched regions.',
                        'This report requires completed review records but does not independently rerun the model input hash audit.']}, indent=2) + '\n')
    return pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--reciprocal-overlap', type=float, default=0.5)
    a = p.parse_args()
    if not 0 < a.reciprocal_overlap <= 1:
        p.error('reciprocal overlap must be in (0,1]')
    print(json.dumps({'pairwise': len(compare(a.manifest.resolve(), a.out.resolve(), a.reciprocal_overlap)), 'state': 'REVIEW_REQUIRED'}))


if __name__ == '__main__':
    main()
