#!/usr/bin/env python3
"""Summarize executed strategy branches on one fixed pooled peak evaluation set."""
import argparse
import csv
import itertools
import json
import math
import statistics
import subprocess
from pathlib import Path

from collect_qc import collect
from consensus import intervals, merge
from compare_peaks import overlap_bp
from decisions import sha256


def read_tsv(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_bed(path, rows):
    with Path(path).open('w') as f:
        for chrom, start, end in rows:
            f.write('%s\t%d\t%d\n' % (chrom, start, end))


def correlation(xs, ys):
    if len(xs) < 2:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    dx, dy = [x-mx for x in xs], [y-my for y in ys]
    den = math.sqrt(sum(x*x for x in dx) * sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy)) / den if den else None


def replicate_correlations(plan_dir):
    counts = read_tsv(plan_dir / 'counts.tsv')
    metadata = read_tsv(plan_dir / 'metadata.tsv')
    by_condition = {}
    for row in metadata:
        by_condition.setdefault(row['condition'], []).append(row['sample_id'])
    output = []
    for condition, samples in sorted(by_condition.items()):
        for left, right in itertools.combinations(sorted(samples), 2):
            xs = [math.log1p(float(row[left])) for row in counts]
            ys = [math.log1p(float(row[right])) for row in counts]
            output.append({'condition': condition, 'left_sample': left, 'right_sample': right,
                           'tested_peaks': len(counts), 'pearson_log1p_raw_counts': correlation(xs, ys)})
    values = [x['pearson_log1p_raw_counts'] for x in output if x['pearson_log1p_raw_counts'] is not None]
    return output, statistics.mean(values) if values else None


def count_common_universe(fragments_table, peaks_bed, out_tsv):
    """Count every candidate's fragments over the same merged evaluation intervals."""
    source = read_tsv(fragments_table)
    peaks = merge(intervals(peaks_bed))
    if not source or not peaks:
        raise ValueError('Common-universe counting needs samples and evaluation peaks')
    columns = []
    peak_file = peaks_bed.resolve()
    for item in source:
        fragment = Path(item['fragments_bed'])
        fragment = fragment if fragment.is_absolute() else fragments_table.parent / fragment
        command = ['bedtools', 'coverage', '-counts', '-a', str(peak_file), '-b', str(fragment.resolve())]
        result = subprocess.run(command, text=True, capture_output=True, check=True)
        rows = [line.split('\t') for line in result.stdout.splitlines()]
        if [(r[0], int(r[1]), int(r[2])) for r in rows] != peaks:
            raise ValueError('Common-universe coverage coordinates/order changed')
        columns.append((item['sample_id'], [int(r[-1]) for r in rows]))
    with out_tsv.open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['peak_id', 'chrom', 'start', 'end'] + [x[0] for x in columns])
        for i, (chrom, start, end) in enumerate(peaks, 1):
            writer.writerow(['evaluation_%08d' % i, chrom, start, end] + [x[1][i-1] for x in columns])
    return out_tsv


def correlations_from_counts(counts_path, metadata_path):
    counts = read_tsv(counts_path)
    metadata = read_tsv(metadata_path)
    by_condition = {}
    for row in metadata:
        by_condition.setdefault(row['condition'], []).append(row['sample_id'])
    output = []
    for condition, samples in sorted(by_condition.items()):
        for left, right in itertools.combinations(sorted(samples), 2):
            xs = [math.log1p(float(row[left])) for row in counts]
            ys = [math.log1p(float(row[right])) for row in counts]
            output.append({'condition': condition, 'left_sample': left, 'right_sample': right,
                           'tested_peaks': len(counts), 'pearson_log1p_raw_counts': correlation(xs, ys)})
    values = [x['pearson_log1p_raw_counts'] for x in output if x['pearson_log1p_raw_counts'] is not None]
    return output, statistics.mean(values) if values else None


def verify_run(status):
    if status.get('state') != 'COMPUTATIONAL_PASS':
        return status.get('state', 'UNKNOWN')
    for step in status.get('steps', []):
        if step.get('state') != 'COMPUTATIONAL_PASS':
            return 'STEP_STATUS_CHANGED'
        for filename, expected in step.get('outputs', {}).items():
            path = Path(filename)
            if not path.is_file() or sha256(path) != expected:
                return 'OUTPUT_HASH_MISMATCH'
    return 'COMPUTATIONAL_PASS'


def compare(matrix_path, out):
    matrix_path = Path(matrix_path).resolve()
    out = Path(out).resolve()
    if out.exists():
        raise ValueError('Output exists')
    matrix = json.loads(matrix_path.read_text())
    records = matrix['candidates']
    passed, summaries, all_peaks = [], [], []
    for row in records:
        plan_dir = Path(row['plan_dir'])
        status_path = Path(row['run_dir']) / 'workflow_status.json'
        state = verify_run(json.loads(status_path.read_text())) if status_path.is_file() else row['state']
        master = plan_dir / 'consensus/master.bed'
        fragments_table = plan_dir / 'fragments.tsv'
        summary = {'candidate_id': row['candidate_id'], 'state': state, 'mapq': row['mapq'],
                   'remove_duplicates': row['remove_duplicates'], 'fraction': row['fraction'],
                   'universe': row['universe'], 'blacklist_mode': row['blacklist_mode'],
                   'artifact_set': row['artifact_set']}
        if state != 'COMPUTATIONAL_PASS' or not master.is_file() or not fragments_table.is_file():
            summary['metrics_state'] = 'UNAVAILABLE'
            summaries.append(summary)
            continue
        peaks = merge(intervals(master))
        all_peaks.extend(peaks)
        widths = [e-s for _, s, e in peaks]
        summary.update(metrics_state='AVAILABLE', peak_count=len(peaks), median_width=statistics.median(widths) if widths else None,
                       union_bp=sum(widths), master_sha256=sha256(master))
        summary['within_condition_replicate_correlations'], summary['mean_replicate_correlation'] = replicate_correlations(plan_dir)
        passed.append((row, summary, peaks, fragments_table))
        summaries.append(summary)
    if not passed:
        raise ValueError('No completed candidate has valid peak and fragment outputs')

    out.mkdir(parents=True)
    evaluation_bed = out / 'pooled_candidate_peak_union.bed'
    pooled = merge(all_peaks)
    write_bed(evaluation_bed, pooled)
    sample_sets = []
    for row, _, _, fragments_table in passed:
        sample_sets.append({x['sample_id'] for x in read_tsv(fragments_table)})
    if any(samples != sample_sets[0] for samples in sample_sets[1:]):
        raise ValueError('Candidates have different biological sample sets; cannot compare FRiP on the same evaluation set')

    frip_by_candidate = {}
    for row, summary, _, fragments_table in passed:
        source = read_tsv(fragments_table)
        common_manifest = out / (row['candidate_id'] + '.evaluation_qc.tsv')
        with common_manifest.open('w') as f:
            writer = csv.DictWriter(f, fieldnames=['sample_id', 'fragments_bed', 'peaks_bed'], delimiter='\t')
            writer.writeheader()
            for item in source:
                writer.writerow({'sample_id': item['sample_id'], 'fragments_bed': str((fragments_table.parent / item['fragments_bed']).resolve()),
                                 'peaks_bed': str(evaluation_bed)})
        qc_out = out / (row['candidate_id'] + '.evaluation_qc')
        collect(common_manifest, qc_out)
        qc_rows = read_tsv(qc_out / 'metrics.tsv')
        values = {x['sample_id']: float(x['frip']) for x in qc_rows}
        frip_by_candidate[row['candidate_id']] = values
        summary['pooled_union_frip_mean'] = statistics.mean(values.values())
        summary['pooled_union_frip_by_sample'] = values
        common_counts = out / (row['candidate_id'] + '.common_universe_counts.tsv')
        count_common_universe(fragments_table, evaluation_bed, common_counts)
        summary['common_universe_within_condition_replicate_correlations'], summary['common_universe_mean_replicate_correlation'] = correlations_from_counts(
            common_counts, Path(row['plan_dir']) / 'metadata.tsv')
        summary['common_universe_counts'] = str(common_counts)
        summary['common_universe_counts_sha256'] = sha256(common_counts)

    pairwise = []
    for (left, left_summary, a, _), (right, right_summary, b, _) in itertools.combinations(passed, 2):
        na, nb = sum(e-s for _, s, e in a), sum(e-s for _, s, e in b)
        overlap = overlap_bp(a, b)
        pairwise.append({'left': left['candidate_id'], 'right': right['candidate_id'],
                         'intersection_bp': overlap, 'union_bp': na + nb - overlap,
                         'jaccard_bp': overlap / (na + nb - overlap) if na + nb - overlap else None,
                         'left_bp_recovered': overlap / na if na else None,
                         'right_bp_recovered': overlap / nb if nb else None})
    with (out / 'candidate_metrics.tsv').open('w') as f:
        keys = ['candidate_id', 'state', 'metrics_state', 'artifact_set', 'mapq', 'remove_duplicates',
                'fraction', 'universe', 'blacklist_mode', 'peak_count', 'median_width', 'union_bp',
                'mean_replicate_correlation', 'common_universe_mean_replicate_correlation', 'pooled_union_frip_mean']
        writer = csv.DictWriter(f, fieldnames=keys, delimiter='\t', extrasaction='ignore')
        writer.writeheader(); writer.writerows(summaries)
    baseline = matrix.get('baseline_candidate')
    data = {'state': 'REVIEW_REQUIRED', 'candidate_metrics': summaries, 'pairwise_peak_agreement': pairwise,
            'evaluation_peak_set': str(evaluation_bed), 'evaluation_peak_set_sha256': sha256(evaluation_bed),
            'baseline_candidate': baseline,
            'recommendation': 'No method is selected automatically. Review QC, replicate agreement, peak coverage and accepted differential-model stability. Pooled-union FRiP and common-universe replicate correlation use the same evaluation intervals for each candidate; candidate-specific correlation uses each candidate own peak set and is descriptive only for within-candidate QC.',
            'limitations': ['This report does not run DESeq2. Add reviewed model outputs before accepting a final strategy.',
                            'Peak overlap and pooled-union FRiP do not establish biological validity.',
                            'Candidates must contain the same biological samples to be compared.',
                            'Common-universe replicate correlation is comparable across candidates; candidate-specific correlation is computed on different regions and should not rank strategies.']}
    (out / 'comparison.json').write_text(json.dumps(data, indent=2) + '\n')
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    result = compare(a.matrix, a.out)
    print(json.dumps({'state': result['state'], 'candidates': len(result['candidate_metrics']),
                      'pairwise_comparisons': len(result['pairwise_peak_agreement'])}))


if __name__ == '__main__':
    main()
