#!/usr/bin/env python3
"""Plan a reproducible BAM/accepted-peak to QC/count/track workflow; stop before model review."""
import argparse
import csv
import json
from pathlib import Path
import sys
from cuttag import load, resolve, validate
from decisions import sha256


def write_table(path, columns, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def plan(config, artifacts, out, mapq, fraction, universe="support_core", blacklist_mode="subtract", min_width=50,
         remove_duplicates=False):
    root, cfg, samples = load(config)
    report = validate(root, cfg, samples, True)
    if report['errors']:
        raise ValueError('; '.join(report['errors']))
    if out.exists():
        raise ValueError('New plan directory required')
    targets = {r['sample_id']: r for r in samples if r['role'] == 'target'}
    with artifacts.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if len(rows) != len(targets) or {r['sample_id'] for r in rows} != set(targets):
        raise ValueError('Exactly one merged target BAM and peak file per declared biological sample required')
    for r in rows:
        if r.get('bam_policy') != 'duplicates_retained':
            raise ValueError('Baseline requires duplicate-retained BAM; use sensitivity runs for deduplicated data')
        for key in ('bam', 'peaks_bed'):
            r[key] = str((artifacts.parent / r[key]).resolve())
            if not Path(r[key]).is_file():
                raise ValueError('Missing ' + key)
    fasta = resolve(root, cfg['reference']['fasta'])
    sizes = Path(str(fasta) + '.fai')
    if not sizes.is_file():
        raise ValueError('FASTA .fai required for tracks')
    bam_hashes=[sha256(Path(r['bam'])) for r in rows]
    if len(set(bam_hashes))!=len(bam_hashes):
        raise ValueError('Identical BAM content assigned to multiple statistical samples')
    out.mkdir(parents=True)
    scripts = Path(__file__).resolve().parent
    steps, peaks, fragments, metrics, metadata = [], [], [], [], []
    lock = root / 'pixi.lock'
    def add(sid, script, args, inputs, outputs, depends=()):
        steps.append({'id': sid, 'depends_on': list(depends), 'command': [sys.executable, str(scripts / script)] + list(map(str, args)),
                      'inputs': list(map(str, [config, artifacts, lock] + sorted(scripts.glob('*.py')) + inputs)), 'outputs': list(map(str, outputs))})
    add('reference', 'reference_audit.py', ['--config', config, '--out', out/'reference.json'],
        [fasta, sizes, resolve(root, cfg['reference']['gtf']), resolve(root, cfg['reference']['blacklist'])], [out/'reference.json'])
    for row in rows:
        sid = row['sample_id']
        sample = targets[sid]
        dest = out/'fragments'/sid
        fragment_args = ['--bam', row['bam'], '--mapq', mapq, '--threads', cfg['resources']['cpus'], '--out', dest]
        if remove_duplicates:
            fragment_args.append('--remove-duplicates')
        add('fragment_'+sid, 'fragments.py', fragment_args,
            [row['bam']], [dest/'fragments.bed', dest/'summary.json', dest/'length_histogram.tsv'], ['reference'])
        peaks.append({'group': sample['condition'], 'biological_sample_id': sample['biological_sample_id'], 'peaks_bed': row['peaks_bed']})
        fragments.append({'sample_id': sid, 'fragments_bed': str(dest/'fragments.bed')})
        metrics.append({'sample_id': sid, 'fragments_bed': str(dest/'fragments.bed'), 'peaks_bed': str(out/'consensus/master.bed')})
        metadata.append({'sample_id': sid, 'biological_sample_id': sample['biological_sample_id'], 'condition': sample['condition'],
                         'replicates_confirmed': str(cfg['analysis']['replicates_confirmed']).upper()})
    for name, items in [('peaks.tsv', peaks), ('fragments.tsv', fragments), ('qc.tsv', metrics), ('metadata.tsv', metadata)]:
        write_table(out/name, list(items[0]), items)
    blacklist = resolve(root, cfg['reference']['blacklist'])
    add('consensus', 'consensus.py', ['--manifest', out/'peaks.tsv', '--fraction', fraction, '--blacklist', blacklist, '--universe', universe, '--blacklist-mode', blacklist_mode, '--min-width', min_width, '--out', out/'consensus'],
        [out/'peaks.tsv', blacklist] + [r['peaks_bed'] for r in rows], [out/'consensus/master.bed', out/'consensus/provenance.json', out/'consensus/interval_decisions.json'], ['reference'])
    dependencies = ['fragment_'+r['sample_id'] for r in rows] + ['consensus']
    fragment_files = [r['fragments_bed'] for r in fragments]
    add('counts', 'count_fragments.py', ['--samples', out/'fragments.tsv', '--peaks', out/'consensus/master.bed', '--out', out/'counts.tsv'],
        [out/'fragments.tsv', out/'consensus/master.bed'] + fragment_files, [out/'counts.tsv'], dependencies)
    add('metrics', 'collect_qc.py', ['--manifest', out/'qc.tsv', '--out', out/'qc_metrics'],
        [out/'qc.tsv', out/'consensus/master.bed', scripts/'consensus.py'] + fragment_files, [out/'qc_metrics/metrics.tsv', out/'qc_metrics/provenance.json'], dependencies)
    add('qc_atlas', 'qc_atlas.py', ['--metrics', out/'qc_metrics/metrics.tsv', '--out', out/'qc_atlas'],
        [out/'qc_metrics/metrics.tsv'], [out/'qc_atlas/index.html', out/'qc_atlas/qc_atlas.pdf', out/'qc_atlas/qc_atlas.png'], ['metrics'])
    add('tracks', 'tracks.py', ['--manifest', out/'fragments.tsv', '--sizes', sizes, '--genome', cfg['reference']['genome'], '--out', out/'tracks'],
        [out/'fragments.tsv', sizes] + fragment_files,
        [out/'tracks'/(r['sample_id']+'.CPM'+extension) for r in rows for extension in ('.bw', '.igv.xml', '.json')], dependencies)
    spec = {'working_directory': str(root), 'steps': steps}
    (out/'workflow.json').write_text(json.dumps(spec, indent=2)+'\n')
    (out/'handoff.json').write_text(json.dumps({'state': 'REVIEW_REQUIRED', 'config_sha256': sha256(config),
        'formal_inference_eligible': report['formal_inference_eligible'], 'warnings': report['warnings'],
        'next': 'Review upstream/peak/independence decisions, then prepare model review and run_differential.py. No formal model is launched by this plan.',
        'fragment_policy': {'mapq_both_mates': mapq, 'duplicates': 'removed' if remove_duplicates else 'retained', 'proper_pairs': True},
        'consensus_support_fraction': fraction,'universe':universe,'blacklist_mode':blacklist_mode,'min_width':min_width}, indent=2)+'\n')
    return out/'workflow.json'


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True, type=Path)
    p.add_argument('--artifacts', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--mapq', type=int, default=20)
    p.add_argument('--fraction', type=float, default=2/3)
    p.add_argument('--universe',choices=['support_core','reproducible_union'],default='support_core')
    p.add_argument('--blacklist-mode',choices=['subtract','drop'],default='subtract')
    p.add_argument('--min-width',type=int,default=50);p.add_argument('--remove-duplicates',action='store_true')
    a = p.parse_args()
    if not 0 <= a.mapq <= 255 or not 0 < a.fraction <= 1 or a.min_width<1:
        p.error('Invalid MAPQ/fraction')
    print(plan(a.config.resolve(), a.artifacts.resolve(), a.out.resolve(), a.mapq, a.fraction, a.universe, a.blacklist_mode, a.min_width, a.remove_duplicates))
