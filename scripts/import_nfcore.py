#!/usr/bin/env python3
"""Import the verified nf-core/cutandrun 3.2.2 duplicate-retained target output layout."""
import argparse
import csv
import json
from pathlib import Path
from cuttag import load, validate
from decisions import sha256


def convert(config, run, out):
    root, cfg, samples = load(config)
    report = validate(root, cfg, samples, True)
    if report['errors']:
        raise ValueError('; '.join(report['errors']))
    if cfg['pipeline']['name'] != 'nf-core/cutandrun' or cfg['pipeline']['version'] != '3.2.2':
        raise ValueError('This adapter supports nf-core/cutandrun 3.2.2 only')
    if json.loads((run/'status.json').read_text()).get('state') != 'COMPUTATIONAL_PASS':
        raise ValueError('Upstream run has not passed computational checks')
    original = json.loads((run/'project.snapshot.json').read_text())
    if any(cfg[key] != original[key] for key in ('pipeline', 'reference', 'nfcore_params', 'spikein')):
        raise ValueError('Configuration does not describe this upstream run')
    if json.loads((run/'samples.snapshot.json').read_text()) != samples:
        raise ValueError('Sample mapping differs from upstream: merge technical units and rerun upstream explicitly')
    params = json.loads((run/'params.json').read_text())
    if params.get('only_filtering') or params.get('dedup_target_reads') or params.get('remove_linear_duplicates') or params.get('peakcaller') != 'macs2' or not params.get('macs2_narrow_peak'):
        raise ValueError('Adapter requires production MACS2 narrow peaks and duplicate-retained target data')
    rows, seen = [], set()
    for sample in samples:
        sid = sample['sample_id']
        if sample['role'] != 'target' or sid in seen:
            continue
        seen.add(sid)
        prefix = sample['group']+'_R'+sample['replicate']
        bam = run/'output/02_alignment/bowtie2/target/markdup'/(prefix+'.target.markdup.sorted.bam')
        peak = run/'output/03_peak_calling/04_called_peaks/macs2'/(prefix+'.macs2_peaks.narrowPeak')
        if not bam.is_file() or not peak.is_file():
            raise ValueError('Expected published outputs absent for '+prefix+'; no fallback to a different BAM/peak policy')
        rows.append({'sample_id': sid, 'bam': str(bam.resolve()), 'peaks_bed': str(peak.resolve()), 'bam_policy': 'duplicates_retained'})
    if out.exists():
        raise ValueError('Output exists')
    out.mkdir(parents=True)
    with (out/'artifacts.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    (out/'provenance.json').write_text(json.dumps({'state': 'REVIEW_REQUIRED', 'adapter': 'nf-core/cutandrun 3.2.2 bowtie2 MACS2 narrow duplicate-retained',
        'config_sha256': sha256(config), 'run': str(run), 'params_sha256': sha256(run/'params.json'),
        'sources': [{**row, 'bam_sha256': sha256(Path(row['bam'])), 'peaks_sha256': sha256(Path(row['peaks_bed']))} for row in rows]}, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--config', required=True, type=Path)
    p.add_argument('--run', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    a = p.parse_args()
    convert(a.config.resolve(), a.run.resolve(), a.out.resolve())
