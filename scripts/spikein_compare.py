#!/usr/bin/env python3
"""Controlled Bowtie2 geometry comparison on identical paired reads; no auto calibration."""
import argparse
import csv
import gzip
import itertools
import json
from pathlib import Path
import re
import subprocess
import sys
from decisions import sha256

MODES = {
    'end_to_end_no_overlap': ['--end-to-end', '--very-sensitive', '--no-overlap', '--no-dovetail'],
    'end_to_end_overlap': ['--end-to-end', '--very-sensitive', '--no-dovetail'],
    'end_to_end_dovetail': ['--end-to-end', '--very-sensitive', '--dovetail'],
    'local_overlap': ['--local', '--very-sensitive-local', '--no-dovetail'],
}


def records(path):
    with gzip.open(path, 'rt') as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            if not header.startswith('@') or not plus.startswith('+') or not seq.strip() or len(seq.strip()) != len(qual.strip()):
                raise ValueError('Invalid FASTQ ' + str(path))
            yield header, seq, plus, qual


def compare(manifest, fasta, out, max_pairs, threads, modes):
    if out.exists():
        raise ValueError('Output exists')
    with manifest.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    ids = [r['sample_id'] for r in rows]
    if not rows or len(ids) != len(set(ids)) or any(not re.fullmatch('[A-Za-z][A-Za-z0-9_]*', sid) for sid in ids):
        raise ValueError('Unique safe sample IDs required; concatenate technical units explicitly before comparison')
    out.mkdir(parents=True)
    index = out / 'index'
    index.mkdir()
    build = ['bowtie2-build', str(fasta.resolve()), str(index / 'spikein')]
    with (out / 'index.log').open('w') as log:
        subprocess.run(build, stdout=log, stderr=subprocess.STDOUT, check=True)
    provenance = {'fasta_sha256': sha256(fasta), 'build_command': build, 'max_pairs': max_pairs,
                  'subset_method': 'first N matched pairs per sample; 0 means full input',
                  'limitations': 'Subset sensitivity audit cannot establish full-data calibration or experimental spike-in scope.',
                  'bowtie2_version': subprocess.check_output(['bowtie2', '--version'], text=True).splitlines()[0],
                  'calibration_accepted': False, 'samples': [], 'runs': []}
    log_rows = []
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    for row in rows:
        sid = row['sample_id']
        directory = out / sid
        directory.mkdir()
        sources = [(manifest.parent / row[k]).resolve() for k in ('fastq_1', 'fastq_2')]
        reads = [directory / 'reads_1.fastq', directory / 'reads_2.fastq']
        n = 0
        pairs = itertools.zip_longest(records(sources[0]), records(sources[1]))
        if max_pairs:
            pairs = itertools.islice(pairs, max_pairs)
        with reads[0].open('w') as one, reads[1].open('w') as two:
            for r1, r2 in pairs:
                if r1 is None or r2 is None or re.sub('/[12]$', '', r1[0].split()[0]) != re.sub('/[12]$', '', r2[0].split()[0]):
                    raise ValueError('FASTQ mates mismatch')
                one.writelines(r1)
                two.writelines(r2)
                n += 1
        if not n:
            raise ValueError('Empty FASTQ')
        hashes = [sha256(path) for path in reads]
        provenance['samples'].append({'sample_id': sid, 'source_paths': list(map(str, sources)), 'pairs_tested': n, 'tested_fastq_sha256': hashes})
        for mode in modes:
            sam, log = directory / (mode + '.sam'), directory / (mode + '.log')
            command = ['bowtie2', '-x', str(index / 'spikein'), '-1', str(reads[0]), '-2', str(reads[1]), '-S', str(sam),
                       '-p', str(threads), '--seed', '42', '--no-mixed', '--no-discordant', '--phred33', '-I', '10', '-X', '700'] + MODES[mode]
            with log.open('w') as f:
                result = subprocess.run(command, stdout=f, stderr=subprocess.STDOUT)
            record = {'sample_id': sid, 'mode': mode, 'command': command, 'returncode': result.returncode, 'tested_fastq_sha256': hashes}
            provenance['runs'].append(record)
            (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
            result.check_returncode()
            high_mapq = total = 0
            with sam.open() as f:
                for line in f:
                    if line.startswith('@'):
                        continue
                    fields = line.split('\t')
                    flag = int(fields[1])
                    if flag & 2 and flag & 64 and not flag & (256 | 2048):
                        total += 1
                        high_mapq += int(fields[4]) >= 20
            record.update(concordant_primary_pairs=total, read1_mapq20_pairs=high_mapq)
            if hashes != [sha256(path) for path in reads]:
                raise ValueError('Reads changed between diagnostic modes')
            log_rows.append((sid, mode, row.get('calibration_group','experiment'), str(log.resolve())))
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    with (out / 'logs.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample_id', 'mode', 'calibration_group', 'bowtie2_log'])
        writer.writerows(log_rows)
    subprocess.run([sys.executable, str(Path(__file__).with_name('spikein_audit.py')), '--manifest', str(out / 'logs.tsv'), '--out', str(out / 'audit')], check=True)
    (out / 'recommendation.md').write_text(
        '# Spike-in parameter review\n\n'
        'All modes used the same recorded paired reads. Compare concordant pair recovery, MAPQ and sample-relative factors in audit/counts.tsv and provenance.json.\n\n'
        'Use end-to-end overlap as the first candidate for short paired fragments. Retain no-overlap as the reference sensitivity run. Dovetail/local require additional alignment/insert/adapter review; higher mapping yield alone is not sufficient.\n\n'
        'No strategy is accepted automatically. Confirm spike-in identity, full-data counts, equal input, addition stage and interpretation scope before calibration. A subset is not a full-data count.\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--fasta', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--max-pairs', type=int, default=20000)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--modes', nargs='+', choices=list(MODES), default=list(MODES))
    a = p.parse_args()
    if a.max_pairs < 0 or a.threads < 1 or len(set(a.modes)) != len(a.modes):
        p.error('Invalid pair limit/threads or duplicate modes')
    compare(a.manifest.resolve(), a.fasta.resolve(), a.out.resolve(), a.max_pairs, a.threads, a.modes)
