#!/usr/bin/env python3
"""Audit spike-in/target cross-mapping by aligning matched read pairs competitively."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
import subprocess
from pathlib import Path


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def fastq(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        while True:
            h = f.readline()
            if not h:
                return
            s, p, q = f.readline(), f.readline(), f.readline()
            if not h.startswith('@') or not p.startswith('+') or not s.strip() or len(s.strip()) != len(q.strip()):
                raise ValueError('Malformed FASTQ ' + str(path))
            yield h, s, p, q


def rid(record):
    return re.sub(r'/[12]$', '', record[0].split()[0])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True, type=Path, help='TSV: sample_id fastq_1 fastq_2')
    p.add_argument('--target-fasta', required=True, type=Path)
    p.add_argument('--spikein-fasta', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--max-pairs', type=int, default=20000, help='First N paired reads; 0 uses all reads')
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--mapq', type=int, default=20)
    a = p.parse_args()
    if a.out.exists():
        p.error('Output exists')
    if a.max_pairs < 0 or a.threads < 1 or not 0 <= a.mapq <= 255:
        p.error('Invalid pair limit, threads or MAPQ')
    rows = list(csv.DictReader(a.manifest.open(), delimiter='\t'))
    if not rows or not {'sample_id', 'fastq_1', 'fastq_2'}.issubset(rows[0]):
        p.error('Manifest requires sample_id fastq_1 fastq_2')
    if len({r['sample_id'] for r in rows}) != len(rows):
        p.error('Sample IDs must be unique; units should be listed separately')
    a.target_fasta = a.target_fasta.resolve(); a.spikein_fasta = a.spikein_fasta.resolve()
    a.out.mkdir(parents=True)
    combined = a.out / 'competitive.fa'
    seen = set()
    with combined.open('w') as out:
        for label, fasta_path in (('target', a.target_fasta), ('spikein', a.spikein_fasta)):
            with fasta_path.open() as f:
                for line in f:
                    if line.startswith('>'):
                        name = line[1:].split()[0]
                        prefixed = label + '__' + name
                        if prefixed in seen:
                            raise ValueError('Duplicate reference ID ' + prefixed)
                        seen.add(prefixed)
                        out.write('>' + prefixed + '\n')
                    else:
                        out.write(line)
    idx = a.out / 'index'; idx.mkdir()
    subprocess.run(['bowtie2-build', str(combined), str(idx / 'competitive')], check=True,
                   stdout=(a.out / 'index.log').open('w'), stderr=subprocess.STDOUT)
    results = []
    for row in rows:
        sid = row['sample_id']
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', sid):
            raise ValueError('Unsafe sample ID')
        r1 = Path(row['fastq_1']); r2 = Path(row['fastq_2'])
        r1 = r1 if r1.is_absolute() else a.manifest.resolve().parent / r1
        r2 = r2 if r2.is_absolute() else a.manifest.resolve().parent / r2
        dest = a.out / sid; dest.mkdir()
        fq1, fq2 = dest / 'tested_R1.fastq', dest / 'tested_R2.fastq'
        tested = 0
        pairs = itertools.zip_longest(fastq(r1), fastq(r2))
        if a.max_pairs:
            pairs = itertools.islice(pairs, a.max_pairs)
        with fq1.open('w') as x, fq2.open('w') as y:
            for one, two in pairs:
                if one is None or two is None or rid(one) != rid(two):
                    raise ValueError('FASTQ mates mismatch in ' + sid)
                x.writelines(one); y.writelines(two); tested += 1
        if not tested:
            raise ValueError('No reads in ' + sid)
        sam = dest / 'competitive.sam'
        cmd = ['bowtie2', '--end-to-end', '--very-sensitive', '--reorder', '--no-mixed', '--no-discordant',
               '-I', '10', '-X', '700', '-x', str(idx / 'competitive'), '-1', str(fq1), '-2', str(fq2),
               '-p', str(a.threads), '--seed', '42', '-S', str(sam)]
        with (dest / 'alignment.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        classes = {'unique_target': 0, 'unique_spikein': 0, 'cross_or_discordant': 0,
                   'ambiguous_or_low_mapq': 0, 'unmapped_or_incomplete': 0}
        current_name, pair = None, {}

        def classify(pair_rows):
            if 'r1' not in pair_rows or 'r2' not in pair_rows:
                classes['unmapped_or_incomplete'] += 1; return
            x, y = pair_rows['r1'], pair_rows['r2']
            if x[0] == '*' or y[0] == '*' or x[2] & 4 or y[2] & 4:
                classes['unmapped_or_incomplete'] += 1; return
            if min(x[1], y[1]) < a.mapq:
                classes['ambiguous_or_low_mapq'] += 1; return
            proper = bool(x[2] & 2) and bool(y[2] & 2)
            if proper and x[0].startswith('target__') and y[0].startswith('target__'):
                classes['unique_target'] += 1
            elif proper and x[0].startswith('spikein__') and y[0].startswith('spikein__'):
                classes['unique_spikein'] += 1
            else:
                classes['cross_or_discordant'] += 1

        with sam.open() as f:
            for line in f:
                if line.startswith('@'):
                    continue
                z = line.rstrip().split('\t'); flag = int(z[1])
                if flag & 64 and not flag & (256 | 2048):
                    mate, key = 'r1', 'r1'
                elif flag & 128 and not flag & (256 | 2048):
                    mate, key = 'r2', 'r2'
                else:
                    continue
                if current_name is not None and z[0] != current_name:
                    classify(pair)
                    pair = {}
                current_name = z[0]
                if key in pair:
                    raise ValueError('Multiple primary alignments for tested QNAME ' + z[0])
                pair[mate] = (z[2], int(z[4]), flag)
            if current_name is not None:
                classify(pair)
        if sum(classes.values()) != tested:
            raise ValueError('SAM pair accounting differs from tested reads for ' + sid)
        results.append({'sample_id': sid, 'pairs_tested': tested, 'mapq_both_mates': a.mapq,
                        **classes, 'target_fasta_sha256': sha(a.target_fasta),
                        'spikein_fasta_sha256': sha(a.spikein_fasta),
                        'source_R1_sha256': sha(r1), 'source_R2_sha256': sha(r2),
                        'tested_R1_sha256': sha(fq1), 'tested_R2_sha256': sha(fq2),
                        'command': cmd, 'interpretation': 'competitive mapping diagnostic; not a formal count'})
    with (a.out / 'crossmap_counts.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(results[0]), delimiter='\t'); w.writeheader(); w.writerows(results)
    (a.out / 'provenance.json').write_text(json.dumps({'max_pairs': a.max_pairs,
        'subset': 'first matched pairs per input; zero means full input', 'mapq_both_mates': a.mapq,
        'state': 'REVIEW_REQUIRED', 'results': results}, indent=2))


if __name__ == '__main__':
    main()
