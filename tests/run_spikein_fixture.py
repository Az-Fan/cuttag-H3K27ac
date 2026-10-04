"""Validate production spike-in pair counting and competitive mapping categories."""
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as d:
    p = Path(d)
    (p / 'lambda.fa').write_text('>lambda1\n' + 'A' * 1000 + '\n')
    (p / 'target.fa').write_text('>chr1\n' + 'C' * 1000 + '\n')
    (p / 'reads1.fq').write_text(''.join('@q%d/1\n%s\n+\n%s\n' % (i, 'C' * 50, 'I' * 50) for i in range(8)))
    (p / 'reads2.fq').write_text(''.join('@q%d/2\n%s\n+\n%s\n' % (i, 'G' * 50, 'I' * 50) for i in range(8)))
    sam = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:lambda1\tLN:1000']
    for i, q in enumerate((60, 60, 60, 60, 25)):
        dup = i == 1
        for flag, pos, mate, tlen, mq in ((99, 101 + i * 100, 151 + i * 100, 100, q),
                                          (147, 151 + i * 100, 101 + i * 100, -100, 0 if i == 4 else q)):
            sam.append('\t'.join(map(str, ['q%d' % i, flag + (1024 if dup else 0), 'lambda1', pos, mq,
                                          '50M', '=', mate, tlen, 'A' * 50, 'I' * 50])))
    (p / 'x.sam').write_text('\n'.join(sam) + '\n')
    subprocess.run(['samtools', 'view', '-b', '-o', str(p / 'x.bam'), str(p / 'x.sam')], check=True)
    with (p / 'manifest.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'biological_sample_id', 'unit_id', 'condition', 'group', 'role', 'bam'])
        w.writerow(['s1', 'bio1', 'u1', 'Control', 'Control', 'target', 'x.bam'])
        w.writerow(['s1', 'bio1', 'u2', 'Control', 'Control', 'target', 'x.bam'])
    # Same BAM cannot be reused as separate sequencing units.
    bad = subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_fragments.py'), '--manifest', str(p / 'manifest.tsv'),
                          '--reference', str(p / 'lambda.fa'), '--out', str(p / 'bad')], capture_output=True)
    assert bad.returncode != 0 and not (p / 'bad').exists()
    with (p / 'manifest.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'biological_sample_id', 'unit_id', 'condition', 'group', 'role', 'bam'])
        w.writerow(['s1', 'bio1', 'u1', 'Control', 'Control', 'target', 'x.bam'])
    subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_fragments.py'), '--manifest', str(p / 'manifest.tsv'),
                    '--reference', str(p / 'lambda.fa'), '--mapq', '20', '30', '--threads', '1', '--out', str(p / 'count')], check=True)
    got = []
    for q in (20, 30):
        with (p / 'count' / ('mapq%d_duplicates_retained' % q) / 'biological_sample_counts.tsv').open() as f:
            got.append(int(next(csv.DictReader(f, delimiter='\t'))['spikein_fragments']))
    assert got == [4, 4], got
    # The BAM/reference dictionary is part of the count contract.
    bad = subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_fragments.py'), '--manifest', str(p / 'manifest.tsv'),
                          '--reference', str(p / 'target.fa'), '--out', str(p / 'mismatch')], capture_output=True)
    assert bad.returncode != 0
    cross = p / 'cross.tsv'
    with cross.open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'fastq_1', 'fastq_2']); w.writerow(['s1', 'reads1.fq', 'reads2.fq'])
    subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_crossmap.py'), '--manifest', str(cross),
                    '--target-fasta', str(p / 'target.fa'), '--spikein-fasta', str(p / 'lambda.fa'),
                    '--max-pairs', '8', '--threads', '1', '--out', str(p / 'cross')], check=True)
    row = json.loads((p / 'cross/provenance.json').read_text())['results'][0]
    assert row['pairs_tested'] == 8 and sum(row[k] for k in ('unique_target', 'unique_spikein', 'cross_or_discordant', 'ambiguous_or_low_mapq', 'unmapped_or_incomplete')) == 8
    print('PASS: MAPQ20/30 spike-in fragment grid, reference guard and competitive paired mapping')
