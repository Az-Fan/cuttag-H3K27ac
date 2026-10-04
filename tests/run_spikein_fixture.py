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
    (p / 'reads1.fq').write_text(''.join('@q%d/1\n%s\n+\n%s\n' % (i, ('ACGT' * 12 + 'AC') if i == 7 else 'C' * 50, 'I' * 50) for i in range(8)))
    (p / 'reads2.fq').write_text(''.join('@q%d/2\n%s\n+\n%s\n' % (i, ('TGCA' * 12 + 'TG') if i == 7 else 'G' * 50, 'I' * 50) for i in range(8)))
    def make_bam(name, accepted):
        sam = ['@HD\tVN:1.6\tSO:unsorted', '@SQ\tSN:lambda1\tLN:1000']
        for i in range(5):
            q = 60 if i < accepted else 5
            dup = i == 1
            for flag, pos, mate, tlen, mq in ((99, 101 + i * 100, 151 + i * 100, 100, q),
                                              (147, 151 + i * 100, 101 + i * 100, -100, q)):
                sam.append('\t'.join(map(str, [name + str(i), flag + (1024 if dup else 0), 'lambda1', pos, mq,
                                              '50M', '=', mate, tlen, 'A' * 50, 'I' * 50])))
        sam_path = p / (name + '.sam')
        bam_path = p / (name + '.bam')
        sam_path.write_text('\n'.join(sam) + '\n')
        subprocess.run(['samtools', 'view', '-b', '-o', str(bam_path), str(sam_path)], check=True)
        return bam_path

    make_bam('x', 4)
    make_bam('x2', 4)
    make_bam('y', 3)
    make_bam('z', 2)
    make_bam('w', 1)
    make_bam('igg', 4)
    with (p / 'manifest.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'biological_sample_id', 'unit_id', 'condition', 'group', 'role', 'calibration_group', 'bam'])
        w.writerow(['s1', 'bio1', 'u1', 'Control', 'Control', 'target', 'batchA', 'x.bam'])
        w.writerow(['s1', 'bio1', 'u2', 'Control', 'Control', 'target', 'batchA', 'x.bam'])
    # Same BAM cannot be reused as separate sequencing units.
    bad = subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_fragments.py'), '--manifest', str(p / 'manifest.tsv'),
                          '--reference', str(p / 'lambda.fa'), '--out', str(p / 'bad')], capture_output=True)
    assert bad.returncode != 0 and not (p / 'bad').exists()
    with (p / 'manifest.tsv').open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'biological_sample_id', 'unit_id', 'condition', 'group', 'role', 'calibration_group', 'bam'])
        w.writerow(['s1', 'bio1', 'u1', 'Control', 'Control', 'target', 'batchA', 'x.bam'])
        w.writerow(['s1', 'bio1', 'u2', 'Control', 'Control', 'target', 'batchA', 'x2.bam'])
        w.writerow(['s2', 'bio2', 'u3', 'KD', 'KD', 'target', 'batchA', 'y.bam'])
        w.writerow(['s3', 'bio3', 'u4', 'KD', 'KD', 'target', 'batchB', 'z.bam'])
        w.writerow(['s4', 'bio4', 'u6', 'Control', 'Control', 'target', 'batchB', 'w.bam'])
        w.writerow(['igg', 'bio_igg', 'u5', 'IgG', 'IgG', 'control', 'batchA', 'igg.bam'])
    subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_fragments.py'), '--manifest', str(p / 'manifest.tsv'),
                    '--reference', str(p / 'lambda.fa'), '--mapq', '20', '30', '--threads', '1', '--out', str(p / 'count')], check=True)
    got = []
    for q in (20, 30):
        with (p / 'count' / ('mapq%d_duplicates_retained' % q) / 'biological_sample_counts.tsv').open() as f:
            got.append(int(next(csv.DictReader(f, delimiter='\t'))['spikein_fragments']))
    assert got == [8, 8], got
    with (p / 'count/mapq20_duplicates_retained/biological_sample_counts.tsv').open() as f:
        counts = {r['sample_id']: r for r in csv.DictReader(f, delimiter='\t')}
    assert counts['s1']['calibration_group'] == counts['s2']['calibration_group'] == 'batchA'
    assert counts['s3']['calibration_group'] == 'batchB'
    assert int(counts['s1']['spikein_fragments']) == 8
    assert int(counts['s2']['spikein_fragments']) == 3
    assert int(counts['s3']['spikein_fragments']) == 2
    assert abs(float(counts['s1']['size_factor']) * float(counts['s2']['size_factor']) - 1) < 1e-9
    assert abs(float(counts['s3']['size_factor']) * float(counts['s4']['size_factor']) - 1) < 1e-9
    assert counts['igg']['calibration_status'] == 'excluded_non_target_role'
    assert counts['igg']['size_factor'] == ''
    with (p / 'count/mapq20_duplicates_retained/target_biological_sample_counts.tsv').open() as f:
        target_rows = list(csv.DictReader(f, delimiter='\t'))
    assert {r['sample_id'] for r in target_rows} == {'s1', 's2', 's3', 's4'}
    # The BAM/reference dictionary is part of the count contract.
    bad = subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_fragments.py'), '--manifest', str(p / 'manifest.tsv'),
                          '--reference', str(p / 'target.fa'), '--out', str(p / 'mismatch')], capture_output=True)
    assert bad.returncode != 0
    cross = p / 'cross.tsv'
    (p / 'reads1_unit2.fq').write_text((p / 'reads1.fq').read_text())
    (p / 'reads2_unit2.fq').write_text((p / 'reads2.fq').read_text())
    with cross.open('w') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'biological_sample_id', 'unit_id', 'calibration_group', 'fastq_1', 'fastq_2'])
        w.writerow(['s1', 'bio1', 'u1', 'batchA', 'reads1.fq', 'reads2.fq'])
        w.writerow(['s1', 'bio1', 'u2', 'batchA', 'reads1_unit2.fq', 'reads2_unit2.fq'])
    subprocess.run([sys.executable, str(ROOT / 'scripts/spikein_crossmap.py'), '--manifest', str(cross),
                    '--target-fasta', str(p / 'target.fa'), '--spikein-fasta', str(p / 'lambda.fa'),
                    '--max-pairs', '8', '--threads', '1', '--out', str(p / 'cross')], check=True)
    rows = json.loads((p / 'cross/provenance.json').read_text())['results']
    assert len(rows) == 2 and {r['unit_id'] for r in rows} == {'u1', 'u2'}
    assert all(r['sample_id'] == 's1' and r['biological_sample_id'] == 'bio1' for r in rows)
    for row in rows:
        assert row['pairs_tested'] == 8 and sum(row[k] for k in ('unique_target', 'unique_spikein', 'cross_or_discordant', 'ambiguous_or_low_mapq', 'unmapped_or_incomplete')) == 8
        assert row['unmapped_or_incomplete'] > 0, row
    print('PASS: MAPQ20/30 spike-in fragment grid, reference guard and unit-aware competitive paired mapping')
