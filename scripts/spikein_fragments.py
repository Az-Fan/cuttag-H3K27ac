#!/usr/bin/env python3
"""Count accepted paired lambda fragments from lambda-only alignment BAMs.

The input manifest has one row per sequencing unit; technical units are summed
to their biological sample only after unit-level validation and counting.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import re
import subprocess
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fasta_sizes(path):
    sizes = {}
    name = None
    length = 0
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    sizes[name] = length
                name = line[1:].split()[0]
                if not name or name in sizes:
                    raise ValueError('Invalid/duplicate FASTA sequence name')
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        sizes[name] = length
    if not sizes or any(n <= 0 for n in sizes.values()):
        raise ValueError('Empty/invalid FASTA')
    return sizes


def read_manifest(path):
    with Path(path).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    required = {'sample_id', 'biological_sample_id', 'unit_id', 'condition', 'group', 'role', 'calibration_group', 'bam'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Manifest requires sample_id biological_sample_id unit_id condition group role calibration_group bam')
    units = set()
    sample_mapping = {}
    biological_mapping = {}
    sample_calibration_groups = {}
    paths = set()
    for r in rows:
        if not all(r[k].strip() for k in required):
            raise ValueError('Empty required manifest value')
        for k in ('sample_id', 'biological_sample_id', 'unit_id', 'condition', 'group'):
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', r[k]):
                raise ValueError('Unsafe identifier: ' + r[k])
        if r['role'] not in ('target', 'control'):
            raise ValueError('role must be target or control')
        if r['unit_id'] in units:
            raise ValueError('Duplicate sequencing unit: ' + r['unit_id'])
        units.add(r['unit_id'])
        calibration_group = r['calibration_group'].strip()
        if not calibration_group:
            raise ValueError('calibration_group is required; declare the experimentally comparable calibration set')
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', calibration_group):
            raise ValueError('Unsafe calibration_group: ' + calibration_group)
        identity = (r['biological_sample_id'], r['condition'], r['group'], r['role'])
        if r['sample_id'] in sample_mapping and sample_mapping[r['sample_id']] != identity:
            raise ValueError('Statistical sample maps to conflicting biological metadata')
        if r['biological_sample_id'] in biological_mapping and biological_mapping[r['biological_sample_id']] != (r['sample_id'], r['condition'], r['group'], r['role']):
            raise ValueError('Biological sample maps to multiple statistical samples or metadata')
        if r['sample_id'] in sample_calibration_groups and sample_calibration_groups[r['sample_id']] != calibration_group:
            raise ValueError('Statistical sample maps to multiple calibration groups')
        sample_mapping[r['sample_id']] = identity
        biological_mapping[r['biological_sample_id']] = (r['sample_id'], r['condition'], r['group'], r['role'])
        sample_calibration_groups[r['sample_id']] = calibration_group
        r['_calibration_group'] = calibration_group
        p = Path(r['bam'])
        p = p if p.is_absolute() else Path(path).resolve().parent / p
        p = p.resolve()
        if p in paths:
            raise ValueError('BAM reused across units')
        paths.add(p)
        if not p.is_file():
            raise ValueError('Missing BAM: ' + str(p))
        r['_bam'] = p
    return rows


def bam_dictionary(bam):
    text = subprocess.check_output(['samtools', 'view', '-H', str(bam)], text=True)
    result = {}
    for line in text.splitlines():
        if line.startswith('@SQ\t'):
            fields = dict(x.split(':', 1) for x in line.split('\t')[1:] if ':' in x)
            result[fields['SN']] = int(fields['LN'])
    return result


def count_one(row, reference, out, mapq, remove_duplicates, threads):
    bam = row['_bam']
    if bam_dictionary(bam) != reference:
        raise ValueError('BAM reference dictionary differs from spike-in FASTA for ' + row['unit_id'] +
                         '; use a spike-in-only alignment, not a mixed/target BAM')
    unit_dir = out / 'units' / row['unit_id']
    unit_dir.mkdir(parents=True)
    exclude = 4 | 8 | 256 | 512 | 2048 | (1024 if remove_duplicates else 0)
    # Collate groups QNAMEs with bounded memory and spill files as needed. MAPQ
    # remains a pair-level decision so one low-quality mate cannot orphan its mate.
    collate_cmd = ['samtools', 'collate', '-@', str(max(0, threads - 1)), '-T', str(unit_dir / 'collate_tmp'),
                   '-u', '-O', str(bam)]
    view_cmd = ['samtools', 'view', '-f', '2', '-F', str(exclude), '-']
    (unit_dir / 'commands.json').write_text(json.dumps({'commands': [collate_cmd, view_cmd],
        'pairing': 'samtools collate with bounded memory and temporary spill; streaming QNAME grouping',
        'MAPQ': 'both mates checked together'}, indent=2))
    collate_log = (unit_dir / 'collate.stderr.log').open('w')
    collate = subprocess.Popen(collate_cmd, text=True, stdout=subprocess.PIPE, stderr=collate_log)
    proc = subprocess.Popen(view_cmd, text=True, stdin=collate.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    collate.stdout.close()
    pair_total = collections.Counter()
    fragments = 0
    with (unit_dir / 'fragments.bed').open('w') as bed:
        def emit(records):
            nonlocal fragments
            if not records:
                return
            qname = records[0][0]
            if len(records) != 2:
                pair_total['invalid_pair_groups'] += 1
                return
            pair_total['proper_pair_groups'] += 1
            one = next((x for x in records if int(x[1]) & 64), None)
            two = next((x for x in records if int(x[1]) & 128), None)
            if one is None or two is None or one is two or one[2] != two[2] or one[2] == '*':
                pair_total['invalid_pair_groups'] += 1
                return
            if int(one[4]) < mapq or int(two[4]) < mapq:
                pair_total['low_mapq_pairs'] += 1
                return
            tlen = int(one[8])
            if tlen == 0 or int(two[8]) != -tlen:
                pair_total['invalid_template_length'] += 1
                return
            start = min(int(one[3]) - 1, int(two[3]) - 1)
            end = start + abs(tlen)
            if start < 0 or end <= start or end > reference[one[2]]:
                pair_total['invalid_template_length'] += 1
                return
            bed.write(f'{one[2]}\t{start}\t{end}\t{qname}\n')
            fragments += 1
        current_name = None
        records = []
        for line in proc.stdout:
            fields = line.rstrip('\n').split('\t')
            if current_name is not None and fields[0] != current_name:
                emit(records)
                records = []
            current_name = fields[0]
            records.append(fields)
            if len(records) > 2:
                raise ValueError('More than two primary alignments for QNAME '+fields[0])
        emit(records)
        stderr = proc.stderr.read()
        view_status = proc.wait()
        collate_status = collate.wait()
        collate_log.close()
        if view_status or collate_status:
            raise RuntimeError('samtools view failed: ' + stderr[-2000:])
    if pair_total['invalid_pair_groups'] or pair_total['invalid_template_length']:
        raise ValueError('Unexpected malformed proper-pair BAM records for ' + row['unit_id'])
    if fragments == 0:
        raise ValueError('No accepted spike-in fragments at MAPQ ' + str(mapq) + ' for ' + row['unit_id'])
    return {'sample_id': row['sample_id'], 'biological_sample_id': row['biological_sample_id'],
            'unit_id': row['unit_id'], 'condition': row['condition'], 'group': row['group'],
            'role': row['role'], 'calibration_group': row['_calibration_group'],
            'accepted_spikein_fragments': fragments, 'low_mapq_pairs': pair_total['low_mapq_pairs'],
            'pair_accounting': json.dumps(dict(pair_total), sort_keys=True),
            'mapq_both_mates': mapq,
            'duplicates_removed': remove_duplicates, 'bam': str(bam), 'bam_sha256': digest(bam),
            'reference_fasta_sha256': digest(reference_path), 'bam_bed_sha256': digest(unit_dir / 'fragments.bed'),
            'excluded_duplicate_records': 'flag 0x400 excluded' if remove_duplicates else 'retained',
            'unit_bed': str((unit_dir / 'fragments.bed').relative_to(out))}


def main():
    global reference_path
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True, type=Path)
    p.add_argument('--reference', required=True, type=Path, help='Exact FASTA used for lambda-only BAM alignment')
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--mapq', type=int, nargs='+', default=[20, 30])
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--remove-duplicates', action='store_true')
    p.add_argument('--equal-amount-confirmed', action='store_true', help='Evidence flag only; never auto-accepts calibration')
    p.add_argument('--added-at', default='unknown')
    p.add_argument('--calibration-scope', default='unknown')
    a = p.parse_args()
    if a.out.exists():
        p.error('Output exists')
    if a.threads < 1 or not a.mapq or len(set(a.mapq)) != len(a.mapq) or any(q < 0 or q > 255 for q in a.mapq):
        p.error('Invalid thread/MAPQ grid')
    a.reference = a.reference.resolve()
    if not a.reference.is_file():
        p.error('Missing reference FASTA')
    reference_path = a.reference
    expected = fasta_sizes(a.reference)
    rows = read_manifest(a.manifest.resolve())
    a.out.mkdir(parents=True)
    all_results = []
    for q in a.mapq:
        policy = f'mapq{q}' + ('_deduplicated' if a.remove_duplicates else '_duplicates_retained')
        policy_dir = a.out / policy
        policy_dir.mkdir()
        unit_rows = [count_one(r, expected, policy_dir, q, a.remove_duplicates, a.threads) for r in rows]
        # Aggregate technical sequencing units into the distinct biological sample only after counting.
        grouped = {}
        for r in unit_rows:
            z = grouped.setdefault(r['sample_id'], dict(sample_id=r['sample_id'],
                biological_sample_id=r['biological_sample_id'], condition=r['condition'], group=r['group'],
                role=r['role'], unit_count=0, spikein_fragments=0, calibration_group=r['calibration_group'],
                size_factor='', track_scale_relative='', calibration_accepted=False))
            if r['calibration_group'] != z['calibration_group']:
                raise ValueError('Statistical sample combines calibration groups')
            if r['biological_sample_id'] != z['biological_sample_id']:
                raise ValueError('Statistical sample combines distinct biological samples')
            z['unit_count'] += 1
            z['spikein_fragments'] += r['accepted_spikein_fragments']
        calibration_groups = {}
        for z in grouped.values():
            if z['role'] != 'target':
                z['calibration_status'] = 'excluded_non_target_role'
                continue
            calibration_groups.setdefault(z['calibration_group'], []).append(z)
        for name, group_rows in calibration_groups.items():
            if len(group_rows) < 2:
                raise ValueError('Calibration group needs at least two target biological samples: ' + name)
            if any(z['spikein_fragments'] <= 0 for z in group_rows):
                raise ValueError('Nonpositive biological-sample count in calibration group ' + name)
            gm = math.exp(sum(math.log(z['spikein_fragments']) for z in group_rows) / len(group_rows))
            for z in group_rows:
                z['size_factor'] = z['spikein_fragments'] / gm
                z['track_scale_relative'] = gm / z['spikein_fragments']
                z['calibration_accepted'] = False
                z['calibration_status'] = 'requires_scientific_review'
        biological_rows = list(grouped.values())
        for table, name in ((unit_rows, 'unit_counts.tsv'), (biological_rows, 'biological_sample_counts.tsv'),
                            ([z for z in biological_rows if z['role'] == 'target'], 'target_biological_sample_counts.tsv')):
            with (policy_dir / name).open('w') as f:
                w = csv.DictWriter(f, fieldnames=list(table[0] if table else biological_rows[0]), delimiter='\t')
                w.writeheader()
                w.writerows(table)
        all_results.append({'policy': policy, 'units': unit_rows, 'biological_samples': list(grouped.values())})
    (a.out / 'review.json').write_text(json.dumps({
        'accepted': False, 'method': 'lambda-only proper primary paired BAM; both mates MAPQ threshold; each pair counted once from TLEN',
        'reference_fasta': str(a.reference), 'reference_fasta_sha256': digest(a.reference),
        'reference_sequence_dictionary': expected, 'manifest_sha256': digest(a.manifest),
        'mapq_grid': a.mapq, 'duplicates_removed': a.remove_duplicates,
        'equal_amount_confirmed': a.equal_amount_confirmed, 'added_at': a.added_at,
        'calibration_scope': a.calibration_scope,
        'required_review': ['same input libraries and read units', 'unique mapping and MAPQ policy',
                            'cross-mapping against target genome', 'equal spike-in input and addition stage',
                            'biological-sample aggregation and calibration group', 'full data and target-library QC'],
        'policies': all_results}, indent=2))


if __name__ == '__main__':
    main()
