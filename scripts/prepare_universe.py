#!/usr/bin/env python3
"""Prepare ORA and motif inputs from the complete adjusted-test universe and annotation."""
import argparse
import csv
import json
import math
from pathlib import Path
from decisions import sha256


def finite(value):
    try:
        return math.isfinite(float(value))
    except (ValueError, TypeError):
        return False


def prepare(results, annotation, out, alpha, lfc):
    if out.exists():
        raise ValueError('Output exists')
    with results.open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    with annotation.open() as f:
        annotations = list(csv.DictReader(f, delimiter='\t'))
    if not rows or len({r['peak_id'] for r in rows}) != len(rows):
        raise ValueError('Complete differential table must have unique peak IDs')
    if len({r['peak_id'] for r in annotations}) != len(annotations):
        raise ValueError('Annotation must have one explicit gene association per peak')
    mapping = {r['peak_id']: r for r in annotations}
    eligible, excluded = [], []
    for row in rows:
        pid = row['peak_id']
        if pid not in mapping:
            raise ValueError('Missing annotation for ' + pid)
        ann = mapping[pid]
        chrom = ann.get('seqnames', ann.get('chrom'))
        if chrom != row['chrom'] or int(ann['start']) != int(row['start']) or int(ann['end']) != int(row['end']):
            raise ValueError('Annotation coordinates do not match result ' + pid)
        if not finite(row.get('padj')) or not finite(row.get('log2FoldChange')):
            excluded.append({'peak_id': pid, 'reason': 'not eligible for adjusted significance test'})
            continue
        padj, fold = float(row['padj']), float(row['log2FoldChange'])
        if not 0 <= padj <= 1:
            raise ValueError('Invalid adjusted p value')
        direction = 'gain' if padj < alpha and fold >= lfc and fold > 0 else 'loss' if padj < alpha and fold <= -lfc and fold < 0 else 'not_significant'
        gene = ann.get('geneId', '')
        if not gene.isdigit():
            gene = None
            excluded.append({'peak_id': pid, 'reason': 'no single Entrez geneId; excluded from gene lists only'})
        region = ann.get('annotation', '')
        region = 'promoter' if region.startswith('Promoter') else 'distal' if region and region != 'NA' else 'unclassified'
        eligible.append(dict(row, direction=direction, gene=gene, region=region))
    if not eligible:
        raise ValueError('No peaks have finite adjusted p value and effect size')
    out.mkdir(parents=True)
    summaries = {}
    for region in ('all', 'promoter', 'distal'):
        subset = [row for row in eligible if region == 'all' or row['region'] == region]
        directory = out / region
        directory.mkdir()
        genes = {row['gene'] for row in subset if row['gene']}
        (directory / 'universe.txt').write_text(''.join(g + '\n' for g in sorted(genes)))
        directions = {}
        for direction in ('gain', 'loss'):
            selected = [row for row in subset if row['direction'] == direction]
            target_genes = {row['gene'] for row in selected if row['gene']}
            directions[direction] = target_genes
            (directory / (direction + '.genes.txt')).write_text(''.join(g + '\n' for g in sorted(target_genes)))
            for label, peaks in [('target', selected), ('background', [row for row in subset if row['direction'] != direction])]:
                (directory / (direction + '.' + label + '.bed')).write_text(''.join('%s\t%s\t%s\t%s\n' % (r['chrom'], r['start'], r['end'], r['peak_id']) for r in peaks))
        summaries[region] = {'eligible_peaks': len(subset), 'universe_genes': len(genes),
                             'gain_genes': len(directions['gain']), 'loss_genes': len(directions['loss']),
                             'genes_with_both_gain_and_loss_peaks': sorted(directions['gain'] & directions['loss'])}
    (out / 'provenance.json').write_text(json.dumps({'results_sha256': sha256(results), 'annotation_sha256': sha256(annotation),
        'alpha': alpha, 'abs_log2fc': lfc, 'universe_definition': 'Peaks with finite padj and log2FoldChange; genes deduplicated within each region.',
        'gene_id_type': 'Entrez', 'regions': summaries, 'exclusions': excluded, 'state': 'REVIEW_REQUIRED',
        'limitations': ['Nearest gene is a positional association; distal is non-promoter annotation, not experimentally established enhancer linkage.',
                        'Genes may have both gain and loss peaks. Motif background contains all other eligible peaks within the same region; window filtering remains necessary.']}, indent=2) + '\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--results', required=True, type=Path)
    p.add_argument('--annotation', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--alpha', type=float, default=.05)
    p.add_argument('--lfc', type=float, default=1)
    a = p.parse_args()
    if not 0 < a.alpha < 1 or not math.isfinite(a.lfc) or a.lfc < 0:
        p.error('Invalid thresholds')
    prepare(a.results, a.annotation, a.out, a.alpha, a.lfc)
