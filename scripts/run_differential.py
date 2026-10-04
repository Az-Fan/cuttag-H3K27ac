#!/usr/bin/env python3
"""Configured differential entry: exact sample mapping, accepted peak QC and one source of thresholds."""
import argparse,csv,os,subprocess
from pathlib import Path
from cuttag import load,validate,resolve
from decisions import require_review
p=argparse.ArgumentParser();p.add_argument('--config',required=True,type=Path);p.add_argument('--counts',required=True,type=Path);p.add_argument('--metadata',required=True,type=Path);p.add_argument('--review',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--spikein',type=Path);a=p.parse_args()
root,cfg,rows=load(a.config);report=validate(root,cfg,rows,missing=True)
if report['errors'] or not report['formal_inference_eligible']:raise ValueError('Design not eligible: '+str(report))
if cfg['qc']['upstream_accepted'] is not True or cfg['qc']['peaks_accepted'] is not True:raise ValueError('Upstream and peak QC must be accepted')
an=cfg['analysis'];target={r['sample_id']:(r['biological_sample_id'],r['condition']) for r in rows if r['role']=='target' and r['condition'] in (an['numerator'],an['denominator'])}
with a.metadata.open() as f:meta=list(csv.DictReader(f,delimiter='\t'))
actual={r['sample_id']:(r['biological_sample_id'],r['condition']) for r in meta}
if actual!=target or len(actual)!=len(meta):raise ValueError('Metadata does not match configured target biological samples')
require_review(a.review,['upstream_accepted','peaks_accepted','replicates_confirmed','simple_design_accepted'],{'config_sha256':a.config,'counts_sha256':a.counts,'metadata_sha256':a.metadata})
cmd=['Rscript',str(Path(__file__).with_name('differential.R')),str(a.counts),str(a.metadata),an['numerator'],an['denominator'],str(a.out),an['normalization']]
if an['normalization']=='spikein':
 if not a.spikein:raise ValueError('Spike-in table required')
 sp=cfg['spikein']
 if sp['enabled'] is not True or sp['calibration_accepted'] is not True or sp['equal_amount_confirmed'] is not True or sp['added_at'] in ('','unknown',None):raise ValueError('Configured spike-in not accepted')
 with a.spikein.open() as f:spike_rows=list(csv.DictReader(f,delimiter='\t'))
 if not spike_rows or not {'sample_id','spikein_fragments','calibration_accepted','calibration_group'}.issubset(spike_rows[0]):raise ValueError('Spike-in table must be target-only and include calibration_group')
 if {r['sample_id'] for r in spike_rows}!=set(target) or len(spike_rows)!=len(target):raise ValueError('Spike-in rows must match configured target biological samples exactly')
 groups={r['calibration_group'] for r in spike_rows}
 if '' in groups or len(groups)!=1:raise ValueError('The current ~ condition model requires one shared calibration_group; multiple relative scales need a reviewed group-aware design')
 if any(r['calibration_accepted'].strip().lower() not in ('true','1') for r in spike_rows):raise ValueError('Every target spike-in row must be accepted in the reviewed table')
 cmd.append(str(a.spikein))
cmd.append(str(a.review));env=os.environ.copy();env.update(CUTTAG_ALPHA=str(an['alpha']),CUTTAG_LFC=str(an['abs_log2fc']),CUTTAG_MIN_COUNT=str(an['min_total_count']))
raise SystemExit(subprocess.run(cmd,env=env).returncode)
