#!/usr/bin/env python3
"""Prepare a content-bound review draft. All scientific acceptance fields remain false."""
import argparse,json
from pathlib import Path
from cuttag import load,resolve
from decisions import sha256
p=argparse.ArgumentParser();p.add_argument('--config',required=True,type=Path);p.add_argument('--kind',required=True,choices=['production','model']);p.add_argument('--counts',type=Path);p.add_argument('--metadata',type=Path);p.add_argument('--spikein',type=Path);p.add_argument('--evidence',required=True,type=Path,nargs='+');p.add_argument('--out',required=True,type=Path);a=p.parse_args()
if a.out.exists():p.error('Review already exists; use a new path')
root,cfg,_=load(a.config);inputs={'config_sha256':sha256(a.config),'samples_sha256':sha256(resolve(root,cfg['samples']))}
keys=['upstream_accepted','control_strategy_accepted'] if a.kind=='production' else ['upstream_accepted','peaks_accepted','replicates_confirmed','simple_design_accepted']
if a.kind=='production':
 inputs['fastq_inventory_sha256']=sha256(root/'results/qc/fastq_inventory.json')
 inputs['reference_audit_sha256']=sha256(root/'results/qc/reference_audit.json')
if a.kind=='model':
 if not a.counts or not a.metadata:p.error('Model review needs counts and metadata')
 inputs.update(counts_sha256=sha256(a.counts),metadata_sha256=sha256(a.metadata))
if cfg['spikein']['enabled'] or cfg['analysis']['normalization']=='spikein':
 keys+=['spikein_counting_accepted','spikein_calibration_accepted']
 if cfg['spikein'].get('alignment_profile'):inputs['alignment_profile_sha256']=sha256(resolve(root,cfg['spikein']['alignment_profile']))
 if a.spikein:inputs['spikein_sha256']=sha256(a.spikein)
review={'reviewer':'','reviewed_at':'','reason':'','decisions':{k:False for k in keys},'inputs':inputs,'model':{k:cfg['analysis'][k] for k in ('numerator','denominator','design','normalization','alpha','abs_log2fc','min_total_count')},'spikein':cfg['spikein'],'evidence':[{'path':str(f.resolve()),'sha256':sha256(f)} for f in a.evidence]}
a.out.parent.mkdir(parents=True,exist_ok=True)
with a.out.open('x') as f:json.dump(review,f,indent=2,ensure_ascii=False);f.write('\n')
