#!/usr/bin/env python3
"""Portable QC summary with explicit unavailable metrics and descriptive plots."""
import argparse,csv,json,html,math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser();p.add_argument('--metrics',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
if a.out.exists():p.error('Output exists')
rows=list(csv.DictReader(a.metrics.open(),delimiter='\t'))
if not rows:raise ValueError('Empty QC table')
a.out.mkdir(parents=True);metrics=['target_fragments','spikein_fraction','duplication_rate','frip','peak_count','median_fragment_length']
fig,axes=plt.subplots(2,3,figsize=(14,8));missing={}
for ax,key in zip(axes.flat,metrics):
 vals=[];labels=[];missing[key]=[]
 for r in rows:
  try:value=float(r.get(key,''))
  except (ValueError,TypeError):missing[key].append(r['sample_id']);continue
  if not math.isfinite(value):missing[key].append(r['sample_id']);continue
  if value<0 or (key in ('spikein_fraction','duplication_rate','frip') and value>1):raise ValueError('QC metric outside valid range: '+key)
  vals.append(value);labels.append(r['sample_id'])
 if vals:ax.bar(labels,vals);ax.tick_params(axis='x',rotation=60)
 else:ax.text(.5,.5,'Unavailable',ha='center')
 ax.set_title(key)
fig.tight_layout();fig.savefig(a.out/'qc_atlas.pdf');fig.savefig(a.out/'qc_atlas.png');plt.close(fig)
(a.out/'missing_metrics.json').write_text(json.dumps(missing,indent=2))
cols=list(rows[0]);table='<table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in cols)+'</tr>'
for r in rows:table+='<tr>'+''.join('<td>'+html.escape(r[k] or '')+'</td>' for k in cols)+'</tr>'
(a.out/'index.html').write_text('<!doctype html><meta charset="utf-8"><h1>QC review</h1><p>Descriptive metrics; sample acceptance requires a documented decision.</p><img src="qc_atlas.png" width="100%">'+table+'</table>')
