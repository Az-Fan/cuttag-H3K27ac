#!/usr/bin/env python3
"""Audit concordant pair counts from Bowtie2 logs; never auto-accept calibration."""
import argparse,csv,json,math,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--equal-amount-confirmed',action='store_true');p.add_argument('--added-at',default='unknown');a=p.parse_args()
if a.out.exists():p.error('Output exists')
rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'));result=[]
for r in rows:
 path=Path(r['bowtie2_log']);path=path if path.is_absolute() else a.manifest.parent/path
 t=path.read_text();total=re.search(r'(\d+) reads; of these:',t);zero=re.search(r'(\d+) \([^\n]+\) aligned concordantly 0 times',t);one=re.search(r'(\d+) \([^\n]+\) aligned concordantly exactly 1 time',t);many=re.search(r'(\d+) \([^\n]+\) aligned concordantly >1 times',t)
 if not all((total,zero,one,many)):raise ValueError('Expected paired-end Bowtie2 summary '+str(path))
 n=int(total[1]);count=int(one[1])+int(many[1])
 if n<=0 or int(zero[1])+count!=n:raise ValueError('Pair totals inconsistent')
 result.append({'sample_id':r['sample_id'],'mode':r['mode'],'calibration_group':r.get('calibration_group','experiment'),'total_pairs':n,'spikein_fragments':count,'fraction':count/n,'calibration_accepted':False})
# Multiple modes intentionally kept separate; normalize only within mode across unique samples.
for mode,group in {(r['mode'],r['calibration_group']) for r in result}:
 z=[r for r in result if r['mode']==mode and r['calibration_group']==group]
 if len({r['sample_id'] for r in z})!=len(z):raise ValueError('Duplicate sample/mode')
 positive=all(r['spikein_fragments']>0 for r in z)
 gm=math.exp(sum(math.log(r['spikein_fragments']) for r in z)/len(z)) if positive else None
 for r in z:r['size_factor']=r['spikein_fragments']/gm if gm else None;r['track_scale_relative']=gm/r['spikein_fragments'] if gm else None
if not result:raise ValueError('Empty manifest')
a.out.mkdir(parents=True)
with (a.out/'counts.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
(a.out/'review.json').write_text(json.dumps({'equal_amount_confirmed':a.equal_amount_confirmed,'added_at':a.added_at,'accepted':False,'required_review':['same reads across diagnostic modes','geometry, MAPQ and multi-mapping','full-data count confirmation','experimental calibration scope'],'rows':result},indent=2))
