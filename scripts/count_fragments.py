#!/usr/bin/env python3
"""Count BED paired fragments against an accepted master BED; each input row is ONE fragment."""
import argparse,csv,json,subprocess,tempfile
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--peaks',required=True,type=Path);p.add_argument('--samples',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
rows=list(csv.DictReader(a.samples.open(),delimiter='\t'))
ids=[r['sample_id'] for r in rows]
if not rows or any(not x for x in ids) or len(ids)!=len(set(ids)):p.error('Unique nonempty biological sample IDs required')
peaks=[l.split() for l in a.peaks.read_text().splitlines() if l.strip()]
if not peaks or any(len(r)!=4 or int(r[1])<0 or int(r[2])<=int(r[1]) for r in peaks) or len({r[3] for r in peaks})!=len(peaks):p.error('BED4 with unique peak IDs required')
if a.out.exists():p.error('Output exists')
columns=[];commands=[]
for r in rows:
 fragment=Path(r['fragments_bed']);fragment=fragment if fragment.is_absolute() else a.samples.resolve().parent/fragment
 cmd=['bedtools','coverage','-counts','-a',str(a.peaks),'-b',str(fragment)];commands.append(cmd)
 result=subprocess.run(cmd,text=True,capture_output=True,check=True)
 lines=[l.split('\t') for l in result.stdout.splitlines()]
 if [l[:4] for l in lines]!=peaks:raise ValueError('Coverage coordinates/order changed')
 columns.append([int(l[-1]) for l in lines])
a.out.parent.mkdir(parents=True,exist_ok=True)
with a.out.open('x') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['peak_id','chrom','start','end']+ids)
 for i,r in enumerate(peaks):w.writerow([r[3],r[0],r[1],r[2]]+[c[i] for c in columns])
a.out.with_suffix('.commands.json').write_text(json.dumps(commands,indent=2))
