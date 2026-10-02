#!/usr/bin/env python3
"""Known HOMER motifs on equal-length nonoverlapping windows and explicit tested background."""
import argparse,json,subprocess
from pathlib import Path
from consensus import intervals
p=argparse.ArgumentParser();p.add_argument('--target',required=True,type=Path);p.add_argument('--background',required=True,type=Path);p.add_argument('--fasta',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--threads',type=int,default=4);p.add_argument('--width',type=int,default=200);p.add_argument('--min-windows',type=int,default=50);p.add_argument('--execute',action='store_true');a=p.parse_args()
if a.out.exists():p.error('Output exists')
if a.width<=0 or a.width%2 or a.threads<1:p.error('Positive even width and positive threads required')
a.out.mkdir(parents=True);excluded=[]
def windows(path,label):
 keep=[]
 for c,s,e in sorted(intervals(path)):
  mid=(s+e)//2;z=(c,mid-a.width//2,mid+a.width//2)
  if z[1]<0 or (keep and keep[-1][0]==c and z[1]<keep[-1][2]):excluded.append((label,*z));continue
  keep.append(z)
 return keep
z=windows(a.target,'target');b=windows(a.background,'background');b=[v for v in b if not any(v[0]==t[0] and v[1]<t[2] and t[1]<v[2] for t in z)]
for label,v in [('target',z),('background',b)]:
 (a.out/(label+'.bed')).write_text(''.join(f'{c}\t{s}\t{e}\t{label}_{i}\n' for i,(c,s,e) in enumerate(v)))
(a.out/'excluded.json').write_text(json.dumps(excluded))
cmd=['findMotifsGenome.pl',str((a.out/'target.bed').resolve()),str(a.fasta.resolve()),str((a.out/'homer').resolve()),'-bg',str((a.out/'background.bed').resolve()),'-size','given','-gc','-nomotif','-mset','vertebrates','-p',str(a.threads)]
state='PLANNED' if min(len(z),len(b))>=a.min_windows else 'SKIPPED_TOO_FEW_WINDOWS'
report={'state':state,'target_windows':len(z),'background_windows':len(b),'command':cmd}
(a.out/'status.json').write_text(json.dumps(report,indent=2))
if a.execute and state=='PLANNED':
 if not a.fasta.is_file():raise ValueError('Missing reference FASTA')
 with (a.out/'execution.log').open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 output=a.out/'homer/knownResults.txt';report['state']='COMPUTATIONAL_PASS' if r.returncode==0 and output.is_file() and output.stat().st_size>0 else 'FAILED';report['returncode']=r.returncode
 (a.out/'status.json').write_text(json.dumps(report,indent=2))
 if report['state']=='FAILED':raise RuntimeError('HOMER failed or missing knownResults.txt')
