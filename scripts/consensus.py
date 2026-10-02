#!/usr/bin/env python3
"""Exact basewise support, counting distinct biological samples, then union across groups."""
import argparse,csv,json,collections,hashlib,math
from pathlib import Path

def intervals(path):
 out=[]
 for line in Path(path).read_text().splitlines():
  if not line.strip() or line.startswith(('#','track','browser')):continue
  r=line.split();c,s,e=r[0],int(r[1]),int(r[2])
  if s<0 or e<=s:raise ValueError('Invalid interval: '+line)
  out.append((c,s,e))
 return out

def merge(z):
 out=[]
 for c,s,e in sorted(z):
  if out and out[-1][0]==c and s<=out[-1][2]:out[-1]=(c,out[-1][1],max(e,out[-1][2]))
  else:out.append((c,s,e))
 return out

def supported(sample_intervals,threshold):
 events=collections.defaultdict(lambda:collections.defaultdict(int))
 for z in sample_intervals:
  for c,s,e in merge(z):events[c][s]+=1;events[c][e]-=1
 result=[]
 for c,changes in sorted(events.items()):
  active=0;previous=None
  for pos,delta in sorted(changes.items()):
   if previous is not None and active>=threshold and previous<pos:result.append((c,previous,pos))
   active+=delta;previous=pos
 return merge(result)

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--fraction',type=float,default=2/3);p.add_argument('--blacklist',type=Path);a=p.parse_args()
 if not 0<a.fraction<=1:p.error('fraction must be >0 and <=1')
 if a.out.exists():p.error('Output directory exists')
 rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'));groups={};bios={};sources=[]
 for r in rows:
  g,b=r['group'],r['biological_sample_id']
  if b in bios and bios[b]!=g:raise ValueError('Biological sample spans groups')
  bios[b]=g;path=Path(r['peaks_bed']);path=path if path.is_absolute() else a.manifest.parent/path
  groups.setdefault(g,{}).setdefault(b,[]).extend(intervals(path))
  sources.append(dict(r,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 if not groups:raise ValueError('Empty peak manifest')
 black=intervals(a.blacklist) if a.blacklist else [];report={};master=[]
 a.out.mkdir(parents=True)
 for g,samples in groups.items():
  if not g.replace('_','').isalnum():raise ValueError('Unsafe group')
  threshold=math.ceil(len(samples)*a.fraction-1e-12);z=supported(list(samples.values()),threshold)
  z=[(c,s,e) for c,s,e in z if not any(c==bc and s<be and bs<e for bc,bs,be in black)]
  (a.out/(g+'.consensus.bed')).write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in z))
  report[g]={'independent_sample_count':len(samples),'support_required':threshold,'interval_count':len(z)};master.extend(z)
 z=merge(master)
 (a.out/'master.bed').write_text(''.join(f'{c}\t{s}\t{e}\tpeak_{i:08d}\n' for i,(c,s,e) in enumerate(z,1)))
 (a.out/'provenance.json').write_text(json.dumps({'groups':report,'sources':sources,'fraction':a.fraction,'rule':'basewise biological sample support; blacklist-overlap intervals discarded; union across groups','master_intervals':len(z)},indent=2))
if __name__=='__main__':main()
