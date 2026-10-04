#!/usr/bin/env python3
"""Exact basewise support, counting distinct biological samples, then union across groups."""
import argparse,csv,json,collections,hashlib,math,bisect
from decisions import sha256
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

def blacklist_index(black):
 """Build once for all peaks."""
 index={}
 for c,s,e in merge(black):
  starts,ends=index.setdefault(c,([],[]));starts.append(s);ends.append(e)
 for starts,ends in index.values():
  paired=sorted(zip(starts,ends));starts[:]=[x[0] for x in paired];ends[:]=[x[1] for x in paired]
 return index

def subtract(z, black, index=None):
 """Subtract half-open intervals; optionally reuse a precomputed index."""
 if index is None:index=blacklist_index(black)
 result=[]
 for c,s,e in z:
  starts,ends=index.get(c,([],[]));i=bisect.bisect_right(ends,s);cursor=s
  while i<len(starts) and starts[i]<e:
   if cursor<starts[i]:result.append((c,cursor,min(starts[i],e)))
   cursor=max(cursor,ends[i]);i+=1
  if cursor<e:result.append((c,cursor,e))
 return result

def reproducible_union(sample_intervals,threshold):
 core=supported(sample_intervals,threshold)
 index={}
 for c,s,e in core:
  starts,ends=index.setdefault(c,([],[]));starts.append(s);ends.append(e)
 for starts,ends in index.values():
  paired=sorted(zip(starts,ends));starts[:]=[x[0] for x in paired];ends[:]=[x[1] for x in paired]
 selected=[]
 for sample in sample_intervals:
  for c,s,e in merge(sample):
   starts,ends=index.get(c,([],[]));i=bisect.bisect_right(ends,s)
   if i<len(starts) and starts[i]<e:selected.append((c,s,e))
 return merge(selected)

def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--fraction',type=float,default=2/3);p.add_argument('--blacklist',type=Path);p.add_argument('--blacklist-mode',choices=['subtract','drop'],default='subtract');p.add_argument('--universe',choices=['support_core','reproducible_union'],default='support_core');p.add_argument('--min-width',type=int,default=50);p.add_argument('--sizes',type=Path,help='Target FASTA .fai; required by planned workflows');a=p.parse_args()
 if a.min_width<1:p.error('min-width must be >=1')
 if not 0<a.fraction<=1:p.error('fraction must be >0 and <=1')
 if a.out.exists():p.error('Output directory exists')
 sizes=None
 if a.sizes:
  entries=[line.split() for line in a.sizes.read_text().splitlines() if line.strip()]
  sizes={r[0]:int(r[1]) for r in entries}
  if not sizes or len(sizes)!=len(entries) or any(n<=0 for n in sizes.values()):raise ValueError('Invalid target sequence dictionary')
 def checked_intervals(path):
  z=intervals(path)
  if sizes is not None and any(c not in sizes or e>sizes[c] for c,s,e in z):raise ValueError('Peak/blacklist interval outside target reference: '+str(path))
  return z
 rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'));groups={};bios={};sources=[]
 for r in rows:
  g,b=r['group'],r['biological_sample_id']
  if b in bios and bios[b]!=g:raise ValueError('Biological sample spans groups')
  bios[b]=g;path=Path(r['peaks_bed']);path=path if path.is_absolute() else a.manifest.parent/path
  groups.setdefault(g,{}).setdefault(b,[]).extend(checked_intervals(path))
  sources.append(dict(r,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
 if not groups:raise ValueError('Empty peak manifest')
 black=checked_intervals(a.blacklist) if a.blacklist else [];report={};master=[];decisions=[]
 index=blacklist_index(black)
 a.out.mkdir(parents=True)
 for g,samples in groups.items():
  if not g.replace('_','').isalnum():raise ValueError('Unsafe group')
  threshold=math.ceil(len(samples)*a.fraction-1e-12)
  z=(supported if a.universe=='support_core' else reproducible_union)(list(samples.values()),threshold)
  before_bp=sum(e-s for c,s,e in z);before_n=len(z);kept=[]
  for c,s,e in z:
   pieces=subtract([(c,s,e)],black,index)
   if a.blacklist_mode=='drop' and pieces!=[(c,s,e)]:pieces=[]
   accepted=[v for v in pieces if v[2]-v[1]>=a.min_width]
   decisions.append({'group':g,'source_interval':[c,s,e],'retained':accepted,'short_pieces_discarded':[v for v in pieces if v not in accepted]})
   kept.extend(accepted)
  z=merge(kept)
  (a.out/(g+'.consensus.bed')).write_text(''.join(f'{c}\t{s}\t{e}\n' for c,s,e in z))
  report[g]={'independent_sample_count':len(samples),'support_required':threshold,'interval_count':len(z),'before_blacklist_intervals':before_n,'before_blacklist_bp':before_bp,'retained_bp':sum(e-s for c,s,e in z)};master.extend(z)
 z=merge(master)
 (a.out/'master.bed').write_text(''.join(f'{c}\t{s}\t{e}\tpeak_{i:08d}\n' for i,(c,s,e) in enumerate(z,1)))
 (a.out/'provenance.json').write_text(json.dumps({'groups':report,'sources':sources,'fraction':a.fraction,'universe':a.universe,'blacklist_mode':a.blacklist_mode,'min_width':a.min_width,'blacklist_sha256':sha256(a.blacklist) if a.blacklist else None,'rule':'support core or full sample intervals intersecting support core; blacklist bases subtracted by default (drop mode retained only for sensitivity); union across groups','master_intervals':len(z)},indent=2))
 (a.out/'interval_decisions.json').write_text(json.dumps(decisions,indent=2))
if __name__=='__main__':main()
