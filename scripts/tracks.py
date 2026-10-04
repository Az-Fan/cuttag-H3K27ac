#!/usr/bin/env python3
"""Fragment-depth tracks; spike-in mode needs externally reviewed factor per target sample."""
import argparse,csv,json,subprocess,math,xml.etree.ElementTree as ET
from decisions import require_review,sha256,known_text
from pathlib import Path
import pyBigWig
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--sizes',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--genome',required=True);p.add_argument('--review',type=Path);a=p.parse_args()
if a.out.exists():p.error('Output exists')
header=[(r[0],int(r[1])) for r in (l.split() for l in a.sizes.read_text().splitlines())];sizes=dict(header)
rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'))
if not rows:raise ValueError('Empty track manifest')
if len(sizes)!=len(header) or any(n<=0 for _,n in header):raise ValueError('Invalid/duplicate chromosome sizes')
ids=[r['sample_id'] for r in rows]
if len(ids)!=len(set(ids)):raise ValueError('Duplicate track sample IDs')
for r in rows:r['spikein_scale']=(r.get('spikein_scale') or '').strip()
has_scale=[bool((r.get('spikein_scale') or '').strip()) for r in rows]
if any(has_scale) and not all(has_scale):raise ValueError('All samples must have spikein_scale, or all must omit it')
if all(has_scale):
 if not a.review:raise ValueError('Spike-in tracks require calibration review')
 review=require_review(a.review,['spikein_counting_accepted','spikein_calibration_accepted'],{'track_manifest_sha256':a.manifest,'sizes_sha256':a.sizes})
 sp=review.get('spikein',{})
 if sp.get('equal_amount_confirmed') is not True or not known_text(sp.get('added_at')) or not known_text(sp.get('calibration_scope')):raise ValueError('Spike-in experiment scope incomplete')
 for r in rows:
  if r.get('spikein_scale') and (not math.isfinite(float(r['spikein_scale'])) or float(r['spikein_scale'])<=0):raise ValueError('Nonfinite/invalid spike-in scale')
a.out.mkdir(parents=True)
for r in rows:
 sample=r['sample_id']
 if not sample.replace('_','').isalnum():raise ValueError('Unsafe sample ID')
 fragments=Path(r['fragments_bed']);fragments=fragments if fragments.is_absolute() else a.manifest.parent/fragments
 count=0
 with fragments.open() as f:
  for line in f:
   z=line.split();c,s,e=z[0],int(z[1]),int(z[2])
   if c not in sizes or s<0 or e<=s or e>sizes[c]:raise ValueError('Fragment outside target reference')
   count+=1
 if count==0:raise ValueError('Empty fragments '+sample)
 sortedbed=a.out/(sample+'.sorted.bed')
 with sortedbed.open('w') as f:subprocess.run(['bedtools','sort','-i',str(fragments.resolve()),'-g',str(a.sizes.resolve())],stdout=f,check=True)
 factors={'CPM':1e6/count}
 if r.get('spikein_scale'):
  sf=float(r['spikein_scale'])
  if not math.isfinite(sf) or sf<=0:raise ValueError('spikein_scale must be positive relative track multiplier')
  factors['spikein']=sf
 for mode,factor in factors.items():
  cmd=['bedtools','genomecov','-bg','-i',str(sortedbed),'-g',str(a.sizes.resolve()),'-scale',str(factor)]
  bg=a.out/(sample+'.'+mode+'.bedGraph')
  with bg.open('w') as f:subprocess.run(cmd,stdout=f,check=True)
  bwpath=a.out/(sample+'.'+mode+'.bw');bw=pyBigWig.open(str(bwpath),'w');bw.addHeader(header)
  try:
   with bg.open() as f:
    for line in f:
     c,s,e,v=line.split();bw.addEntries([c],[int(s)],ends=[int(e)],values=[float(v)])
  finally:bw.close()
  session=ET.Element('Session',genome=a.genome,version='8');resources=ET.SubElement(session,'Resources');ET.SubElement(resources,'Resource',path=bwpath.name);ET.ElementTree(session).write(a.out/(sample+'.'+mode+'.igv.xml'),encoding='UTF-8',xml_declaration=True)
  (a.out/(sample+'.'+mode+'.json')).write_text(json.dumps({'sample_id':sample,'normalization':mode,'retained_target_fragments':count,'multiplier':factor,'unit':'fragment coverage per base','input_sha256':sha256(fragments),'command':cmd},indent=2))
