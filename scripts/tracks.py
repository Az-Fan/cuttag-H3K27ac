#!/usr/bin/env python3
"""Fragment-depth tracks; spike-in mode needs externally reviewed factor per target sample."""
import argparse,csv,json,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
import pyBigWig
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--sizes',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
if a.out.exists():p.error('Output exists')
header=[(r[0],int(r[1])) for r in (l.split() for l in a.sizes.read_text().splitlines())];sizes=dict(header)
rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'))
if not rows:raise ValueError('Empty track manifest')
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
  if sf<=0:raise ValueError('spikein_scale must be positive relative track multiplier')
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
  session=ET.Element('Session',genome='unknown',version='8');resources=ET.SubElement(session,'Resources');ET.SubElement(resources,'Resource',path=bwpath.name);ET.ElementTree(session).write(a.out/(sample+'.'+mode+'.igv.xml'),encoding='UTF-8',xml_declaration=True)
  (a.out/(sample+'.'+mode+'.json')).write_text(json.dumps({'sample_id':sample,'normalization':mode,'retained_target_fragments':count,'multiplier':factor,'unit':'fragment coverage per base','input_sha256':__import__('hashlib').sha256(fragments.read_bytes()).hexdigest(),'command':cmd},indent=2))
