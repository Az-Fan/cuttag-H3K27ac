#!/usr/bin/env python3
"""Pre-register MACS2 paired-end control sensitivity commands; never selects a winner."""
import argparse,csv,json,subprocess
from decisions import sha256
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--gsize',required=True);p.add_argument('--execute',action='store_true');p.add_argument('--macs2-container',type=Path);p.add_argument('--native-macs2',action='store_true');p.add_argument('--qvalue',type=float,default=.05);p.add_argument('--broad-cutoff',type=float,default=.1);a=p.parse_args()
if not 0<a.qvalue<=a.broad_cutoff<1:p.error('Require 0 < qvalue <= broad-cutoff < 1')
if a.native_macs2 and a.macs2_container:p.error('Choose native or container')
if not a.native_macs2 and not a.macs2_container:
 from fetch_tools import check_tool
 a.macs2_container=check_tool('macs2')
if a.out.exists():p.error('Output exists')
rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'))
if not rows or len({r['sample_id'] for r in rows})!=len(rows):raise ValueError('Unique nonempty sample manifest required')
for r in rows:
 for key in ('target_bam','control_bam'):
  if not (a.manifest.resolve().parent/r[key]).is_file():raise ValueError('Missing '+key)
a.out.mkdir(parents=True);reports=[]
for r in rows:
 sid=r['sample_id']
 if not sid.replace('_','').isalnum():raise ValueError('Unsafe ID')
 for shape,mode in __import__('itertools').product(('narrow','broad'),('igg_default','igg_scale_large','no_igg_diagnostic')):
  dest=a.out/(sid+'_'+shape+'_'+mode);dest.mkdir()
  cmd=['macs2','callpeak','-t',str((a.manifest.resolve().parent/Path(r['target_bam'])).resolve()),'-f','BAMPE','-g',a.gsize,'-n',sid,'--outdir',str(dest.resolve()),'--keep-dup','all','-q',str(a.qvalue)]
  if shape=='broad':cmd+=['--broad','--broad-cutoff',str(a.broad_cutoff)]
  if mode!='no_igg_diagnostic':cmd+=['-c',str((a.manifest.resolve().parent/Path(r['control_bam'])).resolve())]
  if mode=='igg_scale_large':cmd+=['--scale-to','large']
  if a.macs2_container:
   binds=sorted({str(a.out.resolve()),str((a.manifest.resolve().parent/r['target_bam']).resolve().parent),str((a.manifest.resolve().parent/r['control_bam']).resolve().parent)})
   cmd=['apptainer','exec','--cleanenv','--no-home']+[v for bind in binds for v in ('--bind',bind)]+[str(a.macs2_container.resolve())]+cmd
  report={'sample_id':sid,'strategy':shape+'_'+mode,'shape':shape,'background':mode,'qvalue':a.qvalue,'broad_cutoff':a.broad_cutoff if shape=='broad' else None,'target_sha256':sha256(a.manifest.resolve().parent/r['target_bam']),'control_sha256':sha256(a.manifest.resolve().parent/r['control_bam']) if mode!='no_igg_diagnostic' else None,'command':cmd,'state':'PLANNED','interpretation':'diagnostic sensitivity comparison; production acceptance required','container_sha256':sha256(a.macs2_container) if a.macs2_container else None}
  if a.execute:
   with (dest/'execution.log').open('w') as f:result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
   output=dest/(sid+'_peaks.'+shape+'Peak');report.update(state='COMPUTATIONAL_PASS' if result.returncode==0 and output.is_file() else 'FAILED',returncode=result.returncode,peak_count=len(output.read_text().splitlines()) if output.is_file() else None)
  reports.append(report);(a.out/'diagnostics.json').write_text(json.dumps(reports,indent=2))
if any(r['state']=='FAILED' for r in reports):raise RuntimeError('One or more MACS2 diagnostics failed')
