#!/usr/bin/env python3
"""Pre-register MACS2 paired-end control sensitivity commands; never selects a winner."""
import argparse,csv,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--gsize',required=True);p.add_argument('--execute',action='store_true');a=p.parse_args()
if a.out.exists():p.error('Output exists')
rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'));a.out.mkdir(parents=True);reports=[]
for r in rows:
 sid=r['sample_id']
 if not sid.replace('_','').isalnum():raise ValueError('Unsafe ID')
 for mode in ('igg_default','igg_scale_large','no_igg_diagnostic'):
  dest=a.out/(sid+'_'+mode);dest.mkdir()
  cmd=['macs2','callpeak','-t',str(Path(r['target_bam']).resolve()),'-f','BAMPE','-g',a.gsize,'-n',sid,'--outdir',str(dest.resolve()),'--keep-dup','all','-q','0.05']
  if mode!='no_igg_diagnostic':cmd+=['-c',str(Path(r['control_bam']).resolve())]
  if mode=='igg_scale_large':cmd+=['--scale-to','large']
  report={'sample_id':sid,'strategy':mode,'command':cmd,'state':'PLANNED','interpretation':'diagnostic sensitivity comparison; production acceptance required'}
  if a.execute:
   with (dest/'execution.log').open('w') as f:result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
   output=dest/(sid+'_peaks.narrowPeak');report.update(state='COMPUTATIONAL_PASS' if result.returncode==0 and output.is_file() else 'FAILED',returncode=result.returncode,peak_count=len(output.read_text().splitlines()) if output.is_file() else None)
  reports.append(report);(a.out/'diagnostics.json').write_text(json.dumps(reports,indent=2))
if any(r['state']=='FAILED' for r in reports):raise RuntimeError('One or more MACS2 diagnostics failed')
