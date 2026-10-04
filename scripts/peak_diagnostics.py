#!/usr/bin/env python3
"""Pre-register MACS2 paired-end control sensitivity commands; never selects a winner."""
import argparse,csv,json,subprocess,itertools
from decisions import sha256
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);p.add_argument('--gsize',required=True);p.add_argument('--execute',action='store_true');p.add_argument('--macs2-container',type=Path);p.add_argument('--native-macs2',action='store_true');p.add_argument('--qvalue',type=float,default=.05);p.add_argument('--broad-cutoff',type=float,default=.1);p.add_argument('--shapes',nargs='+',choices=['narrow','broad'],default=['narrow','broad']);p.add_argument('--backgrounds',nargs='+',choices=['igg_default','igg_scale_large','no_igg_diagnostic'],default=['igg_default','igg_scale_large','no_igg_diagnostic']);p.add_argument('--duplicate-modes',nargs='+',choices=['all','auto'],default=['all']);a=p.parse_args()
if not 0<a.qvalue<=a.broad_cutoff<1:p.error('Require 0 < qvalue <= broad-cutoff < 1')
if a.native_macs2 and a.macs2_container:p.error('Choose native or container')
if not a.native_macs2 and not a.macs2_container:
 from fetch_tools import check_tool
 a.macs2_container=check_tool('macs2')
if a.out.exists():p.error('Output exists')
rows=list(csv.DictReader(a.manifest.open(),delimiter='\t'))
if not rows or len({r['sample_id'] for r in rows})!=len(rows):raise ValueError('Unique nonempty sample manifest required')
for r in rows:
 if not r.get('sample_id') or not (a.manifest.resolve().parent/r.get('target_bam','')).is_file():raise ValueError('Missing sample_id or target_bam')
 if r.get('bam_policy')!='duplicates_retained':raise ValueError('Peak strategy artifact export requires duplicate-retained target BAMs')
 if any(mode!='no_igg_diagnostic' for mode in a.backgrounds) and not r.get('control_bam','').strip():raise ValueError('Selected IgG strategy requires control_bam')
 if r.get('control_bam','').strip() and not (a.manifest.resolve().parent/r['control_bam']).is_file():raise ValueError('Missing control_bam')
a.out.mkdir(parents=True);reports=[]
artifact_manifests={}
for r in rows:
 sid=r['sample_id']
 if not sid.replace('_','').isalnum():raise ValueError('Unsafe ID')
 for shape,mode,duplicate_mode in itertools.product(a.shapes,a.backgrounds,a.duplicate_modes):
  strategy=shape+'_'+mode+'_dup'+duplicate_mode
  dest=a.out/(sid+'_'+strategy);dest.mkdir()
  target=(a.manifest.resolve().parent/Path(r['target_bam'])).resolve()
  control=(a.manifest.resolve().parent/Path(r['control_bam'])).resolve() if r.get('control_bam','').strip() else None
  cmd=['macs2','callpeak','-t',str(target),'-f','BAMPE','-g',a.gsize,'-n',sid,'--outdir',str(dest.resolve()),'--keep-dup',duplicate_mode,'-q',str(a.qvalue)]
  if shape=='broad':cmd+=['--broad','--broad-cutoff',str(a.broad_cutoff)]
  if mode!='no_igg_diagnostic':cmd+=['-c',str(control)]
  if mode=='igg_scale_large':cmd+=['--scale-to','large']
  if a.macs2_container:
   bind_paths={str(a.out.resolve()),str(target.parent)}
   if control:bind_paths.add(str(control.parent))
   binds=sorted(bind_paths)
   cmd=['apptainer','exec','--cleanenv','--no-home']+[v for bind in binds for v in ('--bind',bind)]+[str(a.macs2_container.resolve())]+cmd
  report={'sample_id':sid,'strategy':strategy,'shape':shape,'background':mode,'duplicate_mode':duplicate_mode,'qvalue':a.qvalue,'broad_cutoff':a.broad_cutoff if shape=='broad' else None,'target_sha256':sha256(target),'control_sha256':sha256(control) if control and mode!='no_igg_diagnostic' else None,'command':cmd,'state':'PLANNED','interpretation':'diagnostic sensitivity comparison; production acceptance required','container_sha256':sha256(a.macs2_container) if a.macs2_container else None}
  if a.execute:
   with (dest/'execution.log').open('w') as f:result=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
   output=dest/(sid+'_peaks.'+shape+'Peak');report.update(state='COMPUTATIONAL_PASS' if result.returncode==0 and output.is_file() else 'FAILED',returncode=result.returncode,peak_count=len(output.read_text().splitlines()) if output.is_file() else None)
  reports.append(report)
if a.execute:
 for strategy in sorted({r['strategy'] for r in reports}):
  selected=[r for r in reports if r['strategy']==strategy]
  if len(selected)!=len(rows) or any(r['state']!='COMPUTATIONAL_PASS' for r in selected):
   continue
  artifact_dir=a.out/'artifact_sets';artifact_dir.mkdir(exist_ok=True)
  artifact_path=artifact_dir/(strategy+'.tsv')
  with artifact_path.open('w') as f:
   writer=csv.writer(f,delimiter='\t');writer.writerow(['sample_id','bam','peaks_bed','bam_policy'])
   for report in selected:
    source=next(r for r in rows if r['sample_id']==report['sample_id'])
    bam=(a.manifest.resolve().parent/source['target_bam']).resolve()
    peak=(a.out/(report['sample_id']+'_'+strategy)/(report['sample_id']+'_peaks.'+report['shape']+'Peak')).resolve()
    writer.writerow([report['sample_id'],str(bam),str(peak),'duplicates_retained'])
  artifact_manifests[strategy]=str(artifact_path.resolve())
(a.out/'diagnostics.json').write_text(json.dumps({'state':'REVIEW_REQUIRED','candidates':reports,
 'artifact_manifests':artifact_manifests,'interpretation':'Artifact manifests are generated only when every sample in a strategy completed; scientific review remains required.'},indent=2)+'\n')
if any(r['state']=='FAILED' for r in reports):raise RuntimeError('One or more MACS2 diagnostics failed')
