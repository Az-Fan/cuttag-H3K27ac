#!/usr/bin/env python3
"""Build paired-fragment BED from BAM while preserving duplicates by default."""
import argparse,collections,csv,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--bam',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--mapq',type=int,default=20);p.add_argument('--threads',type=int,default=4);p.add_argument('--remove-duplicates',action='store_true');a=p.parse_args()
if a.out.exists():p.error('Output exists')
if not 0<=a.mapq<=255 or a.threads<1:p.error('Invalid MAPQ/thread count')
a.out.mkdir(parents=True);filtered=a.out/'namesorted.bam';
# Exclude unmapped, mate unmapped, secondary, QC fail, supplementary; retain duplicate-marked reads unless requested.
exclude=4|8|256|512|2048|(1024 if a.remove_duplicates else 0)
cmd1=['samtools','view','-b','-f','2','-F',str(exclude),'-q',str(a.mapq),str(a.bam)]
cmd2=['samtools','sort','-n','-@',str(a.threads),'-o',str(filtered),'-']
with (a.out/'filter_sort.log').open('w') as log:
 first=subprocess.Popen(cmd1,stdout=subprocess.PIPE,stderr=log);second=subprocess.run(cmd2,stdin=first.stdout,stderr=log);first.stdout.close();rc=first.wait()
if rc or second.returncode:raise RuntimeError('samtools filter/name sort failed')
cmd=['samtools','view',str(filtered)];(a.out/'samtools_view.command.json').write_text(json.dumps(cmd))
lengths=collections.Counter();n=0
def emit(pair,out):
 global n
 if len(pair)!=2:return
 one=next((f for f in pair if int(f[1])&0x40),None);two=next((f for f in pair if int(f[1])&0x80),None)
 if not one or not two or one[2]!=two[2] or one[2]=='*':return
 if int(one[4])<a.mapq or int(two[4])<a.mapq:return
 tlen=int(one[8]);start=min(int(one[3])-1,int(two[3])-1);end=start+abs(tlen)
 if tlen==0 or start<0 or end<=start:return
 out.write(f'{one[2]}\t{start}\t{end}\t{one[0]}\n');n+=1;lengths[end-start]+=1
proc=subprocess.Popen(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
with (a.out/'fragments.bed').open('w') as out:
 current=None;pair=[]
 for line in proc.stdout:
  f=line.rstrip('\n').split('\t')
  if current is not None and f[0]!=current:emit(pair,out);pair=[]
  current=f[0];pair.append(f)
 if pair:emit(pair,out)
 stderr=proc.stderr.read();rc=proc.wait()
 if rc:raise RuntimeError('samtools view failed: '+stderr[-2000:])
if not n:raise ValueError('No usable paired fragments')
with (a.out/'length_histogram.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['length','fragments']);w.writerows(sorted(lengths.items()))
(a.out/'summary.json').write_text(json.dumps({'fragments':n,'minimum_MAPQ_both_mates':a.mapq,'duplicates_removed':a.remove_duplicates,'unit':'paired fragment; TLEN from primary proper pair'},indent=2))
