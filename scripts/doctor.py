#!/usr/bin/env python3
"""Probe executable runtime paths, not only installation metadata/version strings."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
from decisions import sha256

p=argparse.ArgumentParser()
p.add_argument('--out',required=True,type=Path)
p.add_argument('--macs2-container',type=Path)
a=p.parse_args()
if not a.macs2_container:
 from fetch_tools import check_tool
 a.macs2_container=check_tool('macs2')
if a.out.exists():p.error('Report exists')
commands={
 'samtools':['samtools','--version'], 'bedtools':['bedtools','--version'], 'bowtie2':['bowtie2','--version'],
 'python_plot_tracks':[sys.executable,'-c','import matplotlib,pyBigWig; print(matplotlib.__version__,pyBigWig.__version__)'],
 'R_analysis':['Rscript','-e','for(p in c("DESeq2","ChIPseeker","clusterProfiler","org.Hs.eg.db","jsonlite","digest")){if(!requireNamespace(p,quietly=TRUE))stop(p);cat(p,as.character(packageVersion(p)),"\\n")}'],
 'macs2_runtime':['macs2','callpeak','--help'],
}
if a.macs2_container:commands['macs2_runtime']=['apptainer','exec','--cleanenv','--no-home',str(a.macs2_container.resolve()),'macs2','callpeak','--help']
records=[]
for name,command in commands.items():
 try:
  result=subprocess.run(command,capture_output=True,text=True,timeout=120)
  records.append({'component':name,'command':command,'executable':shutil.which(command[0]),'state':'PASS' if result.returncode==0 else 'FAILED','returncode':result.returncode,'output':(result.stdout+result.stderr)[-4000:]})
 except (OSError,subprocess.TimeoutExpired) as error:
  records.append({'component':name,'command':command,'state':'FAILED','error':str(error)})
report={'state':'PASS' if all(r['state']=='PASS' for r in records) else 'FAILED','checks':records,
 'macs2_container_sha256':sha256(a.macs2_container) if a.macs2_container else None,
 'limitations':['Runtime probes do not replace official nf-core and synthetic fixture acceptance.','R and Python executable paths/versions are recorded; a system R run does not validate the locked R environment.']}
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2)+'\n')
print(report['state']);raise SystemExit(0 if report['state']=='PASS' else 1)
