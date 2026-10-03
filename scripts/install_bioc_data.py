#!/usr/bin/env python3
"""Complete missing Bioconda data-package post-link installs in the project environment."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
from decisions import sha256

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--prefix',type=Path,default=ROOT/'.pixi/envs/analysis');p.add_argument('--out',required=True,type=Path);a=p.parse_args()
prefix=a.prefix.resolve()
if a.out.exists():p.error('New log directory required')
registry=json.loads((prefix/'share/bioconductor-data-packages/dataURLs.json').read_text())
packages=[]
for script in sorted((prefix/'bin').glob('.bioconductor-*-post-link.sh')):
    match=re.search(r'installBiocDataPackage.sh "([a-z0-9._-]+)"',script.read_text())
    if match:packages.append(match[1])
a.out.mkdir(parents=True);records=[]
for key in packages:
    spec=registry[key];name=spec['fn'].split('_')[0]
    if not re.fullmatch(r'[A-Za-z][A-Za-z0-9.]+',name):raise ValueError('Invalid R package name')
    probe=[str(prefix/'bin/Rscript'),'-e','quit(status=if(requireNamespace("'+name+'",quietly=TRUE))0 else 1)']
    if subprocess.run(probe,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode==0:
        records.append({'package':name,'state':'ALREADY_LOADABLE'});continue
    tarball=a.out/spec['fn'];downloaded=False
    for url in spec['urls']:
        with (a.out/(key+'.download.log')).open('a') as log:
            result=subprocess.run(['curl','-L','--fail','--connect-timeout','20','--max-time','300','--output',str(tarball),url],stdout=log,stderr=subprocess.STDOUT)
        if result.returncode==0 and hashlib.md5(tarball.read_bytes()).hexdigest()==spec['md5']:
            downloaded=True;break
    if not downloaded:raise RuntimeError('Data package download/checksum failed: '+name)
    command=[str(prefix/'bin/R'),'CMD','INSTALL','--library='+str(prefix/'lib/R/library'),str(tarball.resolve())]
    with (a.out/(key+'.install.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    subprocess.run(probe,check=True)
    records.append({'package':name,'state':'INSTALLED_AND_LOADABLE','url':url,'upstream_md5':spec['md5'],'sha256':sha256(tarball),'command':command})
    (a.out/'status.json').write_text(json.dumps({'state':'RUNNING','packages':records},indent=2))
(a.out/'status.json').write_text(json.dumps({'state':'PASS','packages':records,'registry_sha256':sha256(prefix/'share/bioconductor-data-packages/dataURLs.json')},indent=2))
print('Bioconductor data dependencies PASS')
