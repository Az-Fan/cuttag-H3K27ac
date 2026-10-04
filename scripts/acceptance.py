#!/usr/bin/env python3
"""Repeatable local acceptance suite; retain every subprocess log and exact source hashes."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from decisions import sha256
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);p.add_argument('--quick',action='store_true');a=p.parse_args();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
checks=[('contracts',[sys.executable,'-m','unittest','discover','-s','tests','-v'])]
if not a.quick:
 checks += [
  ('runtime',[sys.executable,'scripts/doctor.py','--out',str(a.out/'runtime.json')]),
  ('fragments',[sys.executable,'tests/run_fragment_fixture.py']),
  ('spikein_fragments',[sys.executable,'tests/run_spikein_fixture.py']),
  ('project',[sys.executable,'tests/run_project_fixture.py','--out',str(a.out/'project')]),
  ('options',[sys.executable,'tests/run_options_fixture.py','--out',str(a.out/'options')]),
  ('R_models_ORA',[sys.executable,'tests/run_r_fixture.py','--out',str(a.out/'R')]),
  ('annotation',[sys.executable,'tests/run_annotation_fixture.py','--out',str(a.out/'annotation')]),
  ('motif',[sys.executable,'tests/run_motif_fixture.py','--out',str(a.out/'motif')]),
 ]
env=os.environ.copy();env['MPLCONFIGDIR']=str(a.out/'matplotlib_cache')
records=[]
sources={str(path.relative_to(ROOT)):sha256(path) for directory in ('scripts','tests','config') for path in (ROOT/directory).rglob('*') if path.is_file() and '__pycache__' not in path.parts}
sources['pixi.lock']=sha256(ROOT/'pixi.lock')
for name,command in checks:
 with (a.out/(name+'.log')).open('w') as log:
  result=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
 records.append({'check':name,'command':command,'returncode':result.returncode,'state':'PASS' if result.returncode==0 else 'FAILED'})
 print(name+': '+records[-1]['state'],flush=True)
 (a.out/'status.json').write_text(json.dumps({'state':'RUNNING','checks':records},indent=2))
passed=all(r['state']=='PASS' for r in records)
(a.out/'status.json').write_text(json.dumps({'state':'PASS' if passed else 'FAILED','scope':'contracts only' if a.quick else 'local synthetic integration + runtime',
 'checks':records,'source_sha256':sources,'limitations':['Official nf-core test is separate: scripts/cuttag.py test.','Synthetic tests do not approve scientific decisions or validate biological performance on a new experiment.']},indent=2))
raise SystemExit(0 if passed else 1)
