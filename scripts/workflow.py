#!/usr/bin/env python3
"""Sequential explicit workflow runner: each step is a recorded argv array, no shell interpolation."""
import argparse,json,re,datetime,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--workflow',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
if a.out.exists():p.error('Run output already exists')
spec=json.loads(a.workflow.read_text())
if not spec.get('steps'):p.error('Workflow has no steps')
ids=[s['id'] for s in spec['steps']]
if len(ids)!=len(set(ids)):p.error('Workflow step IDs must be unique')
a.out.mkdir(parents=True);states={};summary=[]
for step in spec['steps']:
 sid=step['id']
 if not re.fullmatch(r'[A-Za-z0-9_-]+',sid):raise ValueError('Unsafe step ID')
 unmet=[x for x in step.get('depends_on',[]) if states.get(x)!='COMPUTATIONAL_PASS']
 if unmet:
  states[sid]='BLOCKED';summary.append({'id':sid,'state':'BLOCKED','unmet':unmet});continue
 command=step['command']
 if not isinstance(command,list) or not command or not all(isinstance(x,str) for x in command):raise ValueError('Step command must be argv array')
 directory=a.out/sid;directory.mkdir()
 (directory/'command.json').write_text(json.dumps(command,indent=2));started=datetime.datetime.now(datetime.timezone.utc).isoformat()
 try:
  with (directory/'execution.log').open('w') as log:r=subprocess.run(command,cwd=spec.get('working_directory'),stdout=log,stderr=subprocess.STDOUT)
  state='COMPUTATIONAL_PASS' if r.returncode==0 else 'FAILED';detail={'returncode':r.returncode}
 except Exception as e:state='FAILED';detail={'message':str(e)}
 states[sid]=state;record={'id':sid,'state':state,'started':started,**detail};summary.append(record);(directory/'status.json').write_text(json.dumps(record,indent=2))
 if state=='FAILED' and not step.get('continue_on_error',False):
  for later in spec['steps'][len(summary):]:
   states[later['id']]='BLOCKED';summary.append({'id':later['id'],'state':'BLOCKED','reason':'upstream failure'})
  break
(a.out/'workflow_status.json').write_text(json.dumps({'state':'COMPUTATIONAL_PASS' if all(v=='COMPUTATIONAL_PASS' for v in states.values()) else 'REVIEW_REQUIRED','steps':summary},indent=2))
if any(v in ('FAILED','BLOCKED') for v in states.values()):raise SystemExit(1)
