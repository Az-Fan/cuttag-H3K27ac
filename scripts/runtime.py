"""Shared immutable run provenance for external-tool modules."""
import hashlib,json,os,subprocess,datetime
from pathlib import Path

def sha256(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def execute(command,folder,env=None):
 folder=Path(folder);folder.mkdir(parents=True,exist_ok=False)
 (folder/'command.json').write_text(json.dumps(command,indent=2))
 started=datetime.datetime.now(datetime.timezone.utc).isoformat()
 try:
  with (folder/'execution.log').open('w') as f:
   r=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,env=env)
  status={'state':'COMPUTATIONAL_PASS' if r.returncode==0 else 'FAILED','returncode':r.returncode,'started':started}
 except Exception as e:
  status={'state':'FAILED','message':str(e),'started':started}
  (folder/'status.json').write_text(json.dumps(status));raise
 (folder/'status.json').write_text(json.dumps(status))
 if r.returncode:raise RuntimeError('Command failed: '+str(command)+'; log: '+str(folder/'execution.log'))
 return r
