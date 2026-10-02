#!/usr/bin/env python3
"""Create inspectable release manifest and HTML navigation, checking links and required artifacts."""
import argparse,hashlib,html,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True,type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
if a.out.exists():p.error('Output exists')
spec=json.loads(a.manifest.read_text());records=[];errors=[]
for item in spec['artifacts']:
 path=Path(item['path']);path=path if path.is_absolute() else a.manifest.parent/path
 if not path.is_file():
  if item.get('required',True):errors.append('Missing/broken artifact '+str(path))
  records.append(dict(item,status='MISSING'));continue
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 records.append(dict(item,status='PRESENT',resolved=str(path.resolve()),bytes=path.stat().st_size,sha256=h.hexdigest()))
if errors:raise ValueError('; '.join(errors))
a.out.mkdir(parents=True)
(a.out/'manifest.json').write_text(json.dumps({'project':spec['project'],'limitations':spec.get('limitations',[]),'artifacts':records,'state':'REVIEW_REQUIRED'},indent=2,ensure_ascii=False))
body='<h1>'+html.escape(spec['project'])+'</h1><p>Release review required. Manifest preserves resolved source paths; package files separately for transfer.</p>'
body+='<ul>'+''.join('<li>'+html.escape(x)+'</li>' for x in spec.get('limitations',[]))+'</ul><table><tr><th>Artifact</th><th>Status</th><th>SHA256</th></tr>'
for r in records:body+='<tr><td>'+html.escape(r.get('label',r['path']))+'</td><td>'+r['status']+'</td><td>'+r.get('sha256','')+'</td></tr>'
(a.out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>CUTTag release review</title><style>body{font:16px sans-serif;max-width:1100px;margin:40px auto}td,th{padding:8px;border:1px solid #ddd;word-break:break-all}table{border-collapse:collapse}</style>'+body+'</table>')
