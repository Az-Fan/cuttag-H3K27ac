"""Review records are explicit decisions with content-bound evidence, not software success flags."""
import hashlib,json
from pathlib import Path

def known_text(value):
    return isinstance(value,str) and bool(value.strip()) and value.strip().lower() not in {'unknown','unspecified','not specified','na','n/a','待确认','未确认'}

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def require_review(path,required_decisions,inputs,contrast=None):
    path=Path(path).resolve();review=json.loads(path.read_text())
    for key in ('reviewer','reviewed_at','reason'):
        if not isinstance(review.get(key),str) or not review[key].strip():raise ValueError('Review requires '+key)
    for key in required_decisions:
        if review.get('decisions',{}).get(key) is not True:raise ValueError('Review has not accepted '+key)
    for name,file in inputs.items():
        if review.get('inputs',{}).get(name)!=sha256(file):raise ValueError('Review input missing or changed: '+name)
    if contrast:
        for key,value in contrast.items():
            if review.get('model',{}).get(key)!=value:raise ValueError('Review model mismatch: '+key)
    evidence=review.get('evidence',[])
    if not evidence:raise ValueError('Review requires evidence files')
    for item in evidence:
        file=Path(item['path']);file=file if file.is_absolute() else path.parent/file
        if item.get('sha256')!=sha256(file):raise ValueError('Review evidence missing or changed: '+str(file))
    return review
