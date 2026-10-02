#!/usr/bin/env python3
"""Portable CUT&Tag controller; configuration only uses Python's standard library."""
import argparse, csv, datetime, gzip, hashlib, itertools, json, os, re, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n")
def read_rows(path):
    with path.open() as f: return list(csv.DictReader(f, delimiter="\t"))
def resolve(root, value):
    p = Path(value); return p if p.is_absolute() else root / p
def load(path):
    path = path.resolve(); root = path.parent.parent
    cfg = json.loads(path.read_text()); rows = read_rows(resolve(root, cfg["samples"]))
    return root, cfg, rows
def validate(root, cfg, rows, missing=False):
    errors=[]; warnings=[]; bios={}; mapping={}; units=set(); files=set(); libs={}; samples={}
    required={'sample_id','biological_sample_id','library_id','unit_id','condition','group','replicate','role','control_group','fastq_1','fastq_2'}
    if not rows: errors.append('Sample table is empty')
    controls={r.get('group') for r in rows if r.get('role')=='control'}
    for r in rows:
        if not required.issubset(r): errors.append('Missing sample table columns'); break
        if any(not r[k] for k in required-{'control_group'}): errors.append('Empty required sample fields')
        if r['role'] not in ('target','control'): errors.append('Unknown sample role')
        if not r['replicate'].isdigit() or int(r['replicate'])<1: errors.append('Replicate must be positive integer')
        for key in ('sample_id','group'):
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',r[key]): errors.append('Unsafe identifier: '+r[key])
        if r['unit_id'] in units: errors.append('Duplicate sequencing unit '+r['unit_id'])
        units.add(r['unit_id'])
        identity=(r['biological_sample_id'],r['condition'],r['group'],r['replicate'],r['role'],r['control_group'])
        for field, registry, val in [('sample',samples,r['sample_id']),('library',libs,r['library_id'])]:
            if val in registry and registry[val]!=identity: errors.append('Conflicting '+field+' mapping: '+val)
            registry[val]=identity
        b=r['biological_sample_id']; key=(r['group'],r['replicate'])
        if b in bios and bios[b]!=(r['condition'],key,r['role']): errors.append('Biological sample assigned multiple conditions/replicates: '+b)
        bios[b]=(r['condition'],key,r['role'])
        if key in mapping and mapping[key]!=b: errors.append('Distinct biological samples share group/replicate')
        mapping[key]=b
        if r['role']=='target' and cfg['nfcore_params'].get('use_control') and r['control_group'] not in controls: errors.append('Missing matched control for '+r['sample_id'])
        for k in ('fastq_1','fastq_2'):
            p=resolve(root,r[k]).resolve()
            if str(p) in files: errors.append('FASTQ reused: '+str(p))
            files.add(str(p))
            if not str(p).endswith(('.fastq.gz','.fq.gz')): errors.append('Expected gzipped FASTQ: '+str(p))
            if not p.is_file(): (warnings if missing else errors).append('Missing '+str(p))
    for key in ('fasta','gtf','blacklist'):
        value=cfg['reference'].get(key)
        if not value: errors.append('Reference requires '+key)
        elif not resolve(root,value).is_file(): (warnings if missing else errors).append('Missing reference '+value)
    sp=cfg['spikein']
    if sp['enabled']:
        if sp['identity']=='unknown' or not sp['fasta']: errors.append('Spike-in identity and FASTA required')
        elif not resolve(root,sp['fasta']).is_file(): (warnings if missing else errors).append('Missing spike-in FASTA')
    counts={c:len({r['biological_sample_id'] for r in rows if r['role']=='target' and r['condition']==c}) for c in (cfg['analysis']['numerator'],cfg['analysis']['denominator'])}
    inference=not errors and cfg['analysis']['replicates_confirmed'] and all(n>=2 for n in counts.values())
    if not inference: warnings.append('Formal differential inference unavailable: confirm independent replicates and >=2 samples per contrast group')
    return {'errors':errors,'warnings':warnings,'biological_samples_per_condition':counts,'formal_inference_eligible':bool(inference)}
def fastq_records(path):
    with gzip.open(path,'rt') as f:
        while True:
            h=f.readline()
            if not h: break
            s,p,q=f.readline().rstrip(),f.readline().rstrip(),f.readline().rstrip()
            if not h.startswith('@') or not p.startswith('+') or not s or len(s)!=len(q): raise ValueError('Malformed FASTQ '+str(path))
            yield re.sub(r'/[12]$', '', h.split()[0]),len(s)
def inventory(root, rows, out):
    result=[]
    for r in rows:
        p1,p2=[resolve(root,r[k]) for k in ('fastq_1','fastq_2')]; n=0; lengths=set()
        for a,b in itertools.zip_longest(fastq_records(p1),fastq_records(p2)):
            if a is None or b is None or a[0]!=b[0]: raise ValueError('Mate mismatch in '+r['unit_id'])
            n+=1; lengths.update((a[1],b[1]))
        if not n: raise ValueError('Empty FASTQ '+r['unit_id'])
        for p in (p1,p2):
            h=hashlib.sha256()
            with p.open('rb') as f:
                for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
            result.append({'unit':r['unit_id'],'path':str(p.resolve()),'bytes':p.stat().st_size,'sha256':h.hexdigest(),'pairs':n,'read_lengths':sorted(lengths)})
    dump(out,result)
def plan(root,cfg,rows,stage,out):
    out.mkdir(parents=True,exist_ok=False)
    with (out/'samplesheet.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['group','replicate','fastq_1','fastq_2','control'])
        for r in rows: w.writerow([r['group'],r['replicate'],str(resolve(root,r['fastq_1']).resolve()),str(resolve(root,r['fastq_2']).resolve()),r['control_group']])
    params=dict(cfg['nfcore_params']); ref=cfg['reference']
    params.update(input=str(out/'samplesheet.csv'),outdir=str(out/'output'),igenomes_ignore=True,fasta=str(resolve(root,ref['fasta']).resolve()),gtf=str(resolve(root,ref['gtf']).resolve()),blacklist=str(resolve(root,ref['blacklist']).resolve()),macs_gsize=ref['macs_gsize'],mito_name=ref['mito_name'],max_cpus=cfg['resources']['cpus'],max_memory=str(cfg['resources']['memory_gb'])+'.GB',normalisation_mode='CPM',only_filtering=stage=='alignment',publish_dir_mode='copy')
    if cfg['spikein']['enabled']:
        params.update(spikein_fasta=str(resolve(root,cfg['spikein']['fasta']).resolve()))
        if stage=='production': params['normalisation_mode']='Spikein'
    dump(out/'params.json',params);dump(out/'project.snapshot.json',cfg);dump(out/'samples.snapshot.json',rows)
    cmd=['nextflow','run',cfg['pipeline']['name'],'-r',cfg['pipeline']['version'],'-profile',cfg['pipeline']['profile'],'-params-file',str(out/'params.json'),'-c',str(root/'config/profiles/local.config'),'-work-dir',str(root/'work'/out.name),'-with-trace',str(out/'trace.tsv'),'-with-report',str(out/'report.html'),'-with-timeline',str(out/'timeline.html')]
    dump(out/'command.json',cmd);return cmd

def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['validate','inventory','plan','run','test']);p.add_argument('--config',type=Path,default=ROOT/'config/project.json');p.add_argument('--allow-missing',action='store_true');p.add_argument('--stage',choices=['alignment','production'],default='alignment');p.add_argument('--run-id');a=p.parse_args()
    root,cfg,rows=load(a.config);report=validate(root,cfg,rows,a.allow_missing or a.command=='plan' or a.command=='test')
    if a.command=='validate': print(json.dumps(report,indent=2,ensure_ascii=False));return int(bool(report['errors']))
    if report['errors']: raise ValueError('; '.join(report['errors']))
    if a.command=='inventory': inventory(root,rows,root/'results/qc/fastq_inventory.json');return 0
    runid=a.run_id or datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    if not re.fullmatch(r'[A-Za-z0-9_-]+',runid): raise ValueError('Invalid run ID')
    out=(root/'results/runs'/runid).resolve()
    if a.command=='run':
        if a.allow_missing: raise ValueError('--allow-missing cannot be used for execution')
        if a.stage=='production' and not cfg['qc']['upstream_accepted']: raise ValueError('Accept upstream QC before production')
        if a.stage=='production' and cfg['spikein']['enabled'] and not cfg['spikein']['calibration_accepted']: raise ValueError('Accept spike-in calibration before production')
        if shutil.disk_usage(root).free < cfg['resources']['min_free_gb']*1024**3: raise ValueError('Insufficient free storage')
    if a.command=='test':
        out.mkdir(parents=True,exist_ok=False)
        cmd=['nextflow','run',cfg['pipeline']['name'],'-r',cfg['pipeline']['version'],'-profile','test,'+cfg['pipeline']['profile'],'--outdir',str(out/'output'),'-work-dir',str(root/'work'/runid)]
        dump(out/'command.json',cmd)
    else: cmd=plan(root,cfg,rows,a.stage,out)
    print(json.dumps({'run_dir':str(out),'command':cmd,'warnings':report['warnings']},indent=2))
    if a.command in ('run','test'):
        env=os.environ.copy();cache=root/'shared_cache';cache.mkdir(exist_ok=True)
        env.update(NXF_HOME=str(cache/'nextflow'),NXF_APPTAINER_CACHEDIR=str(cache/'apptainer'),APPTAINER_CACHEDIR=str(cache/'apptainer'))
        dump(out/'status.json',{'state':'RUNNING','started':datetime.datetime.now(datetime.timezone.utc).isoformat()})
        try:
            with (out/'execution.log').open('w') as f: result=subprocess.run(cmd,cwd=out,env=env,stdout=f,stderr=subprocess.STDOUT)
            dump(out/'status.json',{'state':'COMPUTATIONAL_PASS' if result.returncode==0 else 'FAILED','returncode':result.returncode})
            return result.returncode
        except Exception as e:
            dump(out/'status.json',{'state':'FAILED','message':str(e)});raise
    return 0
if __name__=='__main__':
    try: sys.exit(main())
    except (ValueError,OSError,KeyError) as e: print('ERROR: '+str(e),file=sys.stderr);sys.exit(2)
