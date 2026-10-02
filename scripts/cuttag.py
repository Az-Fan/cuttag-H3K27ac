#!/usr/bin/env python3
"""Portable CUT&Tag controller; configuration only uses Python's standard library."""
import argparse, csv, datetime, gzip, hashlib, itertools, json, os, re, shutil, subprocess, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from decisions import require_review, sha256, known_text
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
    errors=[]; warnings=[]; bios={}; mapping={}; units=set(); files=set(); libs={}; samples={}; group_definitions={}; bio_sample_ids={}
    required={'sample_id','biological_sample_id','library_id','unit_id','condition','group','replicate','role','control_group','fastq_1','fastq_2'}
    for section,keys in {'spikein':['enabled','equal_amount_confirmed','calibration_accepted'], 'analysis':['replicates_confirmed'], 'qc':['upstream_accepted','peaks_accepted'], 'nfcore_params':['use_control','dedup_target_reads','remove_linear_duplicates']}.items():
        for key in keys:
            if type(cfg.get(section,{}).get(key)) is not bool:errors.append(section+'.'+key+' must be JSON boolean')
    analysis=cfg['analysis']
    if analysis['numerator']==analysis['denominator']:errors.append('Contrast numerator and denominator must differ')
    if analysis.get('design')!='~ condition':errors.append('Only ~ condition supported; batch/paired designs require an explicit model implementation')
    if analysis.get('normalization') not in ('conventional','spikein'):errors.append('Unsupported normalization')
    for key,low,high in [('alpha',0,1),('abs_log2fc',0,float('inf')),('min_total_count',0,float('inf'))]:
        val=analysis.get(key)
        if type(val) not in (int,float) or not __import__('math').isfinite(val) or val<low or val>=high or (key=='alpha' and val==0):errors.append('Invalid analysis.'+key)
    if not rows: errors.append('Sample table is empty')
    controls={r.get('group') for r in rows if r.get('role')=='control'}
    for r in rows:
        if not required.issubset(r) or any(r[k] is None for k in required):
            return {'errors':['Missing sample table columns/values'],'warnings':[],'biological_samples_per_condition':{},'formal_inference_eligible':False}
        if any(not r[k] for k in required-{'control_group'}): errors.append('Empty required sample fields')
        if r['role'] not in ('target','control'): errors.append('Unknown sample role')
        if not r['replicate'].isdigit() or int(r['replicate'])<1 or str(int(r['replicate']))!=r['replicate']: errors.append('Replicate must be canonical positive integer')
        for key in ('sample_id','group'):
            if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*',r[key]): errors.append('Unsafe identifier: '+r[key])
        if r['unit_id'] in units: errors.append('Duplicate sequencing unit '+r['unit_id'])
        units.add(r['unit_id'])
        identity=(r['biological_sample_id'],r['condition'],r['group'],r['replicate'],r['role'],r['control_group'])
        for field, registry, val in [('sample',samples,r['sample_id']),('library',libs,r['library_id'])]:
            if val in registry and registry[val]!=identity: errors.append('Conflicting '+field+' mapping: '+val)
            registry[val]=identity
        group_identity=(r['condition'],r['role'],r['control_group'])
        if r['group'] in group_definitions and group_definitions[r['group']]!=group_identity:errors.append('Pipeline group mixes conditions/roles/control strategy')
        group_definitions[r['group']]=group_identity
        b=r['biological_sample_id']; key=(r['group'],r['replicate'])
        if b in bio_sample_ids and bio_sample_ids[b]!=r['sample_id']:errors.append('Biological sample has multiple statistical sample IDs')
        bio_sample_ids[b]=r['sample_id']
        if b in bios and bios[b]!=(r['condition'],key,r['role']): errors.append('Biological sample assigned multiple conditions/replicates: '+b)
        bios[b]=(r['condition'],key,r['role'])
        if key in mapping and mapping[key]!=b: errors.append('Distinct biological samples share group/replicate')
        mapping[key]=b
        if r['role']=='control' and r['control_group']:errors.append('Control cannot name another control group')
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
    if sp['enabled'] is True:
        if not known_text(sp['identity']) or not isinstance(sp['fasta'],str) or not sp['fasta']: errors.append('Spike-in identity and FASTA required')
        elif not resolve(root,sp['fasta']).is_file(): (warnings if missing else errors).append('Missing spike-in FASTA')
        if sp.get('alignment_profile') and not resolve(root,sp['alignment_profile']).is_file():(warnings if missing else errors).append('Missing spike-in alignment profile')
    counts={c:len({r['biological_sample_id'] for r in rows if r['role']=='target' and r['condition']==c}) for c in (cfg['analysis']['numerator'],cfg['analysis']['denominator'])}
    if cfg['analysis']['normalization']=='spikein' and sp['enabled'] is not True:errors.append('Spike-in analysis requires enabled spike-in')
    inference=not errors and cfg['analysis']['replicates_confirmed'] is True and all(n>=2 for n in counts.values())
    if any(n==0 for n in counts.values()):errors.append('Missing contrast group in target samples')
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
        params['normalisation_mode']='Spikein'  # also align spike-in in filtering-only stage
        params['spikein_bowtie2']=None  # build index from explicit FASTA; never inherit a default index
    dump(out/'params.json',params);dump(out/'project.snapshot.json',cfg);dump(out/'samples.snapshot.json',rows)
    cmd=['nextflow','run',cfg['pipeline']['name'],'-r',cfg['pipeline']['version'],'-profile',cfg['pipeline']['profile'],'-params-file',str(out/'params.json'),'-c',str(root/'config/profiles/local.config'),'-work-dir',str(root/'work'/out.name),'-with-trace',str(out/'trace.tsv'),'-with-report',str(out/'report.html'),'-with-timeline',str(out/'timeline.html')]
    if cfg['spikein']['enabled'] and cfg['spikein'].get('alignment_profile'):
        cmd+=['-c',str(resolve(root,cfg['spikein']['alignment_profile']).resolve())]
    dump(out/'command.json',cmd);return cmd

def production_gate(root,cfg,config_path):
    if cfg['qc']['upstream_accepted'] is not True:raise ValueError('Accept upstream QC before production')
    decisions=['upstream_accepted','control_strategy_accepted']
    inputs={'config_sha256':config_path,'samples_sha256':resolve(root,cfg['samples']),'fastq_inventory_sha256':root/'results/qc/fastq_inventory.json'}
    sp=cfg['spikein']
    if sp['enabled']:
        if sp['calibration_accepted'] is not True or sp['equal_amount_confirmed'] is not True or not known_text(sp['added_at']):raise ValueError('Spike-in needs equal-input confirmation, known addition stage and accepted calibration')
        if not known_text(sp.get('calibration_scope')):raise ValueError('Record experimental scope of spike-in calibration')
        decisions+=['spikein_counting_accepted','spikein_calibration_accepted']
        if sp.get('alignment_profile'):inputs['alignment_profile_sha256']=resolve(root,sp['alignment_profile'])
    review=cfg['qc'].get('production_review')
    if not review:raise ValueError('Production requires a documented QC/control-strategy review')
    return require_review(resolve(root,review),decisions,inputs)

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
        test_id=cfg['qc'].get('official_test_run_id','')
        if not re.fullmatch(r'[A-Za-z0-9_-]+',test_id): raise ValueError('Record a successful official_test_run_id before real execution')
        test_dir=root/'results/runs'/test_id
        if not (test_dir/'status.json').is_file() or json.loads((test_dir/'status.json').read_text()).get('state')!='COMPUTATIONAL_PASS': raise ValueError('Official test has not passed')
        test_cmd=json.loads((test_dir/'command.json').read_text())
        if test_cmd[:3]!=['nextflow','run',cfg['pipeline']['name']] or test_cmd[test_cmd.index('-profile')+1]!='test,'+cfg['pipeline']['profile'] or test_cmd[test_cmd.index('-r')+1]!=cfg['pipeline']['version']: raise ValueError('Test profile/version mismatch')
        if a.allow_missing: raise ValueError('--allow-missing cannot be used for execution')
        if a.stage=='production':production_gate(root,cfg,a.config)
        if shutil.disk_usage(root).free < cfg['resources']['min_free_gb']*1024**3: raise ValueError('Insufficient free storage')
    if a.command=='test':
        out.mkdir(parents=True,exist_ok=False)
        cmd=['nextflow','run',cfg['pipeline']['name'],'-r',cfg['pipeline']['version'],'-profile','test,'+cfg['pipeline']['profile'],'--outdir',str(out/'output'),'-work-dir',str(root/'work'/runid)]
        dump(out/'command.json',cmd)
    else: cmd=plan(root,cfg,rows,a.stage,out)
    print(json.dumps({'run_dir':str(out),'command':cmd,'warnings':report['warnings']},indent=2))
    if a.command=='run':
        try:
            inventory(root,rows,out/'fastq_inventory.json')
            if a.stage=='production':
                require_review(resolve(root,cfg['qc']['production_review']),['upstream_accepted'],{'fastq_inventory_sha256':out/'fastq_inventory.json'})
        except Exception as e:
            dump(out/'status.json',{'state':'FAILED','stage':'input_integrity','message':str(e)});raise
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
