"""Real fragment-policy and three MACS2 background strategies on a fixed simulated experiment."""
import argparse
import csv
import json
import random
import subprocess
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--out', type=Path, required=True)
p.add_argument('--macs2-container',type=Path)
a = p.parse_args()
a.out = a.out.resolve()
a.out.mkdir(parents=True, exist_ok=False)
rng = random.Random(749)

def bam(name, starts):
    sam = a.out/(name+'.sam')
    with sam.open('w') as f:
        f.write('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000000\n')
        for i, (start, mapq, duplicate) in enumerate(starts):
            for flag, offset, mate, tlen in [(99,0,100,150),(147,100,0,-150)]:
                f.write('\t'.join(map(str, [name+str(i),flag+(1024 if duplicate else 0),'chr1',start+offset,mapq,'50M','=',start+mate,tlen,'ACGT'*12+'AC','I'*50]))+'\n')
    path = a.out/(name+'.bam')
    subprocess.run(['samtools','view','-b','-o',str(path),str(sam)], check=True)
    return path

probe = bam('policy', [(100,60,False),(300,60,True),(500,25,False),(700,10,False)])
policies = []
for q, remove, expected in [(20,False,3),(20,True,2),(30,False,2),(30,True,1)]:
    dest = a.out/('mapq%d_%s' % (q,'dedup' if remove else 'keep'))
    cmd = [sys.executable,str(ROOT/'scripts/fragments.py'),'--bam',str(probe),'--out',str(dest),'--mapq',str(q)]
    if remove: cmd += ['--remove-duplicates']
    subprocess.run(cmd,check=True)
    n=json.loads((dest/'summary.json').read_text())['fragments']
    assert n==expected,(q,remove,n)
    policies.append({'mapq':q,'remove_duplicates':remove,'expected':expected,'observed':n,'state':'PASS'})
target = bam('target', [(rng.randrange(1,999800),60,False) for _ in range(2000)] + [(center+rng.randrange(-100,100),60,False) for center in (10000,30000,50000) for _ in range(1000)])
control = bam('control', [(rng.randrange(1,999800),60,False) for _ in range(500)])
manifest = a.out/'peaks.tsv'
manifest.write_text('sample_id\ttarget_bam\tcontrol_bam\tbam_policy\ns1\ttarget.bam\tcontrol.bam\tduplicates_retained\n')
subprocess.run([sys.executable,str(ROOT/'scripts/peak_diagnostics.py'),'--manifest',str(manifest),'--out',str(a.out/'peaks'),'--gsize','1000000','--duplicate-modes','all','auto','--execute']+(['--macs2-container',str(a.macs2_container.resolve())] if a.macs2_container else []),check=True,cwd='/tmp')
reports=json.loads((a.out/'peaks/diagnostics.json').read_text())
reports=reports['candidates']
assert len(reports)==12 and all(r['state']=='COMPUTATIONAL_PASS' for r in reports)
assert {r['duplicate_mode'] for r in reports}=={'all','auto'}
assert len(json.loads((a.out/'peaks/diagnostics.json').read_text())['artifact_manifests'])==12
# Every method must recover at least one planted enriched locus, not merely exit zero.
for report in reports:
    path=a.out/'peaks'/('s1_'+report['strategy'])/('s1_peaks.'+report['shape']+'Peak')
    peaks=[line.split() for line in path.read_text().splitlines()]
    assert any(int(r[1])<10050 and int(r[2])>10000 for r in peaks),report
with (a.out/'comparison.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','strategy','peaks_bed'])
    for r in reports:w.writerow(['s1',r['strategy'],'peaks/s1_'+r['strategy']+'/s1_peaks.'+r['shape']+'Peak'])
subprocess.run([sys.executable,str(ROOT/'scripts/compare_peaks.py'),'--manifest',str(a.out/'comparison.tsv'),'--out',str(a.out/'comparison')],check=True)
(a.out/'acceptance.json').write_text(json.dumps({'state':'PASS','scope':'synthetic software correctness, not real-sample strategy acceptance','fragment_policies':policies,'peak_methods':reports},indent=2))
print('PASS: 4 MAPQ/duplicate policies and 6 MACS2 narrow/broad × background strategies; planted signal recovered')
