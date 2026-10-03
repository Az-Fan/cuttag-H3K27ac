"""Initialize -> reference check -> paired fragments -> consensus -> counts/QC/tracks -> portable delivery."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();a.out=a.out.resolve()
subprocess.run([sys.executable,str(ROOT/'scripts/init_project.py'),'--out',str(a.out),'--project-id','synthetic_acceptance'],check=True)
cfg=json.loads((a.out/'config/project.json').read_text());cfg['reference']['mito_name']='';cfg['resources']['cpus']=1
(a.out/'config/project.json').write_text(json.dumps(cfg))
external=a.out/'data/external'
(external/'genome.fa').write_text('>chr1\n'+'A'*1000+'\n')
(external/'genes.gtf').write_text('chr1\tfixture\tgene\t1\t1000\t.\t+\t.\tgene_id "g1";\n')
(external/'blacklist.bed').write_text('chr1\t900\t950\n')
subprocess.run(['samtools','faidx',str(external/'genome.fa')],check=True)
fixture=a.out/'fixture';fixture.mkdir()
lines=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:1000']
for i,start in enumerate((11,101,501)):
 for flag,offset,mate,length in [(99,0,20,30),(147,20,0,-30)]:
  lines.append('\t'.join(map(str,['p'+str(i),flag,'chr1',start+offset,60,'10M','=',start+mate,length,'A'*10,'I'*10])))
(fixture/'x.sam').write_text('\n'.join(lines)+'\n')
subprocess.run(['samtools','view','-b','-o',str(fixture/'x.bam'),str(fixture/'x.sam')],check=True)
(fixture/'y.sam').write_text('\n'.join(lines).replace('p0','q0').replace('p1','q1').replace('p2','q2')+'\n')
subprocess.run(['samtools','view','-b','-o',str(fixture/'y.bam'),str(fixture/'y.sam')],check=True)
(fixture/'peaks.bed').write_text('chr1\t0\t50\nchr1\t100\t150\n')
(fixture/'artifacts.tsv').write_text('sample_id\tbam\tpeaks_bed\tbam_policy\nControl1\tx.bam\tpeaks.bed\tduplicates_retained\nKD1\ty.bam\tpeaks.bed\tduplicates_retained\n')
# Independently named synthetic pairs; fixture is not biological evidence.
plan=a.out/'results/downstream'
subprocess.run([sys.executable,str(a.out/'scripts/plan_downstream.py'),'--config',str(a.out/'config/project.json'),'--artifacts',str(fixture/'artifacts.tsv'),'--out',str(plan)],check=True)
assert json.loads((plan/'handoff.json').read_text())['formal_inference_eligible'] is False
command=[sys.executable,str(a.out/'scripts/workflow.py'),'--workflow',str(plan/'workflow.json'),'--out',str(a.out/'results/execution')]
env=os.environ.copy();env['MPLCONFIGDIR']=str(a.out/'shared_cache/matplotlib')
subprocess.run(command,check=True,env=env)
subprocess.run(command+['--resume'],check=True,env=env)
assert all(r['resumed'] for r in json.loads((a.out/'results/execution/workflow_status.json').read_text())['steps'])
release={'project':'synthetic_acceptance','limitations':['Artificial paired fragments; no biological conclusions.'], 'artifacts':[{'path':'results/downstream/counts.tsv','destination':'tables/counts.tsv'},{'path':'results/downstream/qc_atlas','destination':'qc'},{'path':'results/downstream/tracks','destination':'tracks'},{'path':'results/downstream/handoff.json','destination':'handoff.json'}]}
(a.out/'release.json').write_text(json.dumps(release))
subprocess.run([sys.executable,str(a.out/'scripts/finalize.py'),'--manifest',str(a.out/'release.json'),'--out',str(a.out/'results/delivery')],check=True)
shutil.move(str(a.out/'results/delivery'),str(a.out/'moved_delivery'))
subprocess.run([sys.executable,str(a.out/'scripts/verify_release.py'),'--package',str(a.out/'moved_delivery')],check=True)
print('PASS: initialized project, planned DAG, real tools, automatic QC, exact resume, relocated delivery; inference blocked for 1 vs 1')
