"""Read-only metadata replay of the two source designs; never reads sequencing data."""
import csv,json,importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];s=importlib.util.spec_from_file_location('c',ROOT/'scripts/cuttag.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
base=ROOT.parent;cfg=json.loads((ROOT/'config/project.json').read_text())
for project in ['ACLY-CUTTAG','ACLY-CUTTAG-batch2']:
 rows=[]
 with (base/project/'01_config/samplesheet.csv').open() as f:
  source=list(csv.DictReader(f))
 for i,r in enumerate(source):
  control=r['group']=='IgG';condition='IgG' if control else ('Control' if r['group'].startswith(('NC','siControl')) else 'KD')
  bio=r['group']+'_sample1' if project.endswith('batch2') else r['group']+'_'+r['replicate']
  rows.append(dict(sample_id=bio,biological_sample_id=bio,library_id=bio+'_lib',unit_id='unit'+str(i),condition=condition,group=r['group'],replicate='1' if project.endswith('batch2') else r['replicate'],role='control' if control else 'target',control_group=r['control'],fastq_1=r['fastq_1'],fastq_2=r['fastq_2']))
 cfg['analysis']['replicates_confirmed']=False
 result=m.validate(base/project,cfg,rows,True)
 assert not result['errors'],result['errors'];assert not result['formal_inference_eligible']
 if project.endswith('batch2'):assert result['biological_samples_per_condition']=={'KD':1,'Control':1}
 print(project,json.dumps({'counts':result['biological_samples_per_condition'],'formal_inference_eligible':False}))
