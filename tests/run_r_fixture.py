"""Exercise real DESeq2 on small deterministic overdispersed synthetic counts."""
import csv,random,subprocess,tempfile,json,hashlib,os,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];rng=random.Random(513)
with tempfile.TemporaryDirectory() as d:
 p=Path(d);ids=['C1','C2','C3','K1','K2','K3']
 with (p/'counts.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['chrom','end','peak_id','start']+ids)
  for i in range(300):
   mu=rng.uniform(20,500);fold=3 if i<60 else 1
   v=[max(1,round(rng.gammavariate(8,mu*(fold if j>=3 else 1)/8))) for j in range(6)]
   w.writerow(['chr1',i*1000+200,f'p{i}',i*1000]+v)
 with (p/'meta.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','biological_sample_id','condition','replicates_confirmed'])
  for name in ids:w.writerow([name,name,'Control' if name[0]=='C' else 'KD','TRUE'])
 (p/'evidence.txt').write_text('Synthetic fixture: three independent simulated samples per condition; all inputs artificially generated for software validation only.')
 digest=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
 review={'reviewer':'automated synthetic test','reviewed_at':'2026-10-03','reason':'Synthetic fixture only','decisions':{k:True for k in ('upstream_accepted','peaks_accepted','replicates_confirmed','simple_design_accepted')},'inputs':{'counts_sha256':digest(p/'counts.tsv'),'metadata_sha256':digest(p/'meta.tsv')},'model':{'numerator':'KD','denominator':'Control','normalization':'conventional','design':'~ condition','alpha':.05,'abs_log2fc':1,'min_total_count':10},'evidence':[{'path':'evidence.txt','sha256':digest(p/'evidence.txt')}]}
 (p/'config').mkdir()
 cfg=json.loads((ROOT/'config/project.json').read_text());cfg['analysis']['replicates_confirmed']=True;cfg['qc'].update(upstream_accepted=True,peaks_accepted=True);cfg['nfcore_params']['use_control']=False
 with (p/'config/samples.tsv').open('w') as f:
  cols=['sample_id','biological_sample_id','library_id','unit_id','condition','group','replicate','role','control_group','fastq_1','fastq_2'];w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
  for i,name in enumerate(ids):
   condition='Control' if i<3 else 'KD'
   w.writerow(dict(sample_id=name,biological_sample_id=name,library_id=name,unit_id=name,condition=condition,group=condition,replicate=str(i%3+1),role='target',control_group='',fastq_1=name+'_R1.fastq.gz',fastq_2=name+'_R2.fastq.gz'))
 (p/'config/project.json').write_text(json.dumps(cfg));review['inputs']['config_sha256']=digest(p/'config/project.json')
 (p/'review.json').write_text(json.dumps(review))
 subprocess.run([sys.executable,str(ROOT/'scripts/run_differential.py'),'--config',str(p/'config/project.json'),'--counts',str(p/'counts.tsv'),'--metadata',str(p/'meta.tsv'),'--review',str(p/'review.json'),'--out',str(p/'model')],check=True)
 subprocess.run(['Rscript',str(ROOT/'scripts/visualize.R'),str(p/'model'),str(p/'figures')],check=True)
 z=list(csv.DictReader((p/'model/complete.tsv').open(),delimiter='\t'))
 assert len(z)==300 and sum(float(r['log2FoldChange'])>0 for r in z if int(r['peak_id'][1:])<60)>50
 assert (p/'figures/diagnostics.pdf').stat().st_size>0
 assert all(r['chrom']=='chr1' and int(r['end'])-int(r['start'])==200 for r in z)
 print('Synthetic DESeq2 and visualization PASS; contrast direction verified')
