"""Exercise real DESeq2 on small deterministic overdispersed synthetic counts."""
import csv,random,subprocess,tempfile,json,hashlib,os,sys,argparse,contextlib,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];rng=random.Random(513)
parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path);args=parser.parse_args()
if args.out:args.out.mkdir(parents=True,exist_ok=False)
with (contextlib.nullcontext(str(args.out.resolve())) if args.out else tempfile.TemporaryDirectory()) as d:
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

 # Spike-in branch: explicitly synthetic calibration, known factors and expected normalized values.
 factors=[1,2,3,1,2,3]
 with (p/'spikein.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','spikein_fragments','calibration_accepted','calibration_group'])
  for name,factor in zip(ids,factors):w.writerow([name,1000*factor,'TRUE','prepA'])
 cfg['analysis']['normalization']='spikein'
 cfg['spikein'].update(enabled=True,identity='synthetic_lambda',fasta='synthetic_lambda.fa',equal_amount_confirmed=True,added_at='synthetic_preparation',calibration_accepted=True,calibration_scope='synthetic known factors only')
 (p/'config/project.json').write_text(json.dumps(cfg))
 review['inputs'].update(config_sha256=digest(p/'config/project.json'),spikein_sha256=digest(p/'spikein.tsv'))
 review['model']['normalization']='spikein';review['spikein']=cfg['spikein']
 review['decisions'].update(spikein_counting_accepted=True,spikein_calibration_accepted=True)
 (p/'review_spikein.json').write_text(json.dumps(review))
 subprocess.run([sys.executable,str(ROOT/'scripts/run_differential.py'),'--config',str(p/'config/project.json'),'--counts',str(p/'counts.tsv'),'--metadata',str(p/'meta.tsv'),'--review',str(p/'review_spikein.json'),'--spikein',str(p/'spikein.tsv'),'--out',str(p/'model_spikein')],check=True)
 spike_rows=list(csv.DictReader((p/'spikein.tsv').open(),delimiter='\t'))
 spike_rows[-1]['calibration_group']='prepB'
 with (p/'spikein.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(spike_rows[0]),delimiter='\t');w.writeheader();w.writerows(spike_rows)
 blocked=subprocess.run([sys.executable,str(ROOT/'scripts/run_differential.py'),'--config',str(p/'config/project.json'),'--counts',str(p/'counts.tsv'),'--metadata',str(p/'meta.tsv'),'--review',str(p/'review_spikein.json'),'--spikein',str(p/'spikein.tsv'),'--out',str(p/'model_mixed_calibration')],capture_output=True,text=True)
 assert blocked.returncode!=0 and 'one shared calibration_group' in blocked.stderr,blocked.stderr
 spike_rows=list(csv.DictReader((p/'spikein.tsv').open(),delimiter='\t'))
 spike_rows[-1]['calibration_group']='prepB'
 with (p/'spikein.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(spike_rows[0]),delimiter='\t');w.writeheader();w.writerows(spike_rows)
 blocked=subprocess.run([sys.executable,str(ROOT/'scripts/run_differential.py'),'--config',str(p/'config/project.json'),'--counts',str(p/'counts.tsv'),'--metadata',str(p/'meta.tsv'),'--review',str(p/'review_spikein.json'),'--spikein',str(p/'spikein.tsv'),'--out',str(p/'model_mixed_calibration')],capture_output=True,text=True)
 assert blocked.returncode!=0 and 'one shared calibration_group' in blocked.stderr,blocked.stderr
 with (p/'model_spikein/size_factors.tsv').open() as f:scales=list(csv.DictReader(f,delimiter='\t'))
 gm=math.exp(sum(math.log(x) for x in factors)/len(factors))
 assert all(abs(float(row['size_factor'])-factor/gm)<1e-10 for row,factor in zip(scales,factors))
 with (p/'counts.tsv').open() as f:raw=next(csv.DictReader(f,delimiter='\t'))
 with (p/'model_spikein/normalized_counts.tsv').open() as f:normalized=next(csv.DictReader(f,delimiter='\t'))
 assert all(abs(float(normalized[name])-float(raw[name])/(factor/gm))<1e-7 for name,factor in zip(ids,factors))
 (p/'universe.txt').write_text(''.join(str(i)+'\n' for i in range(1,201)))
 (p/'genes.txt').write_text(''.join(str(i)+'\n' for i in range(1,21)))
 (p/'term2gene.tsv').write_text('term\tgene\n'+''.join('planted\t'+str(i)+'\n' for i in range(1,21))+''.join('background\t'+str(i)+'\n' for i in range(81,101)))
 subprocess.run(['Rscript',str(ROOT/'scripts/enrichment.R'),str(p/'genes.txt'),str(p/'universe.txt'),str(p/'term2gene.tsv'),str(p/'ora.tsv')],check=True)
 with (p/'ora.tsv').open() as f:terms=list(csv.DictReader(f,delimiter='\t'))
 assert any(r['ID']=='planted' and float(r['p.adjust'])<.001 for r in terms)

 print('PASS: conventional + spike-in DESeq2, known factors and normalization algebra, offline ORA')

 # Flat profiles / no finite adjusted tests must yield explicit unavailable panels.
 (p/'flat_model').mkdir()
 rcode='suppressPackageStartupMessages(library(DESeq2)); a<-commandArgs(TRUE); d<-readRDS(a[1]); counts(d)[]<-1L; sizeFactors(d)<-rep(1,ncol(d)); saveRDS(d,a[2]); z<-read.delim(a[3]); z$padj<-NA_real_; write.table(z,a[4],sep="\\t",quote=FALSE,row.names=FALSE)'
 subprocess.run(['Rscript','-e',rcode,str(p/'model/model.rds'),str(p/'flat_model/model.rds'),str(p/'model/complete.tsv'),str(p/'flat_model/complete.tsv')],check=True)
 subprocess.run(['Rscript',str(ROOT/'scripts/visualize.R'),str(p/'flat_model'),str(p/'flat_figures')],check=True)
 assert 'PCA unavailable' in (p/'flat_figures/unavailable_diagnostics.txt').read_text()
 assert 'Volcano unavailable' in (p/'flat_figures/unavailable_diagnostics.txt').read_text()
 print('PASS: flat profiles and all-NA adjusted tests produce explicit diagnostic gaps')

 (p/'acceptance.json').write_text(json.dumps({'state':'PASS','tested':['conventional DESeq2 direction','spike-in known size factors','raw counts divided by expected size factors','visualization including flat/NA profiles','offline ORA planted term'],'scope':'synthetic only; no experimental calibration approval'},indent=2))
