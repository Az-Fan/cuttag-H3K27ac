"""Exercise real DESeq2 on small deterministic overdispersed synthetic counts."""
import csv,random,subprocess,tempfile
from pathlib import Path
ROOT=Path(__file__).parents[1];rng=random.Random(513)
with tempfile.TemporaryDirectory() as d:
 p=Path(d);ids=['C1','C2','C3','K1','K2','K3']
 with (p/'counts.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['peak_id','chrom','start','end']+ids)
  for i in range(300):
   mu=rng.uniform(20,500);fold=3 if i<60 else 1
   v=[max(1,round(rng.gammavariate(8,mu*(fold if j>=3 else 1)/8))) for j in range(6)]
   w.writerow([f'p{i}','chr1',i*1000,i*1000+200]+v)
 with (p/'meta.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','biological_sample_id','condition','replicates_confirmed'])
  for name in ids:w.writerow([name,name,'Control' if name[0]=='C' else 'KD','TRUE'])
 subprocess.run(['Rscript',str(ROOT/'scripts/differential.R'),str(p/'counts.tsv'),str(p/'meta.tsv'),'KD','Control',str(p/'model'),'conventional'],check=True)
 subprocess.run(['Rscript',str(ROOT/'scripts/visualize.R'),str(p/'model'),str(p/'figures')],check=True)
 z=list(csv.DictReader((p/'model/complete.tsv').open(),delimiter='\t'))
 assert len(z)==300 and sum(float(r['log2FoldChange'])>0 for r in z if int(r['peak_id'][1:])<60)>50
 assert (p/'figures/diagnostics.pdf').stat().st_size>0
 print('Synthetic DESeq2 and visualization PASS; contrast direction verified')
