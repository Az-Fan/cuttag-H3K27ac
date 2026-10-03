"""Exercise known-motif execution and the small-input skip gate using local synthetic DNA."""
import argparse
import json
import random
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
rng=random.Random(872)
sequence=''.join(rng.choice('ACGT') for _ in range(240000))
(a.out/'genome.fa').write_text('>chr1\n'+sequence+'\n')
subprocess.run(['samtools','faidx',str(a.out/'genome.fa')],check=True)
(a.out/'target.bed').write_text(''.join('chr1\t%d\t%d\n'%(i*1000+200,i*1000+400) for i in range(60)))
(a.out/'background.bed').write_text(''.join('chr1\t%d\t%d\n'%(i*1000+200,i*1000+400) for i in range(60,240)))
command=[sys.executable,str(ROOT/'scripts/motif.py'),'--target',str(a.out/'target.bed'),'--background',str(a.out/'background.bed'),'--fasta',str(a.out/'genome.fa'),'--threads','1','--execute']
subprocess.run(command+['--out',str(a.out/'known')],check=True)
assert json.loads((a.out/'known/status.json').read_text())['state']=='COMPUTATIONAL_PASS'
subprocess.run(command+['--out',str(a.out/'small'),'--min-windows','100'],check=True)
assert json.loads((a.out/'small/status.json').read_text())['state']=='SKIPPED_TOO_FEW_WINDOWS'
(a.out/'acceptance.json').write_text(json.dumps({'state':'PASS','tested':['known motif runtime','insufficient windows skip'],'scope':'random DNA; no asserted motif enrichment'},indent=2))
print('PASS: known HOMER motifs and small-window skip gate')
