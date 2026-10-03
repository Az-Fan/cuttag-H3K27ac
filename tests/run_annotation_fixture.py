"""Check peak ID/coordinate round-tripping using an explicit locked human TxDb."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=False)
# This fixture uses hg19 deliberately and records it; production accepts an explicit assembly-specific SQLite.
script='suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene)); a<-commandArgs(TRUE); AnnotationDbi::saveDb(TxDb.Hsapiens.UCSC.hg19.knownGene,a[1])'
subprocess.run(['Rscript','-e',script,str(a.out/'hg19.sqlite')],check=True)
(a.out/'peaks.bed').write_text('chr6\t122965000\t122965200\tz_peak\nchr1\t11800\t12000\ta_peak\nchr6\t122970000\t122970200\tm_peak\n')
subprocess.run(['Rscript',str(ROOT/'scripts/annotate.R'),str(a.out/'peaks.bed'),str(a.out/'hg19.sqlite'),'org.Hs.eg.db',str(a.out/'annotation.tsv')],check=True)
with (a.out/'annotation.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert [(r['peak_id'],int(r['start']),int(r['end'])) for r in rows]==[('z_peak',122965000,122965200),('a_peak',11800,12000),('m_peak',122970000,122970200)]
(a.out/'acceptance.json').write_text(json.dumps({'state':'PASS','scope':'arbitrary test intervals using explicit hg19 TxDb and locked OrgDb; not a biological validation','tested':['peak ID order','BED to GRanges and back','explicit TxDb and OrgDb loading']},indent=2))
print('PASS: annotation IDs/order/coordinates with explicit hg19 TxDb and locked OrgDb')
