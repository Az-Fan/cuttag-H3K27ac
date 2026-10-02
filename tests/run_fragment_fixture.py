"""Real samtools/bedtools/pyBigWig integration on tiny synthetic proper pairs."""
import csv,json,os,subprocess,sys,tempfile
from pathlib import Path
import pyBigWig
ROOT=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory() as d:
 p=Path(d);sam=['@HD\tVN:1.6\tSO:unsorted','@SQ\tSN:chr1\tLN:1000']
 for name,pos,dup,low in [('regular',11,False,False),('duplicate',101,True,False),('lowmate',201,False,True)]:
  for flag,offset,mate,tlen,mapq in [(99,0,20,30,60),(147,20,0,-30,0 if low else 60)]:
   sam.append('\t'.join(map(str,[name,flag+(1024 if dup else 0),'chr1',pos+offset,mapq,'10M','=',pos+mate,tlen,'ACGTACGTAC','IIIIIIIIII'])))
 (p/'x.sam').write_text('\n'.join(sam)+'\n');subprocess.run(['samtools','view','-b','-o',str(p/'x.bam'),str(p/'x.sam')],check=True)
 subprocess.run([sys.executable,str(ROOT/'scripts/fragments.py'),'--bam',str(p/'x.bam'),'--out',str(p/'frags')],check=True)
 assert (p/'frags/fragments.bed').read_text().splitlines()==['chr1\t100\t130\tduplicate','chr1\t10\t40\tregular']
 subprocess.run([sys.executable,str(ROOT/'scripts/fragments.py'),'--bam',str(p/'x.bam'),'--out',str(p/'dedup'),'--remove-duplicates'],check=True)
 assert json.loads((p/'dedup/summary.json').read_text())['fragments']==1
 (p/'peaks.bed').write_text('chr1\t0\t50\tp1\nchr1\t100\t150\tp2\n');(p/'fragments.tsv').write_text('sample_id\tfragments_bed\ns1\tfrags/fragments.bed\n')
 subprocess.run([sys.executable,str(ROOT/'scripts/count_fragments.py'),'--peaks',str(p/'peaks.bed'),'--samples',str(p/'fragments.tsv'),'--out',str(p/'counts.tsv')],check=True,cwd='/tmp')
 with (p/'counts.tsv').open() as f:z=list(csv.DictReader(f,delimiter='\t'))
 assert [r['s1'] for r in z]==['1','1']
 (p/'sizes').write_text('chr1\t1000\n')
 subprocess.run([sys.executable,str(ROOT/'scripts/tracks.py'),'--manifest',str(p/'fragments.tsv'),'--sizes',str(p/'sizes'),'--genome','synthetic','--out',str(p/'tracks')],check=True)
 with pyBigWig.open(str(p/'tracks/s1.CPM.bw')) as bw:assert bw.values('chr1',10,11)==[500000.0]
 assert 'genome="synthetic"' in (p/'tracks/s1.CPM.igv.xml').read_text()
 (p/'badtracks.tsv').write_text('sample_id\tfragments_bed\tspikein_scale\ns1\tfrags/fragments.bed\tnan\n')
 bad=subprocess.run([sys.executable,str(ROOT/'scripts/tracks.py'),'--manifest',str(p/'badtracks.tsv'),'--sizes',str(p/'sizes'),'--genome','synthetic','--out',str(p/'badtracks')],capture_output=True)
 assert bad.returncode!=0 and not (p/'badtracks').exists()
 (p/'metrics.tsv').write_text('sample_id\tfrip\tspikein_fraction\ns1\tnan\t0.25\n')
 env=os.environ.copy();env['MPLCONFIGDIR']=str(p/'mpl')
 subprocess.run([sys.executable,str(ROOT/'scripts/qc_atlas.py'),'--metrics',str(p/'metrics.tsv'),'--out',str(p/'qc')],env=env,check=True)
 assert json.loads((p/'qc/missing_metrics.json').read_text())['frip']==['s1']
 print('Real fragment/count/CPM BigWig/QC fixture PASS; both-mate MAPQ and duplicate policy verified')
