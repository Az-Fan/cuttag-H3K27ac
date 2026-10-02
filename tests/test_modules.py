import unittest,importlib.util,tempfile,subprocess,json,csv,gzip
from pathlib import Path
ROOT=Path(__file__).parents[1]
s=importlib.util.spec_from_file_location('consensus',ROOT/'scripts/consensus.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Modules(unittest.TestCase):
 def test_same_sample_overlap_not_double_counted(self):
  self.assertEqual(m.supported([[('chr1',0,20),('chr1',5,30)],[('chr1',10,25)]],2),[('chr1',10,25)])
 def test_touching_intervals_are_not_overlap(self):
  self.assertEqual(m.supported([[('chr1',0,10)],[('chr1',10,20)]],2),[])
 def test_consensus_blacklist_and_duplicate_units(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'a.bed').write_text('chr1\t0\t30\n');(p/'b.bed').write_text('chr1\t10\t40\n');(p/'black.bed').write_text('chr1\t15\t16\n')
   (p/'manifest.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\ng\ta\ta.bed\ng\tb\tb.bed\n')
   subprocess.run(['python3',str(ROOT/'scripts/consensus.py'),'--manifest',str(p/'manifest.tsv'),'--out',str(p/'out'),'--blacklist',str(p/'black.bed')],check=True)
   self.assertEqual((p/'out/master.bed').read_text(),'')
   self.assertEqual(json.loads((p/'out/provenance.json').read_text())['groups']['g']['independent_sample_count'],2)
 def test_release_missing_required_fails(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'release.json').write_text(json.dumps({'project':'x','artifacts':[{'path':'absent'}]}))
   r=subprocess.run(['python3',str(ROOT/'scripts/finalize.py'),'--manifest',str(p/'release.json'),'--out',str(p/'out')],capture_output=True)
   self.assertNotEqual(r.returncode,0);self.assertFalse((p/'out').exists())

class Audit(unittest.TestCase):
 def test_spikein_pair_counts_and_no_auto_acceptance(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'align.log').write_text('100 reads; of these:\n  100 (100.00%) were paired; of these:\n    80 (80.00%) aligned concordantly 0 times\n    15 (15.00%) aligned concordantly exactly 1 time\n    5 (5.00%) aligned concordantly >1 times\n')
   (p/'logs.tsv').write_text('sample_id\tmode\tbowtie2_log\ns1\toverlap\talign.log\n')
   subprocess.run(['python3',str(ROOT/'scripts/spikein_audit.py'),'--manifest',str(p/'logs.tsv'),'--out',str(p/'out')],check=True)
   with (p/'out/counts.tsv').open() as f:z=list(csv.DictReader(f,delimiter='\t'))[0]
   self.assertEqual(z['spikein_fragments'],'20');self.assertEqual(z['calibration_accepted'],'False')
 def test_release_hash(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'x.tsv').write_text('x');(p/'release.json').write_text(json.dumps({'project':'x','artifacts':[{'path':'x.tsv'}]}))
   subprocess.run(['python3',str(ROOT/'scripts/finalize.py'),'--manifest',str(p/'release.json'),'--out',str(p/'out')],check=True)
   z=json.loads((p/'out/manifest.json').read_text());self.assertEqual(len(z['artifacts'][0]['sha256']),64)
