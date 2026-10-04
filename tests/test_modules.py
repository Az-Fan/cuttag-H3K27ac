import unittest,importlib.util,tempfile,subprocess,json,csv,gzip
from pathlib import Path
ROOT=Path(__file__).parents[1]
s=importlib.util.spec_from_file_location('consensus',ROOT/'scripts/consensus.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Modules(unittest.TestCase):
 def test_consensus_rejects_peak_outside_target_reference(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'sizes').write_text('chr1\t1000\n')
   (p/'manifest.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\n')
   for peak in ('chr2\t0\t100\n','chr1\t900\t1100\n'):
    (p/'a.bed').write_text(peak)
    result=subprocess.run(['python3',str(ROOT/'scripts/consensus.py'),'--manifest',str(p/'manifest.tsv'),'--sizes',str(p/'sizes'),'--out',str(p/'out')],capture_output=True,text=True)
    self.assertNotEqual(result.returncode,0);self.assertIn('outside target reference',result.stderr);self.assertFalse((p/'out').exists())
 def test_default_min_width_discards_blacklist_fragments_below_50bp(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'a.bed').write_text('chr1\t100\t180\n');(p/'black.bed').write_text('chr1\t145\t148\n')
   (p/'manifest.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\n')
   subprocess.run(['python3',str(ROOT/'scripts/consensus.py'),'--manifest',str(p/'manifest.tsv'),'--out',str(p/'out'),'--blacklist',str(p/'black.bed')],check=True)
   self.assertEqual((p/'out/master.bed').read_text(),'')
   decisions=json.loads((p/'out/interval_decisions.json').read_text())[0]
   self.assertEqual(decisions['short_pieces_discarded'],[['chr1',100,145],['chr1',148,180]])
   self.assertEqual(json.loads((p/'out/provenance.json').read_text())['min_width'],50)
 def test_same_sample_overlap_not_double_counted(self):
  self.assertEqual(m.supported([[('chr1',0,20),('chr1',5,30)],[('chr1',10,25)]],2),[('chr1',10,25)])
 def test_touching_intervals_are_not_overlap(self):
  self.assertEqual(m.supported([[('chr1',0,10)],[('chr1',10,20)]],2),[])
 def test_consensus_blacklist_and_duplicate_units(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'a.bed').write_text('chr1\t0\t30\n');(p/'b.bed').write_text('chr1\t10\t40\n');(p/'black.bed').write_text('chr1\t15\t16\n')
   (p/'manifest.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\ng\ta\ta.bed\ng\tb\tb.bed\n')
   subprocess.run(['python3',str(ROOT/'scripts/consensus.py'),'--manifest',str(p/'manifest.tsv'),'--out',str(p/'out'),'--blacklist',str(p/'black.bed'),'--blacklist-mode','drop'],check=True)
   self.assertEqual((p/'out/master.bed').read_text(),'')
   self.assertEqual(json.loads((p/'out/provenance.json').read_text())['groups']['g']['independent_sample_count'],2)
 def test_consensus_blacklist_subtracts_only_overlapping_bases(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'a.bed').write_text('chr1\t100\t180\n');(p/'black.bed').write_text('chr1\t145\t148\n')
   (p/'manifest.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\n')
   subprocess.run(['python3',str(ROOT/'scripts/consensus.py'),'--manifest',str(p/'manifest.tsv'),'--out',str(p/'out'),'--blacklist',str(p/'black.bed'),'--min-width','1'],check=True)
   self.assertEqual((p/'out/g.consensus.bed').read_text(),'chr1\t100\t145\nchr1\t148\t180\n')
 def test_reproducible_union_preserves_full_supported_sample_intervals(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'a.bed').write_text('chr1\t100\t180\n');(p/'b.bed').write_text('chr1\t120\t200\n')
   (p/'manifest.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\ng\tb\tb.bed\n')
   subprocess.run(['python3',str(ROOT/'scripts/consensus.py'),'--manifest',str(p/'manifest.tsv'),'--out',str(p/'out'),'--universe','reproducible_union','--min-width','1'],check=True)
   self.assertEqual((p/'out/master.bed').read_text(),'chr1\t100\t200\tpeak_00000001\n')
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
   self.assertEqual(z['concordant_unique_pairs'],'15');self.assertEqual(z['concordant_multi_pairs'],'5');self.assertEqual(z['concordant_pairs_diagnostic_only'],'20');self.assertEqual(z['calibration_accepted'],'False')
 def test_release_hash(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'x.tsv').write_text('x');(p/'release.json').write_text(json.dumps({'project':'x','artifacts':[{'path':'x.tsv'}]}))
   subprocess.run(['python3',str(ROOT/'scripts/finalize.py'),'--manifest',str(p/'release.json'),'--out',str(p/'out')],check=True)
   z=json.loads((p/'out/manifest.json').read_text());self.assertEqual(len(z['artifacts'][0]['sha256']),64)
 def test_workflow_blocks_downstream_after_failure(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'w.json').write_text(json.dumps({'steps':[{'id':'fail','command':['python3','-c','raise SystemExit(3)']},{'id':'downstream','depends_on':['fail'],'command':['python3','-c','print(1)']}]}))
   r=subprocess.run(['python3',str(ROOT/'scripts/workflow.py'),'--workflow',str(p/'w.json'),'--out',str(p/'out')])
   self.assertEqual(r.returncode,1);z=json.loads((p/'out/workflow_status.json').read_text());self.assertEqual(z['steps'][1]['state'],'BLOCKED')
