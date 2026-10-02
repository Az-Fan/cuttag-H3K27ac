import copy,hashlib,importlib.util,json,subprocess,sys,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import cuttag
from decisions import require_review,sha256
class AuditRegressions(unittest.TestCase):
 def setUp(self):self.root,self.cfg,self.rows=cuttag.load(ROOT/'config/project.json')
 def test_string_false_rejected(self):
  self.cfg['analysis']['replicates_confirmed']='false'
  self.assertTrue(cuttag.validate(self.root,self.cfg,self.rows,True)['errors'])
 def test_missing_group_rejected(self):
  self.cfg['analysis']['numerator']='absent'
  self.assertTrue(cuttag.validate(self.root,self.cfg,self.rows,True)['errors'])
 def test_same_contrast_rejected(self):
  self.cfg['analysis']['numerator']='Control'
  self.assertTrue(cuttag.validate(self.root,self.cfg,self.rows,True)['errors'])
 def test_batch_design_not_silently_ignored(self):
  self.cfg['analysis']['design']='~ batch + condition'
  self.assertTrue(cuttag.validate(self.root,self.cfg,self.rows,True)['errors'])
 def test_malformed_table_reports_error(self):
  self.assertTrue(cuttag.validate(self.root,self.cfg,[{'sample_id':'x'}],True)['errors'])
 def test_unknown_spikein_stage_blocks_production(self):
  self.cfg['qc']['upstream_accepted']=True
  self.cfg['spikein'].update(enabled=True,calibration_accepted=True,equal_amount_confirmed=True,added_at='unknown')
  with self.assertRaisesRegex(ValueError,'addition stage'):cuttag.production_gate(ROOT,self.cfg,ROOT/'config/project.json')
 def test_production_needs_control_review(self):
  self.cfg['qc']['upstream_accepted']=True
  with self.assertRaisesRegex(ValueError,'QC/control-strategy review'):cuttag.production_gate(ROOT,self.cfg,ROOT/'config/project.json')
 def test_spikein_aligned_before_calibration(self):
  self.cfg['spikein'].update(enabled=True,fasta='lambda.fa',alignment_profile='config/profiles/spikein_overlap_diagnostic.config')
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'plan';cmd=cuttag.plan(ROOT,self.cfg,self.rows,'alignment',out)
   params=json.loads((out/'params.json').read_text())
   self.assertEqual(params['normalisation_mode'],'Spikein');self.assertTrue(params['only_filtering']);self.assertIsNone(params['spikein_bowtie2']);self.assertIn(str(ROOT/self.cfg['spikein']['alignment_profile']),cmd)
 def test_review_bound_to_input_and_evidence(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'input').write_text('counts');(p/'evidence').write_text('QC decision')
   record={'reviewer':'fixture','reviewed_at':'2026-10-03','reason':'synthetic review','decisions':{'peaks_accepted':True},'inputs':{'counts_sha256':sha256(p/'input')},'evidence':[{'path':'evidence','sha256':sha256(p/'evidence')}]}
   (p/'review.json').write_text(json.dumps(record))
   require_review(p/'review.json',['peaks_accepted'],{'counts_sha256':p/'input'})
   (p/'input').write_text('changed counts')
   with self.assertRaisesRegex(ValueError,'input missing or changed'):require_review(p/'review.json',['peaks_accepted'],{'counts_sha256':p/'input'})
 def test_false_peak_review_blocks_even_with_other_acceptances(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'review.json';p.write_text(json.dumps({'reviewer':'fixture','reviewed_at':'2026','reason':'QC','decisions':{'upstream_accepted':True,'peaks_accepted':False}}))
   with self.assertRaisesRegex(ValueError,'peaks_accepted'):require_review(p,['upstream_accepted','peaks_accepted'],{})
 def test_zero_exit_missing_artifact_is_failure(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'workflow.json').write_text(json.dumps({'working_directory':d,'steps':[{'id':'step','command':[sys.executable,'-c','pass'],'outputs':['missing.tsv']}]}))
   r=subprocess.run([sys.executable,str(ROOT/'scripts/workflow.py'),'--workflow',str(p/'workflow.json'),'--out',str(p/'out')],capture_output=True)
   self.assertEqual(r.returncode,1)
   self.assertEqual(json.loads((p/'out/workflow_status.json').read_text())['steps'][0]['missing_outputs'],['missing.tsv'])
 def test_configured_differential_checks_peak_gate_before_launch(self):
  cfg=copy.deepcopy(self.cfg);cfg['analysis']['replicates_confirmed']=True;cfg['qc'].update(upstream_accepted=True,peaks_accepted=False)
  rows=copy.deepcopy(self.rows)
  for r in self.rows:
   if r['role']=='target':
    z={k:(v+'2' if k in ('sample_id','biological_sample_id','library_id','unit_id','fastq_1','fastq_2') else v) for k,v in r.items()};z['replicate']='2';z['fastq_1']='data/x'+z['sample_id']+'_R1.fastq.gz';z['fastq_2']='data/x'+z['sample_id']+'_R2.fastq.gz';rows.append(z)
  with tempfile.TemporaryDirectory() as d:
   import csv
   p=Path(d);(p/'config').mkdir();cfg['samples']='config/samples.tsv'
   (p/'config/project.json').write_text(json.dumps(cfg))
   with (p/'config/samples.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
   r=subprocess.run([sys.executable,str(ROOT/'scripts/run_differential.py'),'--config',str(p/'config/project.json'),'--counts','unused','--metadata','unused','--review','unused','--out',str(p/'out')],capture_output=True,text=True)
   self.assertNotEqual(r.returncode,0);self.assertIn('Upstream and peak QC must be accepted',r.stderr)
 def test_production_review_rejects_changed_inventory(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'config').mkdir();(p/'results/qc').mkdir(parents=True)
   cfg=copy.deepcopy(self.cfg);cfg['qc'].update(upstream_accepted=True,production_review='review.json')
   config=p/'config/project.json';config.write_text(json.dumps(cfg));(p/'config/samples.tsv').write_text('samples');inventory=p/'results/qc/fastq_inventory.json';inventory.write_text('["original FASTQ hashes"]');(p/'evidence.txt').write_text('Synthetic QC and control review')
   review={'reviewer':'fixture','reviewed_at':'2026-10-03','reason':'synthetic','decisions':{'upstream_accepted':True,'control_strategy_accepted':True},'inputs':{'config_sha256':sha256(config),'samples_sha256':sha256(p/'config/samples.tsv'),'fastq_inventory_sha256':sha256(inventory)},'evidence':[{'path':'evidence.txt','sha256':sha256(p/'evidence.txt')}]}
   (p/'review.json').write_text(json.dumps(review));cuttag.production_gate(p,cfg,config)
   inventory.write_text('["replaced FASTQ hashes"]')
   with self.assertRaisesRegex(ValueError,'fastq_inventory_sha256'):cuttag.production_gate(p,cfg,config)
 def test_pipeline_group_cannot_mix_conditions(self):
  self.rows[1]['group']=self.rows[0]['group'];self.rows[1]['replicate']='2'
  errors=cuttag.validate(self.root,self.cfg,self.rows,True)['errors']
  self.assertIn('Pipeline group mixes conditions/roles/control strategy',errors)
 def test_leading_zero_replicate_cannot_alias_pipeline_integer(self):
  self.rows[0]['replicate']='01'
  self.assertTrue(cuttag.validate(self.root,self.cfg,self.rows,True)['errors'])
