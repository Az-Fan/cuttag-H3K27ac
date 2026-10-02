import importlib.util,unittest,json,copy
from pathlib import Path
s=importlib.util.spec_from_file_location('cuttag',Path(__file__).parents[1]/'scripts/cuttag.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Contracts(unittest.TestCase):
 def setUp(self): self.root,self.cfg,self.rows=m.load(m.ROOT/'config/project.json')
 def test_one_per_group_blocks_inference(self):
  self.cfg['analysis']['replicates_confirmed']=True
  self.assertFalse(m.validate(self.root,self.cfg,self.rows,True)['formal_inference_eligible'])
 def test_technical_split_does_not_add_replicates(self):
  r=copy.deepcopy(self.rows[0]);r.update(unit_id='lane2',fastq_1='data/raw/split_R1.fastq.gz',fastq_2='data/raw/split_R2.fastq.gz');self.rows.append(r)
  v=m.validate(self.root,self.cfg,self.rows,True)
  self.assertEqual(v['biological_samples_per_condition']['Control'],1);self.assertEqual(v['errors'],[])
 def test_relabeling_same_bio_as_new_replicate_fails(self):
  r=copy.deepcopy(self.rows[0]);r.update(unit_id='lane2',replicate='2',fastq_1='data/raw/split_R1.fastq.gz',fastq_2='data/raw/split_R2.fastq.gz');self.rows.append(r)
  self.assertTrue(m.validate(self.root,self.cfg,self.rows,True)['errors'])
 def test_reusing_fastq_fails(self):
  self.rows[1]['fastq_1']=self.rows[0]['fastq_1'];self.assertTrue(m.validate(self.root,self.cfg,self.rows,True)['errors'])
if __name__=='__main__':unittest.main()
