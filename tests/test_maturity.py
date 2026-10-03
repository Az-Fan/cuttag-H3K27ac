import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from collect_qc import collect
from compare_peaks import overlap_bp
from finalize import package
from init_project import initialize
from prepare_universe import prepare
from reference_audit import audit
from verify_release import verify
from workflow import validate

class Maturity(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name)

    def test_init_is_clean_and_does_not_copy_approvals_or_cache(self):
        destination=self.p/'project'
        initialize(ROOT,destination,'study1')
        cfg=json.loads((destination/'config/project.json').read_text())
        self.assertEqual(cfg['project_id'],'study1')
        self.assertIs(cfg['analysis']['replicates_confirmed'],False)
        self.assertEqual(list((destination/'data/raw').iterdir()),[])
        self.assertFalse((destination/'.git').exists())
        self.assertFalse((destination/'.pixi').exists())
        self.assertTrue((destination/'.github/workflows/contracts.yml').is_file())
        self.assertTrue((destination/'docs/validation/2026-10-03/official_test.json').is_file())
        r=subprocess.run([sys.executable,str(destination/'scripts/cuttag.py'),'validate','--allow-missing'],capture_output=True)
        self.assertEqual(r.returncode,0,r.stderr)
        with self.assertRaises(ValueError):initialize(ROOT,destination,'study2')

    def test_qc_counts_each_fragment_once_with_half_open_bounds(self):
        (self.p/'fragments.bed').write_text('chr1\t0\t10\na\t0\t40\nchr1\t10\t30\nchr1\t10\t30\n')
        (self.p/'peaks.bed').write_text('chr1\t10\t20\nchr1\t15\t25\n')
        (self.p/'m.tsv').write_text('sample_id\tfragments_bed\tpeaks_bed\ns1\tfragments.bed\tpeaks.bed\n')
        collect(self.p/'m.tsv',self.p/'qc')
        with (self.p/'qc/metrics.tsv').open() as f:row=next(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(row['target_fragments'],'4')
        self.assertEqual(row['frip'],'0.5')
        self.assertEqual(row['median_fragment_length'],'20.0')
        self.assertEqual(row['peak_union_bp'],'15')
        self.assertEqual(row['duplication_rate'],'NA')

    def test_qc_empty_fragments_fail_without_report(self):
        (self.p/'f').write_text('');(self.p/'b').write_text('')
        (self.p/'m').write_text('sample_id\tfragments_bed\tpeaks_bed\ns1\tf\tb\n')
        with self.assertRaises(ValueError):collect(self.p/'m',self.p/'out')
        self.assertFalse((self.p/'out').exists())

    def reference(self):
        (self.p/'genome.fa').write_text('>chr1\n'+'A'*100+'\n')
        (self.p/'genes.gtf').write_text('chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "g";\n')
        (self.p/'black.bed').write_text('chr1\t99\t100\n')
        return {'reference':{'genome':'fixture','fasta':'genome.fa','gtf':'genes.gtf','blacklist':'black.bed','mito_name':''},'spikein':{'enabled':False}}

    def test_reference_bounds_and_hashes(self):
        cfg=self.reference()
        result=audit(self.p,cfg)
        self.assertEqual(result['target_sizes'],{'chr1':100})
        self.assertEqual(len(result['files']['gtf']['sha256']),64)
        (self.p/'black.bed').write_text('chr1\t99\t101\n')
        with self.assertRaises(ValueError):audit(self.p,cfg)

    def test_reference_wrong_assembly_dictionary_fails(self):
        cfg=self.reference();(self.p/'genes.gtf').write_text('1\ttest\tgene\t1\t90\t.\t+\t.\tx\n')
        with self.assertRaises(ValueError):audit(self.p,cfg)

    def test_reference_stale_fai_fails(self):
        cfg=self.reference();(self.p/'genome.fa.fai').write_text('chr1\t99\t6\t100\t101\n')
        with self.assertRaises(ValueError):audit(self.p,cfg)

    def test_workflow_rejects_unknown_dependency_before_execution(self):
        with self.assertRaises(ValueError):validate({'steps':[{'id':'x','command':['true'],'depends_on':['missing']}]})

    def workflow(self):
        (self.p/'input').write_text('original')
        spec={'working_directory':str(self.p),'steps':[{'id':'copy','command':[sys.executable,'-c','from pathlib import Path; Path("output").write_text(Path("input").read_text())'],'inputs':['input'],'outputs':['output']}]}
        (self.p/'w.json').write_text(json.dumps(spec))
        return [sys.executable,str(ROOT/'scripts/workflow.py'),'--workflow',str(self.p/'w.json'),'--out',str(self.p/'run')]

    def test_resume_checks_content_and_skips_valid_outputs(self):
        command=self.workflow()
        subprocess.run(command,check=True,capture_output=True)
        stamp=(self.p/'output').stat().st_mtime_ns
        subprocess.run(command+['--resume'],check=True,capture_output=True)
        self.assertEqual(stamp,(self.p/'output').stat().st_mtime_ns)
        (self.p/'input').write_text('replaced')
        self.assertNotEqual(subprocess.run(command+['--resume'],capture_output=True).returncode,0)
        self.assertEqual((self.p/'output').read_text(),'original')

    def test_resume_detects_corrupt_output(self):
        command=self.workflow();subprocess.run(command,check=True,capture_output=True)
        (self.p/'output').write_text('corrupt')
        self.assertNotEqual(subprocess.run(command+['--resume'],capture_output=True).returncode,0)

    def test_release_survives_move_and_detects_tampering(self):
        assets=self.p/'assets';assets.mkdir()
        (assets/'s.bw').write_bytes(b'fixture')
        (assets/'s.igv.xml').write_text('<Session genome="hg38"><Resources><Resource path="s.bw"/></Resources></Session>')
        spec=self.p/'release.json';spec.write_text(json.dumps({'project':'x','artifacts':[{'path':'assets','destination':'tracks'}]}))
        package(spec,self.p/'delivery')
        shutil.rmtree(assets)
        shutil.move(str(self.p/'delivery'),str(self.p/'moved'))
        self.assertEqual(verify(self.p/'moved')['files_verified'],2)
        self.assertIn('href="tracks/s.bw"',(self.p/'moved/index.html').read_text())
        (self.p/'moved/tracks/s.bw').write_bytes(b'corrupt')
        with self.assertRaises(ValueError):verify(self.p/'moved')

    def test_missing_igv_resource_leaves_no_package(self):
        (self.p/'s.igv.xml').write_text('<Session><Resources><Resource path="missing.bw"/></Resources></Session>')
        spec=self.p/'r.json';spec.write_text(json.dumps({'project':'x','artifacts':[{'path':'s.igv.xml'}]}))
        with self.assertRaises(ValueError):package(spec,self.p/'out')
        self.assertFalse((self.p/'out').exists())

    def test_package_rejects_parent_traversal(self):
        (self.p/'x').write_text('x')
        spec=self.p/'r.json';spec.write_text(json.dumps({'project':'x','artifacts':[{'path':'x','destination':'../x'}]}))
        with self.assertRaises(ValueError):package(spec,self.p/'out')

    def test_universe_excludes_untested_and_preserves_both_directions(self):
        (self.p/'results').write_text('peak_id\tchrom\tstart\tend\tpadj\tlog2FoldChange\np1\tchr1\t0\t10\t0.01\t2\np2\tchr1\t20\t30\t0.01\t-2\np3\tchr1\t40\t50\tNA\t3\np4\tchr1\t60\t70\t0.8\t0\n')
        (self.p/'ann').write_text('peak_id\tseqnames\tstart\tend\tgeneId\tannotation\np1\tchr1\t0\t10\t11\tPromoter (<=1kb)\np2\tchr1\t20\t30\t11\tDistal Intergenic\np3\tchr1\t40\t50\t12\tPromoter\np4\tchr1\t60\t70\t13\tExon\n')
        prepare(self.p/'results',self.p/'ann',self.p/'universe',.05,1)
        self.assertEqual((self.p/'universe/all/universe.txt').read_text(),'11\n13\n')
        z=json.loads((self.p/'universe/provenance.json').read_text())
        self.assertEqual(z['regions']['all']['genes_with_both_gain_and_loss_peaks'],['11'])
        self.assertNotIn('p1',(self.p/'universe/all/gain.background.bed').read_text())

    def test_peak_intersection_respects_chromosomes_and_touching(self):
        self.assertEqual(overlap_bp([('chr1',0,10),('chr2',0,10)],[('chr1',10,20),('chr2',5,15)]),5)

    def test_consensus_fraction_grid(self):
        for name,end in [('a',30),('b',20),('c',10)]:
            (self.p/(name+'.bed')).write_text('chr1\t0\t%d\n'%end)
        (self.p/'peaks.tsv').write_text('group\tbiological_sample_id\tpeaks_bed\ng\ta\ta.bed\ng\tb\tb.bed\ng\tc\tc.bed\n')
        for fraction,end in [('0.3333333333333333',30),('0.6666666666666666',20),('1',10)]:
            dest=self.p/('fraction_'+fraction)
            subprocess.run([sys.executable,str(ROOT/'scripts/consensus.py'),'--manifest',str(self.p/'peaks.tsv'),'--fraction',fraction,'--out',str(dest)],check=True,capture_output=True)
            self.assertEqual((dest/'master.bed').read_text().split()[2],str(end))

    def test_spikein_calibration_groups_not_pooled(self):
        lines=['sample_id\tmode\tcalibration_group\tbowtie2_log']
        for sample,group,count in [('a','batch1',10),('b','batch1',40),('c','batch2',100)]:
            (self.p/(sample+'.log')).write_text('%d reads; of these:\n  0 (0%%) aligned concordantly 0 times\n  %d (100%%) aligned concordantly exactly 1 time\n  0 (0%%) aligned concordantly >1 times\n'%(count,count))
            lines.append('\t'.join([sample,'overlap',group,sample+'.log']))
        (self.p/'m.tsv').write_text('\n'.join(lines)+'\n')
        subprocess.run([sys.executable,str(ROOT/'scripts/spikein_audit.py'),'--manifest',str(self.p/'m.tsv'),'--out',str(self.p/'audit')],check=True,capture_output=True)
        with (self.p/'audit/counts.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertAlmostEqual(float(rows[0]['size_factor']),.5)
        self.assertAlmostEqual(float(rows[1]['size_factor']),2)
        self.assertAlmostEqual(float(rows[2]['size_factor']),1)

    def test_nextflow_ignored_failure_remains_visible(self):
        from audit_run import audit
        (self.p/'trace.tsv').write_text('task_id\tname\tstatus\texit\n1\talign\tCOMPLETED\t0\n2\tpreseq\tFAILED\t1\n')
        result=audit(self.p)
        self.assertEqual(result['state'],'REVIEW_REQUIRED')
        self.assertEqual(result['nonpassing_tasks'][0]['name'],'preseq')

    def test_nfcore_adapter_preserves_sample_mapping_and_rejects_changed_snapshot(self):
        from import_nfcore import convert
        from cuttag import load
        root=self.p/'project';initialize(ROOT,root,'adapter_fixture')
        _,cfg,samples=load(root/'config/project.json')
        run=root/'results/upstream';run.mkdir(parents=True)
        (run/'status.json').write_text(json.dumps({'state':'COMPUTATIONAL_PASS'}))
        (run/'project.snapshot.json').write_text(json.dumps(cfg))
        (run/'samples.snapshot.json').write_text(json.dumps(samples))
        (run/'params.json').write_text(json.dumps(dict(cfg['nfcore_params'],only_filtering=False)))
        for sample in samples:
            if sample['role']!='target':continue
            prefix=sample['group']+'_R'+sample['replicate']
            bam=run/'output/02_alignment/bowtie2/target/markdup'/(prefix+'.target.markdup.sorted.bam')
            peak=run/'output/03_peak_calling/04_called_peaks/macs2'/(prefix+'.macs2_peaks.narrowPeak')
            bam.parent.mkdir(parents=True,exist_ok=True);peak.parent.mkdir(parents=True,exist_ok=True)
            bam.write_text('layout fixture '+sample['sample_id']);peak.write_text('chr1\t0\t10\n')
        convert(root/'config/project.json',run,root/'results/import')
        with (root/'results/import/artifacts.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual({r['sample_id'] for r in rows},{'Control1','KD1'})
        samples[0]['biological_sample_id']='different'
        (run/'samples.snapshot.json').write_text(json.dumps(samples))
        with self.assertRaises(ValueError):convert(root/'config/project.json',run,root/'results/import_changed')
