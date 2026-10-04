import csv
import json
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from compare_models import compare as compare_models
from compare_strategy_matrix import compare as compare_strategies, verify_run
from strategy_matrix import build as build_matrix, expand


class StrategyTools(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_matrix_expands_only_declared_candidate_combinations(self):
        spec = {'artifact_sets': [{'id': 'narrow_igg', 'manifest': 'narrow.tsv'},
                                  {'id': 'broad_igg', 'manifest': 'broad.tsv'}],
                'axes': {'mapq': [20, 30], 'remove_duplicates': [False, True],
                         'universe': ['support_core', 'reproducible_union']}}
        result = expand(spec)
        self.assertEqual(len(result), 16)
        self.assertEqual({r['min_width'] for r in result}, {50})
        self.assertEqual(len({r['candidate_id'] for r in result}), 16)
        self.assertEqual({r['artifact_set'] for r in result}, {'narrow_igg', 'broad_igg'})
        with self.assertRaises(ValueError):
            expand({**spec, 'axes': {'universe': ['unknown']}})

    def test_template_matrix_compares_minimum_width_sensitivity_and_declares_baseline(self):
        example = json.loads((ROOT / 'examples/strategy_matrix.json').read_text())
        candidates = expand(example)
        self.assertEqual({row['min_width'] for row in candidates}, {1, 50, 100})
        self.assertIn(example['baseline_candidate'], {row['candidate_id'] for row in candidates})
        self.assertEqual(len(candidates), 48)

    def test_matrix_registry_resumes_without_rebuilding_completed_candidates(self):
        config = self.root / 'config.json'; config.write_text('{}')
        artifact = self.root / 'artifacts.tsv'; artifact.write_text('header\nrow\n')
        matrix = self.root / 'matrix.json'
        matrix.write_text(json.dumps({'artifact_sets': [{'id': 'base', 'manifest': 'artifacts.tsv'}]}))

        def fake_plan(config_path, artifact_path, out, *args):
            out.mkdir(parents=True)
            (out / 'workflow.json').write_text('{"steps": []}')
            return out / 'workflow.json'

        with patch('strategy_matrix.plan', side_effect=fake_plan) as plan_mock, \
                patch('strategy_matrix.subprocess.run', return_value=type('Result', (), {'returncode': 0})()) as run_mock:
            first = build_matrix(config, matrix, self.root / 'matrix_run', execute=True)
            self.assertEqual(first['candidates'][0]['state'], 'COMPUTATIONAL_PASS')
            second = build_matrix(config, matrix, self.root / 'matrix_run', execute=True, resume=True)
            self.assertEqual(second['candidates'][0]['state'], 'COMPUTATIONAL_PASS')
            self.assertEqual(plan_mock.call_count, 1)
            self.assertEqual(run_mock.call_count, 2)
            self.assertIn('--resume', run_mock.call_args_list[-1][0][0])

    def test_no_igg_only_peak_diagnostic_does_not_require_control_bam(self):
        target = self.root / 'target.bam'; target.write_text('fixture')
        manifest = self.root / 'peaks.tsv'
        manifest.write_text('sample_id\ttarget_bam\tcontrol_bam\tbam_policy\ns1\ttarget.bam\t\tduplicates_retained\n')
        out = self.root / 'peak_diagnostic'
        result = subprocess.run([sys.executable, str(ROOT / 'scripts/peak_diagnostics.py'), '--manifest', str(manifest),
                                 '--out', str(out), '--gsize', '1000', '--native-macs2', '--shapes', 'narrow',
                                 '--backgrounds', 'no_igg_diagnostic'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads((out / 'diagnostics.json').read_text())['candidates'][0]
        self.assertEqual(report['state'], 'PLANNED')
        self.assertNotIn('-c', report['command'])

    def _candidate(self, name, fragments, peaks, own_counts):
        root = self.root / name
        (root / 'plan/consensus').mkdir(parents=True)
        (root / 'plan').mkdir(exist_ok=True)
        (root / 'run').mkdir()
        (root / 'run/workflow_status.json').write_text(json.dumps({'state': 'COMPUTATIONAL_PASS'}))
        (root / 'plan/consensus/master.bed').write_text(peaks)
        fragment_files = []
        for sid in ('s1', 's2'):
            fragment_file = root / ('plan/fragments/' + sid + '/fragments.bed')
            fragment_file.parent.mkdir(parents=True)
            fragment_file.write_text(fragments[sid])
            fragment_files.append((sid, fragment_file))
        with (root / 'plan/fragments.tsv').open('w') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'fragments_bed'])
            for sid, fragment_file in fragment_files: w.writerow([sid, str(fragment_file)])
        (root / 'plan/counts.tsv').write_text(own_counts)
        (root / 'plan/metadata.tsv').write_text('sample_id\tbiological_sample_id\tcondition\treplicates_confirmed\ns1\tbio1\tControl\tTRUE\ns2\tbio2\tControl\tTRUE\n')
        return {'candidate_id': name, 'state': 'COMPUTATIONAL_PASS', 'mapq': 20, 'remove_duplicates': False,
                'fraction': 2/3, 'universe': 'support_core', 'blacklist_mode': 'subtract',
                'artifact_set': 'test', 'plan_dir': str(root / 'plan'), 'run_dir': str(root / 'run')}

    def test_strategy_comparison_uses_one_pooled_evaluation_peak_set(self):
        fragments = {'s1': 'chr1\t1\t5\nchr1\t41\t45\nchr1\t42\t46\n',
                     's2': 'chr1\t2\t6\nchr1\t61\t65\nchr1\t62\t66\n'}
        a = self._candidate('a', fragments, 'chr1\t0\t10\nchr1\t40\t50\n',
                            'peak_id\tchrom\tstart\tend\ts1\ts2\na1\tchr1\t0\t10\t5\t1\na2\tchr1\t40\t50\t1\t5\n')
        b = self._candidate('b', fragments, 'chr1\t20\t30\nchr1\t60\t70\n',
                            'peak_id\tchrom\tstart\tend\ts1\ts2\nb1\tchr1\t20\t30\t1\t1\nb2\tchr1\t60\t70\t5\t5\n')
        matrix = self.root / 'matrix.json'
        matrix.write_text(json.dumps({'state': 'REVIEW_REQUIRED', 'candidates': [a, b]}))
        result = compare_strategies(matrix, self.root / 'comparison')
        self.assertEqual(len(result['pairwise_peak_agreement']), 1)
        self.assertAlmostEqual(result['pairwise_peak_agreement'][0]['jaccard_bp'], 0.0)
        self.assertTrue((self.root / 'comparison/pooled_candidate_peak_union.bed').is_file())
        metrics = {r['candidate_id']: r for r in result['candidate_metrics']}
        self.assertEqual(metrics['a']['mean_replicate_correlation'], -1.0)
        self.assertEqual(metrics['b']['mean_replicate_correlation'], 1.0)
        self.assertAlmostEqual(metrics['a']['common_universe_mean_replicate_correlation'],
                               metrics['b']['common_universe_mean_replicate_correlation'])
        self.assertTrue((self.root / 'comparison/a.common_universe_counts.tsv').is_file())
        self.assertIn('No method is selected automatically', result['recommendation'])

    def test_strategy_comparison_rejects_tampered_candidate_outputs(self):
        artifact = self.root / 'artifact.tsv'; artifact.write_text('original')
        status = {'state': 'COMPUTATIONAL_PASS', 'steps': [{'state': 'COMPUTATIONAL_PASS',
                  'outputs': {str(artifact): __import__('hashlib').sha256(b'original').hexdigest()}}]}
        self.assertEqual(verify_run(status), 'COMPUTATIONAL_PASS')
        artifact.write_text('changed')
        self.assertEqual(verify_run(status), 'OUTPUT_HASH_MISMATCH')

    def test_model_comparison_requires_review_and_matches_intervals(self):
        results = []
        for i, (lfc, direction) in enumerate(((2, 'gain'), (1, 'not_significant'))):
            result = self.root / ('result%d.tsv' % i)
            result.write_text('peak_id\tchrom\tstart\tend\tpadj\tlog2FoldChange\tdirection\n'
                              'p1\tchr1\t0\t100\t0.01\t%s\t%s\np2\tchr1\t200\t300\t0.8\t-1\tnot_significant\n' % (lfc, direction))
            review = self.root / ('review%d.json' % i)
            review.write_text(json.dumps({'reviewer': 'reviewer', 'reviewed_at': '2026-10-04', 'reason': 'fixture',
                'decisions': {'upstream_accepted': True, 'peaks_accepted': True, 'replicates_confirmed': True,
                              'simple_design_accepted': True},
                'inputs': {'metadata_sha256': 'a' * 64},
                'model': {'normalization': 'conventional', 'numerator': 'KD', 'denominator': 'Control',
                          'design': '~ condition', 'alpha': .05, 'abs_log2fc': 1, 'min_total_count': 10}}))
            results.append((result, review))
        manifest = self.root / 'models.tsv'
        with manifest.open('w') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['candidate_id', 'normalization', 'result_tsv', 'review_json'])
            for i, (result, review) in enumerate(results): w.writerow(['m%d' % i, 'conventional', result, review])
        pairs = compare_models(manifest, self.root / 'model_comparison')
        self.assertEqual(pairs[0]['matched_regions'], 2)
        self.assertEqual(pairs[0]['significance_status_disagreements'], 1)


if __name__ == '__main__':
    unittest.main()
