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
        self.assertEqual(len({r['candidate_id'] for r in result}), 16)
        self.assertEqual({r['artifact_set'] for r in result}, {'narrow_igg', 'broad_igg'})
        with self.assertRaises(ValueError):
            expand({**spec, 'axes': {'universe': ['unknown']}})

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

    def _candidate(self, name, fragments, peaks):
        root = self.root / name
        (root / 'plan/consensus').mkdir(parents=True)
        (root / 'plan').mkdir(exist_ok=True)
        (root / 'run').mkdir()
        (root / 'run/workflow_status.json').write_text(json.dumps({'state': 'COMPUTATIONAL_PASS'}))
        (root / 'plan/consensus/master.bed').write_text(peaks)
        fragment_file = root / 'plan/fragments/s1/fragments.bed'
        fragment_file.parent.mkdir(parents=True)
        fragment_file.write_text(fragments)
        with (root / 'plan/fragments.tsv').open('w') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['sample_id', 'fragments_bed']); w.writerow(['s1', str(fragment_file)])
        (root / 'plan/counts.tsv').write_text('peak_id\tchrom\tstart\tend\ts1\np1\tchr1\t0\t30\t2\n')
        (root / 'plan/metadata.tsv').write_text('sample_id\tbiological_sample_id\tcondition\treplicates_confirmed\ns1\tbio1\tControl\tFALSE\n')
        (root / 'plan/counts.tsv').write_text('peak_id\tchrom\tstart\tend\ts1\np1\tchr1\t0\t30\t2\n')
        (root / 'plan/metadata.tsv').write_text('sample_id\tbiological_sample_id\tcondition\treplicates_confirmed\ns1\tbio1\tControl\tFALSE\n')
        return {'candidate_id': name, 'state': 'COMPUTATIONAL_PASS', 'mapq': 20, 'remove_duplicates': False,
                'fraction': 2/3, 'universe': 'support_core', 'blacklist_mode': 'subtract',
                'artifact_set': 'test', 'plan_dir': str(root / 'plan'), 'run_dir': str(root / 'run')}

    def test_strategy_comparison_uses_one_pooled_evaluation_peak_set(self):
        a = self._candidate('a', 'chr1\t0\t10\nchr1\t20\t30\n', 'chr1\t0\t20\n')
        b = self._candidate('b', 'chr1\t0\t10\nchr1\t20\t30\n', 'chr1\t10\t30\n')
        matrix = self.root / 'matrix.json'
        matrix.write_text(json.dumps({'state': 'REVIEW_REQUIRED', 'candidates': [a, b]}))
        result = compare_strategies(matrix, self.root / 'comparison')
        self.assertEqual(len(result['pairwise_peak_agreement']), 1)
        self.assertAlmostEqual(result['pairwise_peak_agreement'][0]['jaccard_bp'], 1/3)
        self.assertTrue((self.root / 'comparison/pooled_candidate_peak_union.bed').is_file())
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
