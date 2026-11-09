"""Test production collector/journal functions with synthetic fit inputs, without ROOT.

AST loading excludes site imports and unrelated ROOT readers, not function bodies.
"""
import ast
import datetime
import hashlib
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[1]


def functions_and_constants(path, namespace, functions=None):
    tree = ast.parse(path.read_text())
    tree.body = [node for node in tree.body
                 if isinstance(node, ast.FunctionDef) and (functions is None or node.name in functions)
                 or isinstance(node, ast.Assign) and not any(
                     isinstance(target, ast.Name) and target.id == 'COLLECT' for target in node.targets)]
    exec(compile(tree, str(path), 'exec'), namespace)
    return namespace


class ResolutionJournalTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        def fail(message):
            raise RuntimeError(message)
        runner = SimpleNamespace(EXIT_OK=0, EXIT_FAIL=1, EXIT_NOT_CONVERGED=2,
                                 EXIT_CHECK=3, EXIT_SKIP=4, warn=lambda *a: None, fail=fail)
        namespace = dict(config=SimpleNamespace(OUTPUT_DIR=self.root), runner=runner,
                         Path=Path, datetime=datetime, hashlib=hashlib, json=json)
        self.conv = SimpleNamespace(**functions_and_constants(
            REPO / 'calibration/scripts/calib_convergence.py', namespace))
        self.decode = self.root / 'decode.root'
        self.decode.write_text('synthetic decode identity')
        self.target = self.root / 'param'
        self.target.write_text('before')
        self.current = {str(i): ['0'] * 8 + ['0.2'] for i in range(1, 17)}
        self.pf = SimpleNamespace(path=self.target, get_data_dict=lambda **kw: self.current)
        self.args = SimpleNamespace(all=False, chamber='BcIn', run=1, suffix='Pi')
        self.ns = functions_and_constants(
            REPO / 'calibration/scripts/update_param.py',
            dict(np=np, math=math, conv=self.conv, runner=runner),
            {'collect_resolution', 'resolution_apply_key', 'pull_width_ok'})
        self.ns.update(decode=lambda *a: self.decode,
                       pull_groups=lambda a: [('BLC1', 1.2, 1.0, 1)],
                       report_resolution=lambda *a, **kw: None)

    def collect(self):
        return self.ns['collect_resolution'](self.args, self.pf)

    def apply(self, state='written'):
        data, states, pending = self.collect()
        self.assertEqual(data['1'], ['0.24'])
        before = self.conv.file_sha256(self.target)
        self.target.write_text('after')
        after = self.conv.file_sha256(self.target)
        rec = pending[0]
        self.conv.record(rec['path'], rec['group'], rec['value'], state,
                         self.conv.details_with_apply(rec['apply_key'], before, after))
        for key, value in data.items():
            self.current[key][8] = value[0]
        return rec

    def test_completed_write_is_not_applied_twice(self):
        self.apply()
        self.assertEqual(self.collect(), ({}, [self.conv.SKIPPED], []))

    def test_interruption_after_write_is_not_applied_twice(self):
        self.apply(self.conv.PENDING)
        self.assertEqual(self.collect(), ({}, [self.conv.SKIPPED], []))

    def test_new_decode_allows_next_pass(self):
        self.apply()
        self.decode.write_text('a different decode with a different size')
        data, states, pending = self.collect()
        self.assertEqual(data['1'], ['0.288'])
        self.assertEqual(states, [self.conv.WRITTEN])

    def test_legacy_journal_requires_new_decode(self):
        rec = self.apply()
        text = rec['path'].read_text().replace('reso-v2:', '')
        rec['path'].write_text(text)
        with self.assertRaisesRegex(RuntimeError, 'legacy apply journal'):
            self.collect()

    def test_pending_before_write_can_retry(self):
        data, _, pending = self.collect()
        rec = pending[0]
        self.conv.record(rec['path'], rec['group'], rec['value'], self.conv.PENDING,
                         self.conv.details_with_apply(rec['apply_key'],
                             self.conv.file_sha256(self.target), 'expected-after'))
        self.assertEqual(self.collect()[0], data)


if __name__ == '__main__':
    unittest.main()
