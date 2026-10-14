"""Negative checks for the gate's failure contract."""
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from solver_gate import require_close, run_reference


class GateFailureTests(unittest.TestCase):
    def test_nonfinite_and_wrong_values_cannot_pass(self):
        for value in (math.nan, math.inf, -math.inf, 1.01):
            with self.subTest(value=value), self.assertRaises(AssertionError):
                require_close(value, 1., .001, 'coefficient')

    def test_successful_process_with_no_converged_row_fails(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            geometry = root / 'input.dat'
            geometry.write_text('test\n1 .01\n0 0\n1 -.01\n')
            (root / 'polar.txt').write_text('1 1 Reynolds number fixed\nxtrf = 1.000 (top) 1.000 (bottom)\nMach = 0.000 Re = 1.000 e 6 Ncrit = 9.000 9.000\nalpha CL CD\n------ ----- -----\n')
            with patch('solver_gate.subprocess.run', return_value=subprocess.CompletedProcess([], 0, '', '')):
                with self.assertRaisesRegex(AssertionError, 'exactly one converged'):
                    run_reference(Path('unused'), geometry, 4., 1e6, 0., 9., (1., 1.), root)

    def test_reference_process_failure_is_fatal(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            geometry = root / 'input.dat'
            geometry.write_text('test\n1 .01\n0 0\n1 -.01\n')
            with patch('solver_gate.subprocess.run', return_value=subprocess.CompletedProcess([], 2, '', 'failed')):
                with self.assertRaisesRegex(AssertionError, 'XFOIL exit 2'):
                    run_reference(Path('unused'), geometry, 4., 1e6, 0., 9., (1., 1.), root)


if __name__ == '__main__':
    unittest.main()
