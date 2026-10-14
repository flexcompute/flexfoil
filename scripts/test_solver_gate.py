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

    def test_drag_components_use_final_iteration_not_polar_pressure(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            geometry = root / 'input.dat'
            geometry.write_text('test\n1 0\n0 0\n1 -.1\n')
            (root / 'polar.txt').write_text('1 1 Reynolds number fixed\nxtrf = 1.000 (top) 1.000 (bottom)\nMach = 0.000 Re = 1.000 e 6 Ncrit = 9.000 9.000\nalpha CL CD CDp CM XtrT XtrB\n------ ----- -----\n4 .5 .00716 -.00010 -.05 .3 1 0 0\n')
            (root / 'bl.txt').write_text('s x y ue ds th cf hk\n0 1 0 1 .001 .0001 .002 2.5\n1 0 0 1 .001 .0001 .002 2.5\n2 1 -.1 1 .001 .0001 .002 2.5\n')
            with patch('solver_gate.subprocess.run', return_value=subprocess.CompletedProcess([], 0, 'CDf = 0.00600 CDp = 0.00116\nCDf = 0.00465 CDp = 0.00250', '')):
                result, _ = run_reference(Path('unused'), geometry, 4., 1e6, 0., 9., (1., 1.), root)
            self.assertEqual(result['cdp'], -.00010)
            self.assertEqual(result['cd_friction'], .00465)
            self.assertEqual(result['cd_pressure'], .00250)
            with patch('solver_gate.subprocess.run', return_value=subprocess.CompletedProcess([], 0, '', '')):
                with self.assertRaisesRegex(AssertionError, 'drag components missing'):
                    run_reference(Path('unused'), geometry, 4., 1e6, 0., 9., (1., 1.), root)


if __name__ == '__main__':
    unittest.main()
