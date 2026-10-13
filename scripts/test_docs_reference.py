"""Documented names must exist and excerpts must come from source."""
import unittest
from build_docs_reference import ROOT, extract


class SourceExtractionTests(unittest.TestCase):
    def test_python_function_and_method_are_real_source(self):
        for path, name in [('scripts/build_solver_reference.py', 'build'),
                           ('packages/flexfoil-python/src/flexfoil/airfoil.py', 'Airfoil.solve')]:
            snippet = extract(path, name)
            lines = (ROOT / path).read_text().splitlines()[snippet['start'] - 1:snippet['end']]
            self.assertIn('def ' + name.split('.')[-1] + '(', '\n'.join(lines))
            self.assertTrue(snippet['code'].startswith('def '))

    def test_missing_documented_name_is_fatal(self):
        with self.assertRaises(StopIteration):
            extract('scripts/build_solver_reference.py', 'removed_primitive')
        with self.assertRaises(ValueError):
            extract('crates/rustfoil-xfoil/src/forces.rs', 'removed_primitive')

    def test_rust_and_widget_have_located_code(self):
        for path, name in [('crates/rustfoil-xfoil/src/forces.rs', 'compute_panel_forces_from_gamma'),
                           ('docs-site/src/components/SolverGate/Inspector.jsx', 'GateInspector')]:
            snippet = extract(path, name)
            self.assertIn(name, snippet['code'])
            self.assertEqual(snippet['code'].splitlines()[-1], '}')
            self.assertGreater(snippet['end'], snippet['start'])


if __name__ == '__main__':
    unittest.main()
