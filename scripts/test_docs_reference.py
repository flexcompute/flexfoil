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

class CargoTargetTests(unittest.TestCase):
    def test_build_copies_resolved_cargo_target_not_root_target(self):
        import json
        import tempfile
        from pathlib import Path
        from unittest.mock import patch
        import build_docs_reference as docs
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            site = root / 'docs-site'
            (site / 'docs').mkdir(parents=True)
            (root / 'packages/flexfoil-python/src/flexfoil').mkdir(parents=True)
            (root / 'crates').mkdir()
            target = root / 'custom-target/doc/rustfoil_xfoil'
            target.mkdir(parents=True)
            (target / 'index.html').write_text('current Rust API')
            with patch.multiple(docs, ROOT=root, SITE=site, GENERATED=site/'src/generated', API=site/'docs/api'), patch.object(docs.subprocess, 'run'), patch.object(docs.subprocess, 'check_output', return_value=json.dumps({'target_directory': str(root/'custom-target')})):
                docs.build()
            self.assertEqual((site/'static/api/rust/rustfoil_xfoil/index.html').read_text(), 'current Rust API')


if __name__ == '__main__':
    unittest.main()
