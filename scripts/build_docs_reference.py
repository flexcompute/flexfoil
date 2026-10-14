"""Build-time source extraction, public Python API and documentation coverage.

Python definitions use ast (as in FlightModule); Rust/JS snippets use this repo's
column-zero closing-brace convention. Rustdoc generates the canonical Rust API.
"""
import ast
import json
from pathlib import Path
import re
import shutil
import subprocess
import textwrap

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / 'docs-site'
GENERATED = SITE / 'src/generated'
API = SITE / 'docs/api'


# shortcut: Rust/JS extraction uses column-zero braces; use compiler parsing if that source convention changes.
def extract(path, name):
    text = (ROOT / path).read_text()
    lines = text.splitlines()
    if path.endswith('.py'):
        node = ast.parse(text)
        for part in name.split('.'):
            node = next(n for n in ast.iter_child_nodes(node)
                        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and n.name == part)
        start, end = node.lineno, node.end_lineno
    else:
        keyword = r'(?:fn|struct|enum)' if path.endswith('.rs') else r'function'
        pattern = rf'^(?:pub(?:\([^)]*\))?\s+|export\s+(?:default\s+)?)?{keyword}\s+{re.escape(name)}\b'
        matches = [i for i, line in enumerate(lines) if re.match(pattern, line)]
        if len(matches) != 1:
            raise ValueError(f'{path}: expected one documented definition {name}, found {len(matches)}')
        start = matches[0] + 1
        end = next(i + 1 for i in range(start, len(lines)) if lines[i] == '}')
    return dict(path=path, name=name, start=start, end=end, code=textwrap.dedent('\n'.join(lines[start - 1:end])))


def definitions(tree, prefix=''):
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and (not node.name.startswith('_') or node.name == '__init__'):
            name = prefix + node.name
            yield name, node
            if isinstance(node, ast.ClassDef):
                yield from definitions(node, name + '.')


def build():
    GENERATED.mkdir(parents=True, exist_ok=True)
    API.mkdir(parents=True, exist_ok=True)
    for old in API.glob('*.md'):
        old.unlink()
    snippets = {}
    pages = {}
    for page in (SITE / 'docs').rglob('*.mdx'):
        for path, name in re.findall(r'<SourceCode\s+path="([^"]+)"\s+name="([^"]+)"\s*/>', page.read_text()):
            snippets[path + ':' + name] = extract(path, name)
            pages[(path, name)] = '/' + page.relative_to(SITE / 'docs').with_suffix('').as_posix()
    (GENERATED / 'source.json').write_text(json.dumps(snippets, indent=2))
    coverage = []
    modules = []
    package = ROOT / 'packages/flexfoil-python/src/flexfoil'
    for source in sorted(package.rglob('*.py')):
        relative = source.relative_to(package)
        if any(part.startswith('_') and part not in ('__init__.py', '__main__.py') for part in relative.parts):
            continue
        path = source.relative_to(ROOT).as_posix()
        module = 'flexfoil.' + relative.with_suffix('').as_posix().replace('/', '.')
        module = module.removesuffix('.__init__').removesuffix('.')
        tree = ast.parse(source.read_text())
        items = list(definitions(tree))
        if not items:
            continue
        modules.append(module)
        body = [f'---\ntitle: {json.dumps(module)}\n---\n', 'Generated from AST signatures and docstrings; no package code is imported.\n']
        for name, node in items:
            if isinstance(node, ast.ClassDef):
                signature = 'class ' + node.name + '(' + ', '.join(ast.unparse(base) for base in node.bases) + ')'
            else:
                signature = ('async def ' if isinstance(node, ast.AsyncFunctionDef) else 'def ') + node.name + '(' + ast.unparse(node.args) + ')'
                if node.returns:
                    signature += ' -> ' + ast.unparse(node.returns)
            page = pages.get((path, name))
            label = module + '.' + name
            coverage.append(dict(name=label, path=path, line=node.lineno, page=page))
            body += [f'## `{name}`\n', f'```python\n{signature}\n```\n', f'Source: `{path}:{node.lineno}`\n']
            if page:
                body.append(f'[Hand-written primitive page]({page})\n')
            body.append('~~~~text\n' + (ast.get_docstring(node) or 'No docstring.') + '\n~~~~\n')
            if isinstance(node, ast.ClassDef):
                for field in node.body:
                    if isinstance(field, ast.AnnAssign) and isinstance(field.target, ast.Name) and not field.target.id.startswith('_'):
                        body.append(f'```python\n{ast.unparse(field)}\n```\n')
        (API / (module + '.md')).write_text('\n'.join(body))
    rust = []
    for source in sorted((ROOT / 'crates').rglob('*.rs')):
        path = source.relative_to(ROOT).as_posix()
        for line, text in enumerate(source.read_text().splitlines(), 1):
            match = re.match(r'\s*pub\s+(?:(?:async|unsafe|const)\s+)?(?:fn|struct|enum|trait|type)\s+(\w+)', text)
            if match:
                name = match.group(1)
                rust.append(dict(name=path + ':' + name, path=path, line=line, page=pages.get((path, name))))
    all_items = coverage + rust
    (GENERATED / 'coverage.json').write_text(json.dumps(all_items, indent=2))
    (API / '_category_.json').write_text(json.dumps({'label': 'Generated API', 'position': 9}))
    (API / 'index.md').write_text('---\ntitle: Generated API\nslug: /api\n---\n\nPython signatures/docstrings are generated on every build with stdlib AST. Rustdoc generates the canonical public Rust API.\n\n' + ''.join(f'- [{module}](./{module}.md)\n' for module in modules) + '\n<a href="/docs/api/rust/rustfoil_xfoil/index.html">Rust API and linked workspace crates</a>\n\n## Primitive pages\n\n' + ''.join(f'- [`{item["name"]}`](../{item["page"].lstrip("/")}.mdx)\n' for item in all_items if item['page']) + '\n[Hand-written page coverage](./coverage.md)\n')
    (API / 'coverage.md').write_text('---\ntitle: Hand-written page coverage\n---\n\nGenerated on every build. Python entries are public AST declarations (including constructors). Rust entries are public source declarations, possibly inside private modules; rustdoc is authoritative for export visibility. A link is counted only when a hand-written page extracts that named primitive. Existing overview pages have not been marked as complete primitive documentation.\n\n' + f'{sum(bool(x["page"]) for x in all_items)} of {len(all_items)} declarations have a mapped hand-written page.\n\n' + ''.join(f'- `{x["name"]}` — `{x["path"]}:{x["line"]}`\n' for x in all_items if not x['page']))
    subprocess.run(['cargo', 'doc', '--workspace', '--no-deps', '--quiet'], cwd=ROOT, check=True)
    destination = SITE / 'static/api/rust'
    if destination.exists():
        shutil.rmtree(destination)
    metadata = json.loads(subprocess.check_output(['cargo', 'metadata', '--no-deps', '--format-version', '1'], cwd=ROOT))
    shutil.copytree(Path(metadata['target_directory']) / 'doc', destination)
    assert (destination / 'rustfoil_xfoil/index.html').exists()
    print(f'{len(snippets)} extracted primitives; {len(coverage)} Python API entries; {len(rust)} Rust coverage candidates')


if __name__ == '__main__':
    build()
