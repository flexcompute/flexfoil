import assert from 'node:assert/strict';
import test from 'node:test';
import equations from './equations.mjs';
import fs from 'node:fs';
import path from 'node:path';
import os from 'node:os';
import {pathToFileURL} from 'node:url';
const page = path.resolve('docs/internals/compressible-boundary-layer.mdx');
const label = /\\label\{([^}]+)\}/.exec(fs.readFileSync(page, 'utf8'))[1];
test('labels produce a numbered equation and a same-page link', () => {
  const tree = {children: [{type: 'math', value: `U=1\\label{${label}}`}, {type: 'inlineMath', value: `\\eqref{${label}}`}]};
  equations()(tree, {path: page});
  assert.equal(tree.children[0].name, 'a');
  assert.match(tree.children[1].value, /\\tag\{1\}/);
  assert.equal(tree.children[2].type, 'link');
  assert.equal(tree.children[2].url, '#'+label.replace(/[^A-Za-z0-9_-]/g, '-'));
});
test('cross-page links include the docs deployment path', () => {
  const tree = {children: [{type:'inlineMath', value:`\\eqref{${label}}`}]};
  equations()(tree, {path:path.resolve('docs/internals/compressible-forces.mdx')});
  assert.match(tree.children[0].url, /^\/docs\/internals\/compressible-boundary-layer#/);
});
test('unknown equation references fail compilation', () => {
  assert.throws(() => equations()({children:[{type:'inlineMath',value:'\\eqref{eq:missing-label}'}]}, {path:page}), /unknown equation label/);
});


test('registry resolves paths containing spaces and percent characters', async () => {
  const root=fs.mkdtempSync(path.join(os.tmpdir(),'flexfoil docs 100% '));
  try {
    const remark=path.join(root,'docs-site/src/remark');
    const docs=path.join(root,'docs-site/docs');
    fs.mkdirSync(remark,{recursive:true});fs.mkdirSync(docs);
    const module=path.join(remark,'equations.mjs');
    fs.copyFileSync(path.resolve('src/remark/equations.mjs'),module);
    const page=path.join(docs,'sample.mdx');
    fs.writeFileSync(page,'---\ntitle: Sample\n---\n\n$$\nU=1\\label{eq:encoded-path}\n$$\n');
    const plugin=(await import(pathToFileURL(module).href)).default;
    const tree={children:[{type:'inlineMath',value:'\\eqref{eq:encoded-path}'}]};
    plugin()(tree,{path:fs.realpathSync(page)});
    assert.equal(tree.children[0].url,'#eq-encoded-path');
  } finally {fs.rmSync(root,{recursive:true,force:true});}
});
