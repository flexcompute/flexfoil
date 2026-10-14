import assert from 'node:assert/strict';
import test from 'node:test';
import equations from './equations.mjs';
import fs from 'node:fs';
import path from 'node:path';
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
