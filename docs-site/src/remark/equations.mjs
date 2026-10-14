// Numbered, linkable equations for every page (remark plugin, runs after remark-math).
//
//   $$ ... \label{eq:air-states} $$   every display equation is numbered (1), (2), ... per page; \label names it
//   $\eqref{eq:air-states}$            inline: becomes a link "(n)" to that equation, on this page or any other
//
// Numbers come from a registry built by scanning every page in docs/ with the same rule, so a reference to
// another page resolves to that page's number. An unknown label fails the build, as a missing code snippet does.
import fs from 'node:fs';
import path from 'node:path';

const DOCS = path.resolve(path.dirname(new URL(import.meta.url).pathname), '../../docs');
const BLOCK = /^\$\$\s*\n([\s\S]*?)\n\$\$/gm;
const LABEL = /\\label\{([^}]+)\}/;

function walk(dir) {
  return fs.readdirSync(dir, {withFileTypes: true}).flatMap(e =>
    e.isDirectory() ? walk(path.join(dir, e.name)) : /\.mdx?$/.test(e.name) ? [path.join(dir, e.name)] : []);
}

function route(file, text) {
  const slug = /^---[\s\S]*?^slug:\s*(\S+)/m.exec(text);
  if (slug) return slug[1];
  const rel = path.relative(DOCS, file).replace(/\.mdx?$/, '').split(path.sep).join('/');
  return '/' + rel.replace(/(^|\/)index$/, '');
}

let REGISTRY = null;
function registry() {
  if (REGISTRY) return REGISTRY;
  REGISTRY = {};
  for (const f of walk(DOCS)) {
    const text = fs.readFileSync(f, 'utf8');
    const title = (/^title:\s*"?(.+?)"?\s*$/m.exec(text) || [])[1] || path.basename(f);
    let n = 0, m;
    BLOCK.lastIndex = 0;
    while ((m = BLOCK.exec(text))) {
      n += 1;
      const label = (LABEL.exec(m[1]) || [])[1];
      if (label) {
        if (REGISTRY[label]) throw new Error(`equation label ${label} defined twice (${REGISTRY[label].file}, ${f})`);
        REGISTRY[label] = {n, route: route(f, text), title, file: f};
      }
    }
  }
  return REGISTRY;
}

const anchorId = (label, n) => (label || `eq-${n}`).replace(/[^A-Za-z0-9_-]/g, '-');

export default function equations() {
  return (tree, file) => {
    const reg = registry();
    const here = file.path ? route(file.path, fs.readFileSync(file.path, 'utf8')) : '';
    const mdx = /\.mdx$/.test(file.path || '');
    let n = 0;
    const visit = (node, parent) => {
      if (node.type === 'math') {
        n += 1;
        const label = (LABEL.exec(node.value) || [])[1];
        node.value = node.value.replace(LABEL, '') + ` \\tag{${n}}`;
        // remark-math copies the TeX into data.hChildren at parse time, and that copy is what gets rendered
        // (pre > code.math-display > text): replace the text, keep the elements rehype-mathjax looks for
        const retext = h => { if (h.type === 'text') h.value = node.value; (h.children || []).forEach(retext); };
        (node.data && node.data.hChildren || []).forEach(retext);
        const id = anchorId(label, n);
        const anchor = mdx
          ? {type: 'mdxJsxFlowElement', name: 'a', attributes: [{type: 'mdxJsxAttribute', name: 'id', value: id},
              {type: 'mdxJsxAttribute', name: 'className', value: 'eq-anchor'}], children: []}
          : {type: 'html', value: `<a id="${id}" class="eq-anchor"></a>`};
        parent.children.splice(parent.children.indexOf(node), 0, anchor);
        return 1; // skip the inserted anchor
      }
      if (node.type === 'inlineMath') {
        const ref = /^\s*\\(?:eq)?ref\{([^}]+)\}\s*$/.exec(node.value);
        if (ref) {
          const e = reg[ref[1]];
          if (!e) throw new Error(`${file.path}: unknown equation label ${ref[1]}`);
          const same = e.route === here;
          Object.assign(node, {type: 'link', url: `${same ? '' : '/docs' + e.route}#${anchorId(ref[1], e.n)}`,
            title: same ? `equation ${e.n}` : `${e.title}, equation ${e.n}`,
            children: [{type: 'text', value: same ? `(${e.n})` : `(${e.n}) in ${e.title}`}]});
          delete node.value;
          delete node.data; // remark-math hints (hName, hChildren) would otherwise still render it as maths
        }
      }
      if (node.children) for (let i = 0; i < node.children.length; i++) i += visit(node.children[i], node) || 0;
      return 0;
    };
    visit(tree, null);
  };
}
