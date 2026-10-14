import React from 'react';
import CodeBlock from '@theme/CodeBlock';
import snippets from '@site/src/generated/source.json';

export default function SourceCode({ path, name }) {
  const snippet = snippets[`${path}:${name}`];
  if (!snippet) throw new Error(`Missing build-extracted source: ${path}:${name}`);
  const language = path.endsWith('.py') ? 'python' : path.endsWith('.rs') ? 'rust' : 'jsx';
  return <CodeBlock language={language} title={`${path}: lines ${snippet.start}–${snippet.end}`} showLineNumbers>{snippet.code}</CodeBlock>;
}
