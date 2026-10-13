import React from 'react';
import katex from 'katex';
import 'katex/dist/katex.min.css';

export default function Equation({ math }) {
  return <div style={{ overflowX: 'auto' }} dangerouslySetInnerHTML={{ __html: katex.renderToString(math, { displayMode: true, throwOnError: true, trust: false }) }} />;
}
