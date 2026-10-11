import { describe, expect, it } from 'vitest';
import { parseInverseTarget } from './inverseTargetImport';

describe('inverse target import', () => {
  it('preserves exact Cp targets and their segment', () => {
    expect(parseInverseTarget('# Upper\nx,cp\n.2,-4\n.6,-.123456789\n', 'cp'))
      .toEqual({ x: [.2, .6], values: [-4, -.123456789] });
  });
  it('accepts headerless whitespace and velocity headers', () => {
    expect(parseInverseTarget('x/c ue\n0 1.2\n1 .8', 'velocity').values).toEqual([1.2, .8]);
    expect(parseInverseTarget('0 -1\n1 0', 'cp').x).toEqual([0, 1]);
  });
  it.each(['x,ue\n0,1\n1,1', '0,NaN\n1,0', '0,Infinity\n1,0',
    '0,1,2\n1,1', '0,1\n0,2', '.5,1\n.2,1', '-.1,1\n1,1',
    '0,1\n1.1,1', '0,1', '', '0,\n1,0'])('rejects invalid input: %s', (text) => {
    expect(() => parseInverseTarget(text, 'cp')).toThrow();
  });
});
