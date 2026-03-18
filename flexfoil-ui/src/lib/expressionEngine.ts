/**
 * Safe algebraic expression evaluator for custom computed columns.
 *
 * Supports: +, -, *, /, ^, % operators, parentheses, unary minus,
 * numeric literals, column references, and a whitelist of math functions.
 *
 * Uses a recursive-descent parser — no eval(), no Function(), no external deps.
 */

import type { RunRow } from '../types';

// ── Whitelisted functions ──

const MATH_FUNCTIONS: Record<string, (...args: number[]) => number> = {
  abs: Math.abs,
  sqrt: Math.sqrt,
  log: Math.log,
  log10: Math.log10,
  log2: Math.log2,
  exp: Math.exp,
  sin: Math.sin,
  cos: Math.cos,
  tan: Math.tan,
  asin: Math.asin,
  acos: Math.acos,
  atan: Math.atan,
  ceil: Math.ceil,
  floor: Math.floor,
  round: Math.round,
  min: Math.min,
  max: Math.max,
  sign: Math.sign,
};

const FUNCTION_NAMES = new Set(Object.keys(MATH_FUNCTIONS));

const CONSTANTS: Record<string, number> = {
  pi: Math.PI,
  e: Math.E,
};

// ── Tokenizer ──

type TokenKind =
  | 'number'
  | 'identifier'
  | 'operator'
  | 'lparen'
  | 'rparen'
  | 'comma'
  | 'eof';

interface Token {
  kind: TokenKind;
  value: string;
  pos: number;
}

function tokenize(input: string): Token[] {
  const tokens: Token[] = [];
  let i = 0;

  while (i < input.length) {
    const ch = input[i];

    if (/\s/.test(ch)) {
      i++;
      continue;
    }

    if (/[0-9.]/.test(ch)) {
      let num = '';
      const start = i;
      while (i < input.length && /[0-9.eE]/.test(input[i])) {
        if ((input[i] === 'e' || input[i] === 'E') && i + 1 < input.length && (input[i + 1] === '+' || input[i + 1] === '-')) {
          num += input[i] + input[i + 1];
          i += 2;
        } else {
          num += input[i];
          i++;
        }
      }
      tokens.push({ kind: 'number', value: num, pos: start });
      continue;
    }

    if (/[a-zA-Z_]/.test(ch)) {
      let id = '';
      const start = i;
      while (i < input.length && /[a-zA-Z0-9_]/.test(input[i])) {
        id += input[i];
        i++;
      }
      tokens.push({ kind: 'identifier', value: id, pos: start });
      continue;
    }

    if ('+-*/%^'.includes(ch)) {
      tokens.push({ kind: 'operator', value: ch, pos: i });
      i++;
      continue;
    }

    if (ch === '(') {
      tokens.push({ kind: 'lparen', value: '(', pos: i });
      i++;
      continue;
    }
    if (ch === ')') {
      tokens.push({ kind: 'rparen', value: ')', pos: i });
      i++;
      continue;
    }
    if (ch === ',') {
      tokens.push({ kind: 'comma', value: ',', pos: i });
      i++;
      continue;
    }

    throw new ExpressionError(`Unexpected character '${ch}'`, i);
  }

  tokens.push({ kind: 'eof', value: '', pos: input.length });
  return tokens;
}

// ── Parser (recursive descent, builds an AST) ──

type ASTNode =
  | { type: 'number'; value: number }
  | { type: 'column'; name: string }
  | { type: 'unary'; op: string; operand: ASTNode }
  | { type: 'binary'; op: string; left: ASTNode; right: ASTNode }
  | { type: 'call'; name: string; args: ASTNode[] };

class ExpressionError extends Error {
  pos: number;
  constructor(message: string, pos: number) {
    super(message);
    this.name = 'ExpressionError';
    this.pos = pos;
  }
}

class Parser {
  private tokens: Token[];
  private pos = 0;
  private columns: Set<string>;

  constructor(tokens: Token[], columns: Set<string>) {
    this.tokens = tokens;
    this.columns = columns;
  }

  private peek(): Token {
    return this.tokens[this.pos];
  }

  private advance(): Token {
    const t = this.tokens[this.pos];
    this.pos++;
    return t;
  }

  private expect(kind: TokenKind, value?: string): Token {
    const t = this.peek();
    if (t.kind !== kind || (value !== undefined && t.value !== value)) {
      throw new ExpressionError(
        `Expected ${value ?? kind}, got '${t.value}'`,
        t.pos,
      );
    }
    return this.advance();
  }

  parse(): ASTNode {
    const node = this.parseExpression();
    if (this.peek().kind !== 'eof') {
      throw new ExpressionError(
        `Unexpected token '${this.peek().value}'`,
        this.peek().pos,
      );
    }
    return node;
  }

  // expression → term (('+' | '-') term)*
  private parseExpression(): ASTNode {
    let left = this.parseTerm();
    while (this.peek().kind === 'operator' && (this.peek().value === '+' || this.peek().value === '-')) {
      const op = this.advance().value;
      const right = this.parseTerm();
      left = { type: 'binary', op, left, right };
    }
    return left;
  }

  // term → power (('*' | '/' | '%') power)*
  private parseTerm(): ASTNode {
    let left = this.parsePower();
    while (this.peek().kind === 'operator' && ('*/%'.includes(this.peek().value))) {
      const op = this.advance().value;
      const right = this.parsePower();
      left = { type: 'binary', op, left, right };
    }
    return left;
  }

  // power → unary ('^' power)?  (right-associative)
  private parsePower(): ASTNode {
    const base = this.parseUnary();
    if (this.peek().kind === 'operator' && this.peek().value === '^') {
      this.advance();
      const exp = this.parsePower();
      return { type: 'binary', op: '^', left: base, right: exp };
    }
    return base;
  }

  // unary → ('-' | '+') unary | atom
  private parseUnary(): ASTNode {
    if (this.peek().kind === 'operator' && (this.peek().value === '-' || this.peek().value === '+')) {
      const op = this.advance().value;
      const operand = this.parseUnary();
      if (op === '+') return operand;
      return { type: 'unary', op: '-', operand };
    }
    return this.parseAtom();
  }

  // atom → NUMBER | IDENTIFIER | function_call | '(' expression ')'
  private parseAtom(): ASTNode {
    const t = this.peek();

    if (t.kind === 'number') {
      this.advance();
      const val = parseFloat(t.value);
      if (isNaN(val)) throw new ExpressionError(`Invalid number '${t.value}'`, t.pos);
      return { type: 'number', value: val };
    }

    if (t.kind === 'identifier') {
      this.advance();
      const name = t.value;

      // Function call
      if (this.peek().kind === 'lparen' && FUNCTION_NAMES.has(name)) {
        this.advance(); // consume '('
        const args: ASTNode[] = [];
        if (this.peek().kind !== 'rparen') {
          args.push(this.parseExpression());
          while (this.peek().kind === 'comma') {
            this.advance();
            args.push(this.parseExpression());
          }
        }
        this.expect('rparen');
        return { type: 'call', name, args };
      }

      // Constant
      if (name in CONSTANTS) {
        return { type: 'number', value: CONSTANTS[name] };
      }

      // Column reference
      if (this.columns.has(name)) {
        return { type: 'column', name };
      }

      throw new ExpressionError(
        `Unknown identifier '${name}'. Available columns: ${[...this.columns].join(', ')}`,
        t.pos,
      );
    }

    if (t.kind === 'lparen') {
      this.advance();
      const node = this.parseExpression();
      this.expect('rparen');
      return node;
    }

    throw new ExpressionError(`Unexpected token '${t.value}'`, t.pos);
  }
}

// ── Evaluator ──

function evalAST(node: ASTNode, vars: Record<string, number>): number {
  switch (node.type) {
    case 'number':
      return node.value;
    case 'column': {
      const v = vars[node.name];
      return v ?? NaN;
    }
    case 'unary':
      return -evalAST(node.operand, vars);
    case 'binary': {
      const l = evalAST(node.left, vars);
      const r = evalAST(node.right, vars);
      switch (node.op) {
        case '+': return l + r;
        case '-': return l - r;
        case '*': return l * r;
        case '/': return r === 0 ? NaN : l / r;
        case '%': return r === 0 ? NaN : l % r;
        case '^': return Math.pow(l, r);
        default: return NaN;
      }
    }
    case 'call': {
      const fn = MATH_FUNCTIONS[node.name];
      if (!fn) return NaN;
      const argVals = node.args.map(a => evalAST(a, vars));
      return fn(...argVals);
    }
  }
}

// ── Public API ──

/** Numeric columns available for expressions (keys from RunRow that hold numbers). */
export const EXPRESSION_COLUMNS = [
  'id', 'alpha', 'reynolds', 'mach', 'ncrit', 'n_panels', 'max_iter',
  'cl', 'cd', 'cm', 'ld', 'iterations', 'residual', 'x_tr_upper', 'x_tr_lower',
] as const;

const COLUMN_SET = new Set<string>(EXPRESSION_COLUMNS);

export interface ValidationResult {
  valid: boolean;
  error?: string;
  errorPos?: number;
}

/**
 * Validate an expression string before evaluation.
 * Returns { valid: true } or { valid: false, error, errorPos }.
 */
export function validateExpression(
  expr: string,
  extraColumns?: string[],
): ValidationResult {
  if (!expr.trim()) {
    return { valid: false, error: 'Expression is empty' };
  }
  try {
    const cols = new Set([...COLUMN_SET, ...(extraColumns ?? [])]);
    const tokens = tokenize(expr);
    const parser = new Parser(tokens, cols);
    parser.parse();
    return { valid: true };
  } catch (err) {
    if (err instanceof ExpressionError) {
      return { valid: false, error: err.message, errorPos: err.pos };
    }
    return { valid: false, error: String(err) };
  }
}

/**
 * Compile an expression into a fast evaluator function.
 * Throws if the expression is invalid.
 */
export function compileExpression(
  expr: string,
  extraColumns?: string[],
): (row: RunRow) => number | null {
  const cols = new Set([...COLUMN_SET, ...(extraColumns ?? [])]);
  const tokens = tokenize(expr);
  const parser = new Parser(tokens, cols);
  const ast = parser.parse();

  return (row: RunRow): number | null => {
    const vars: Record<string, number> = {};
    for (const col of EXPRESSION_COLUMNS) {
      const v = row[col];
      if (typeof v === 'number') vars[col] = v;
      else vars[col] = NaN;
    }
    const result = evalAST(ast, vars);
    return Number.isFinite(result) ? result : null;
  };
}

/**
 * Evaluate an expression for a single row (convenience wrapper).
 * Returns null if the expression is invalid or the result is non-finite.
 */
export function evaluateExpression(expr: string, row: RunRow): number | null {
  try {
    const fn = compileExpression(expr);
    return fn(row);
  } catch {
    return null;
  }
}

/** List of available math function names for UI display. */
export const AVAILABLE_FUNCTIONS = [...FUNCTION_NAMES].sort();
