/** @file Tests the compiled C++ boundary, independently of browser validation. */
import assert from 'node:assert/strict';
import { before, test } from 'node:test';
import { existsSync, readFileSync } from 'node:fs';

const moduleURL = new URL('../../out/web/solver.mjs', import.meta.url);
let createGeode;
let engine;

before(async () => {
  assert.ok(existsSync(moduleURL), 'Build the WebAssembly solver with npm run build:web.');
  ({ default: createGeode } = await import(moduleURL));
  engine = await createGeode();
});

/** @brief Call the C ABI directly, bypassing all JavaScript input checks. */
function solve(coefficients, bits = 256, digits = 30) {
  return JSON.parse(engine.ccall('geode_solve', 'string', ['string', 'number', 'number'],
    [coefficients, bits, digits]));
}

/** @brief Parse a returned decimal as an exact rational, without JavaScript Number rounding. */
function decimalRational(text) {
  const match = /^(-?)(\d+)(?:\.(\d+))?(?:e([+-]?\d+))?$/.exec(text);
  assert.ok(match, 'The result must be decimal text.');
  const fraction = match[3] || '';
  let numerator = BigInt((match[1] || '') + match[2] + fraction);
  let denominator = 1n;
  const exponent = Number(match[4] || 0) - fraction.length;
  if (exponent >= 0) numerator *= 10n ** BigInt(exponent);
  else denominator = 10n ** BigInt(-exponent);
  return { numerator, denominator };
}

test('decimal output retains the absolute bound for a very large stored root', () => {
  const result = solve('-1e100 1', 128, 1);
  assert.equal(result.error, undefined);
  // Exact integer represented by round-to-nearest MPFR 1e100 at 128-bit precision.
  const stored = 9999999999999999999999999999999999999977097647868758532633238403510597078305024613172747796087832576n;
  const { numerator, denominator } = decimalRational(result.roots[0]);
  const difference = numerator - stored * denominator;
  assert.ok((difference < 0n ? -difference : difference) * 10n <= denominator,
    'The printed decimal must be within 0.1 of the stored polynomial root.');
});

test('a small dyadic root is serialized without rounding its decimal expansion', () => {
  const result = solve('-1 1208925819614629174706176', 128, 30);
  const { numerator, denominator } = decimalRational(result.roots[0]);
  assert.equal(numerator * (2n ** 80n), denominator);
});

test('compiled solver finds ordered distinct real roots', () => {
  const result = solve('6, -5, -2, 1');
  assert.equal(result.error, undefined);
  assert.deepEqual(result.roots.map(Number), [-2, 1, 3]);
  assert.equal(result.degree, 3);
  assert.equal(result.bits, 256);
  assert.equal(result.digits, 30);
});

test('repeated roots occur once; absent real roots and constants return empty lists', () => {
  assert.deepEqual(solve('2,-3,0,1').roots.map(Number), [-2, 1]);
  assert.deepEqual(solve('1 0 1').roots, []);
  assert.deepEqual(solve('7').roots, []);
  assert.deepEqual(solve('0, 0, 1, 0').roots.map(Number), [0]);
  assert.equal(solve('0,0,1,0').degree, 2);
});

test('512-bit sqrt(2) results satisfy a 120-digit absolute error bound', () => {
  const result = solve('-2 0 1', 512, 120);
  assert.equal(result.error, undefined);
  assert.equal(result.roots.length, 2);
  for (const root of result.roots) {
    const { numerator, denominator } = decimalRational(root);
    const scaled = (numerator < 0n ? -numerator : numerator) * 10n ** 120n;
    const target = 2n * (denominator * 10n ** 120n) ** 2n;
    assert.ok((scaled - denominator) ** 2n <= target);
    assert.ok((scaled + denominator) ** 2n >= target);
  }
});

test('C++ rejects malicious, malformed, and excessive inputs without UI checks', () => {
  const rejected = ['', '0 0', '1,,2', ',1', '1,', '1, ,2', '1;alert(1)',
    '<img src=x onerror=alert(1)>', 'NaN 1', 'Infinity 1', '0x1p0 1', '1e101 1',
    '1e-101 1', '1e999999999999999 1', '1e+ 1', '--1 1', '1.2.3 1',
    `${'9'.repeat(65)} 1`, `${'1 '.repeat(18)}`, ' '.repeat(2049)];
  for (const input of rejected) assert.equal(typeof solve(input).error, 'string', input);
  for (const [bits, digits] of [[0, 30], [127, 30], [1024, 30], [512, 121],
    [128, 31], [256, 61], [256, 0], [256, -1]]) {
    assert.equal(typeof solve('-1 1', bits, digits).error, 'string');
  }
  const nullResult = JSON.parse(engine.ccall('geode_solve', 'string',
    ['number', 'number', 'number'], [0, 256, 30]));
  assert.equal(typeof nullResult.error, 'string');
  assert.deepEqual(solve('-1 1').roots.map(Number), [1], 'A rejection must not poison the engine.');
});

test('valid boundary values and degree sixteen remain usable', () => {
  assert.deepEqual(solve('1e100 1e100').roots.map(Number), [-1]);
  assert.deepEqual(solve('1e-100 1e-100').roots.map(Number), [-1]);
  assert.deepEqual(solve('.5 +1.').roots.map(Number), [-0.5]);
  const result = solve(['-1', ...Array(15).fill('0'), '1'].join(' '));
  assert.deepEqual(result.roots.map(Number), [-1, 1]);
  assert.equal(result.degree, 16);
});

test('compiled linear memory cannot grow past 128 MiB', async () => {
  let instance;
  await createGeode({
    instantiateWasm(imports, ready) {
      const module = new WebAssembly.Module(readFileSync(new URL('../../out/web/solver.wasm', import.meta.url)));
      instance = new WebAssembly.Instance(module, imports);
      ready(instance, module);
      return instance.exports;
    },
  });
  const memory = Object.values(instance.exports).find(value => value instanceof WebAssembly.Memory);
  assert.ok(memory);
  const pages = memory.buffer.byteLength / 65536;
  assert.throws(() => memory.grow(2049 - pages), RangeError);
  assert.equal(memory.grow(2048 - pages), pages);
  assert.throws(() => memory.grow(1), RangeError);
});
