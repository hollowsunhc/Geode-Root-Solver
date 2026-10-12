/** @file Independent checks of exact series terms, branches, and truncation bounds. */
import assert from 'node:assert/strict';
import { before, test } from 'node:test';
import createGeode from '../../out/web/solver.mjs';

let engine;
before(async () => { engine = await createGeode(); });

/** @brief Call the real compiled solver, including its analytical output. */
function solve(coefficients, bits = 256, digits = 30) {
  const result = JSON.parse(engine.ccall('geode_solve', 'string', ['string', 'number', 'number'],
    [coefficients, bits, digits]));
  assert.equal(result.error, undefined);
  assert.equal(result.series?.length, result.roots.length, 'Each root needs its own formula.');
  return result;
}

/** @brief Read fractions or decimal output as exact BigInt pairs for independent checks. */
function rational(text) {
  if (text.includes('/')) return text.split('/').map(BigInt);
  const [, sign, whole, fraction = '', exponent = '0'] = /^(-?)(\d+)(?:\.(\d+))?(?:e([+-]?\d+))?$/.exec(text);
  const power = Number(exponent) - fraction.length;
  const numerator = BigInt(sign + whole + fraction);
  return power < 0 ? [numerator, 10n ** BigInt(-power)] : [numerator * 10n ** BigInt(power), 1n];
}
const add = ([a, b], [c, d]) => [a * d + c * b, b * d];
const neg = ([a, b]) => [-a, b];
const mul = ([a, b], [c, d]) => [a * c, b * d];
const cmp = ([a, b], [c, d]) => a * d < c * b ? -1 : a * d > c * b ? 1 : 0;
const eq = (a, b) => assert.equal(cmp(a, b), 0);

test('quadratic series terms agree with the first Catalan numbers', () => {
  for (const report of solve('-2,0,1').series) {
    assert.equal(report.kind, 'series');
    assert.equal(report.parameters.length, 1);
    assert.ok(cmp(rational(report.radius), [1n, 1n]) > 0);
    const t = rational(report.parameters[0]);
    let power = [1n, 1n];
    for (const [order, catalan] of [1n, 1n, 2n, 5n, 14n, 42n, 132n, 429n, 1430n].entries()) {
      const step = report.steps[order];
      assert.equal(step.order, order);
      eq(rational(step.coefficient), mul([catalan, 1n], power));
      eq(rational(step.term), mul(rational(report.scale), rational(step.coefficient)));
      power = mul(power, t);
    }
  }
});

test('every quadratic partial sum has a valid decreasing absolute error bound', () => {
  const result = solve('-2,0,1', 512, 120);
  for (const [index, report] of result.series.entries()) {
    assert.equal(report.kind, 'series');
    let previous;
    for (const step of report.steps) {
      const sum = rational(step.sum);
      const bound = rational(step.bound);
      assert.ok(cmp(bound, [0n, 1n]) > 0);
      if (previous) assert.ok(cmp(bound, previous) < 0);
      previous = bound;
      const lower = add(sum, neg(bound));
      const upper = add(sum, bound);
      if (index === 0) {
        assert.ok(cmp(upper, [0n, 1n]) < 0);
        assert.ok(cmp(mul(lower, lower), [2n, 1n]) >= 0);
        assert.ok(cmp(mul(upper, upper), [2n, 1n]) <= 0);
      } else {
        assert.ok(cmp(lower, [0n, 1n]) > 0);
        assert.ok(cmp(mul(lower, lower), [2n, 1n]) <= 0);
        assert.ok(cmp(mul(upper, upper), [2n, 1n]) >= 0);
      }
    }
  }
});

test('cubic grouped terms include mixed hyper-Catalan contributions and select each root', () => {
  const result = solve('6,-5,-2,1');
  const expectedRoots = [-2n, 1n, 3n];
  for (const [index, report] of result.series.entries()) {
    assert.equal(report.kind, 'series');
    const [t, u] = report.parameters.map(rational);
    eq(rational(report.steps[1].coefficient), add(t, u));
    eq(rational(report.steps[2].coefficient), add(add(mul([2n, 1n], mul(t, t)), mul([5n, 1n], mul(t, u))), mul([3n, 1n], mul(u, u))));
    for (const step of report.steps) {
      const difference = add(rational(step.sum), [-expectedRoots[index], 1n]);
      assert.ok(cmp([difference[0] < 0n ? -difference[0] : difference[0], difference[1]], rational(step.bound)) <= 0);
    }
  }
});

test('repeated roots use the square-free polynomial, and exact linear roots need no series', () => {
  const repeated = solve('2,-3,0,1');
  for (const report of repeated.series) {
    assert.equal(report.reducedDegree, 2);
    assert.equal(report.kind, 'series');
    assert.equal(report.parameters.length, 1);
  }
  for (const coefficients of ['0.01,-0.2,1', '1,-20,100']) {
    const result = solve(coefficients);
    assert.equal(result.roots.length, 1);
    assert.equal(result.series[0].kind, 'exact');
    eq(rational(result.series[0].value), [1n, 10n]);
  }
  eq(rational(solve('-1,3').series[0].value), [1n, 3n]);
  eq(rational(solve('0,-1,0,1').series[1].value), [0n, 1n]);
  assert.deepEqual(solve('1,0,1').series, []);
});

test('a proposed center that hits an exact rational root retains the correct branch', () => {
  const result = solve('2.0078125,-3.00390625,1'); // (x - 257/256)(x - 2)
  assert.equal(result.series[0].kind, 'exact');
  eq(rational(result.series[0].value), [257n, 256n]);
  assert.equal(result.series[1].kind, 'series');
});

test('a fraction-budget failure preserves numerical roots without claiming a formula', () => {
  const result = solve('1e-100,1e100,1', 512, 120);
  assert.equal(result.roots.length, 2);
  for (const report of result.series) {
    assert.equal(report.kind, 'unavailable');
    assert.match(report.reason, /fraction-size limit/);
    assert.equal(report.steps, undefined);
  }
});

test('formula output stays bounded for allowed high-degree and clustered inputs', () => {
  for (const coefficients of [
    ['-1', ...Array(15).fill('0'), '1'].join(','),
    '1.0009765625,-2.0009765625,1',
    '1e-100,1e100,1',
  ]) {
    const result = solve(coefficients, 512, 1);
    assert.ok(JSON.stringify(result).length < 1000000);
    for (const report of result.series) {
      assert.ok(['series', 'exact', 'unavailable'].includes(report.kind));
      if (report.kind === 'series') assert.equal(report.steps.length, 9);
      if (report.kind === 'unavailable') assert.ok(report.reason.length > 0);
    }
  }
});
