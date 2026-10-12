/** @file Tests for approximate graph arithmetic, independent of the precise solver. */
import assert from 'node:assert/strict';
import { spawnSync } from 'node:child_process';
import test from 'node:test';
import { axisTicks, curvePath, evaluatePolynomial, shortRootLabel, zoomView } from '../../web/plot.js';

test('evaluates ascending coefficients, including repeated and nearby roots', () => {
  assert.equal(evaluatePolynomial([6, -5, -2, 1], -2), 0);
  assert.equal(evaluatePolynomial([6, -5, -2, 1], 1), 0);
  assert.equal(evaluatePolynomial([6, -5, -2, 1], 3), 0);
  assert.equal(evaluatePolynomial([2, -3, 0, 1], 1), 0);
  assert.equal(evaluatePolynomial([1.0009765625, -2.0009765625, 1], 1.0009765625), 0);
  assert.equal(evaluatePolynomial([1, 0, 1], 0), 1);
});

test('axis ticks remain readable for normal and tightly zoomed views', () => {
  const ordinary = axisTicks(-4, 4);
  assert.ok(ordinary.some(({ value }) => value === 0));
  assert.ok(ordinary.every(({ value }) => value >= -4 && value <= 4));
  const nearby = axisTicks(0.9995, 1.0015);
  assert.ok(nearby.length >= 3);
  assert.equal(new Set(nearby.map(({ label }) => label)).size, nearby.length);
  assert.ok(axisTicks(-2e-7, 8e-7).every(({ label }) => !label.includes('NaN')));
});

test('zoom preserves an anchor and leaves the original bounds untouched', () => {
  const view = { xMin: -4, xMax: 4, yMin: -10, yMax: 10 };
  const zoomed = zoomView(view, 0.5, 2, 5);
  assert.deepEqual(zoomed, { xMin: -1, xMax: 3, yMin: -2.5, yMax: 7.5 });
  assert.deepEqual(view, { xMin: -4, xMax: 4, yMin: -10, yMax: 10 });
});

test('zoom refuses unusable floating-point scales', () => {
  const view = { xMin: -1, xMax: 1, yMin: -1, yMax: 1 };
  assert.deepEqual(zoomView(view, 1e-20), view);
  assert.deepEqual(zoomView(view, 1e20), view);
});

test('curve coordinates are bounded and nonfinite samples never enter SVG', () => {
  const view = { xMin: -4, xMax: 4, yMin: -10, yMax: 10 };
  const ordinary = curvePath([6, -5, -2, 1], view, 600, 300);
  assert.ok(ordinary.startsWith('M'));
  assert.ok(ordinary.includes('L'));
  const extreme = curvePath([1e308, 1e308, 1e308], view, 600, 300);
  assert.ok(!/NaN|Infinity/.test(extreme));
  assert.ok(!/NaN|Infinity/.test(curvePath([1], view, 0, 0)));
});

test('root previews shorten decimal text without rounding it through Number', () => {
  assert.equal(shortRootLabel('-0.200000000000000000e1'), '−2');
  assert.equal(shortRootLabel('0.100000000000000000000000000000000001e1'), '1.00000000000000…');
  assert.equal(shortRootLabel('0.1414213562373095048801688724209698e1'), '1.41421356237309…');
  assert.equal(shortRootLabel('0.100000000000000000e-49'), '1 × 10^−50');
  assert.equal(shortRootLabel('0'), '0');
});

test('very large and degenerate axis ranges cannot stall the page', () => {
  const moduleUrl = new URL('../../web/plot.js', import.meta.url).href;
  const script = `import assert from 'node:assert/strict';
    import { axisTicks } from ${JSON.stringify(moduleUrl)};
    for (const [min, max] of [[1e16, 1e16 + 2], [1, 1], [0, Infinity], [-Infinity, 0], [NaN, 1]]) {
      const ticks = axisTicks(min, max);
      assert.ok(ticks.length <= 16);
      assert.ok(ticks.every(tick => Number.isFinite(tick.value)));
    }`;
  const run = spawnSync(process.execPath, ['--input-type=module', '-e', script], { timeout: 1500, encoding: 'utf8' });
  assert.equal(run.status, 0, run.error?.message ?? run.stderr);
});

test('zoom rejects bounds that collapse after floating-point rounding', () => {
  const view = { xMin: 1e16, xMax: 1e16 + 2, yMin: 0, yMax: 1 };
  assert.deepEqual(zoomView(view, 0.5), view);
});
