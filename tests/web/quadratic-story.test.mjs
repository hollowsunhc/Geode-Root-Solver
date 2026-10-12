/** @file Mathematical and workload checks for the interactive quadratic lesson. */
import test from 'node:test';
import assert from 'node:assert/strict';
import { quadraticSeries } from '../../web/quadratic-story.js';

test('the rational example uses Catalan terms and reflects the same sum to the other root', () => {
  const lower = quadraticSeries(2 / 9, 3);
  const upper = quadraticSeries(2 / 9, 3, true);
  assert.ok(Math.abs(lower.roots[0] - 1.5) < 1e-14);
  assert.ok(Math.abs(lower.roots[1] - 3) < 1e-14);
  assert.ok(Math.abs(lower.sum - 107 / 81) < 1e-14);
  assert.deepEqual(lower.terms, [1, 2 / 9, 8 / 81]);
  assert.equal(lower.center, 0);
  assert.equal(lower.scale, 1);
  assert.equal(upper.center, 4.5);
  assert.equal(upper.scale, -1);
  assert.equal(upper.sum, lower.sum);
  assert.ok(Math.abs(upper.partial - 515 / 162) < 1e-14);
  assert.ok(Math.abs(lower.tailBound - 320 / 729) < 1e-14);
  assert.ok(1.5 - lower.partial <= lower.tailBound);
  assert.ok(upper.partial - 3 <= upper.tailBound);
  const longer = quadraticSeries(2 / 9, 81);
  assert.equal(longer.terms.length, 81);
  assert.ok(Math.abs(longer.partial - 1.5) < 1e-6);
  assert.ok(longer.tailBound < lower.tailBound);
});

test('the series still converges when the roots meet, with a finite boundary remainder', () => {
  const lower = quadraticSeries(1 / 4, 3);
  const upper = quadraticSeries(1 / 4, 3, true);
  assert.deepEqual(lower.roots, [2, 2]);
  assert.equal(lower.sum, 11 / 8);
  assert.equal(lower.tailBound, 5 / 8);
  assert.equal(upper.partial, 21 / 8);
  assert.equal(upper.partial - 2, lower.tailBound);
  const first = quadraticSeries(1 / 4, 1);
  assert.equal(first.partial, 1);
  assert.equal(first.tailBound, 1);
  const longer = quadraticSeries(1 / 4, 81);
  assert.ok(longer.partial > lower.partial && longer.partial < 2);
  assert.ok(Math.abs((2 - longer.partial) - longer.tailBound) < 1e-14);
  assert.ok(longer.tailBound > quadraticSeries(2 / 9, 81).tailBound);
});

test('the lesson rejects invalid parameters and unbounded term counts before summing', () => {
  for (const t of [0, -1, 0.124, 0.251, Infinity, NaN, '0.2']) {
    assert.throws(() => quadraticSeries(t, 9), RangeError);
  }
  for (const count of [0, -1, 82, 1.5, Number.MAX_SAFE_INTEGER, Infinity, NaN, '9']) {
    assert.throws(() => quadraticSeries(2 / 9, count), RangeError);
  }
  assert.throws(() => quadraticSeries(2 / 9, 9, 'false'), RangeError);
  assert.ok(Number.isFinite(quadraticSeries(1 / 8, 81).partial));
});
