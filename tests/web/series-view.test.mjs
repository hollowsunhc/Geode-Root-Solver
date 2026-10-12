/** @file Check the result boundary before exact series data reaches the browser view. */
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { SolverClient } from '../../web/client.js';

/** @brief Deliver a worker payload through the real client validation boundary. */
async function receive(result) {
  const originalWorker = globalThis.Worker;
  globalThis.Worker = class {
    constructor() { queueMicrotask(() => this.onmessage?.({ data: { type: 'ready' } })); }
    postMessage() { queueMicrotask(() => this.onmessage?.({ data: { type: 'result', result } })); }
    terminate() {}
  };
  try {
    return await new SolverClient().solve({ coefficients: '0.01 -0.2 1', bits: 128, digits: 20 });
  } finally {
    globalThis.Worker = originalWorker;
  }
}

/** @brief Provide a complete result fixture for structural validation only. */
function resultWith(entry) {
  return { degree: 2, bits: 128, digits: 20, roots: ['0.1'], series: [entry] };
}

/** @brief Supply bounded product metadata without treating mock values as a mathematical proof. */
function seriesFixture(parameters = ['1/36']) {
  return {
    kind: 'series', center: '3/2', scale: '-1/12', parameters, radius: '2', reducedDegree: parameters.length + 1,
    steps: Array.from({ length: 9 }, (_, order) => {
      const productCount = parameters.length || order === 0 ? 1 : 0;
      return { order, coefficient: '1', term: '-1/12', sum: '1.41666666667', bound: '0.01', productCount,
        products: productCount ? [{ exponents: parameters.map(() => order), numerator: '1', denominator: '1', coefficient: '1', value: '1' }] : [] };
    }),
  };
}

test('accepts bounded expansions, including exact and linear roots', async () => {
  for (const entry of [
    { ...seriesFixture([]), exactValue: '1/10' },
    { kind: 'unavailable', reason: 'Series exceeded the display limit.', reducedDegree: 2 },
    seriesFixture(),
  ]) {
    const result = resultWith(entry);
    assert.deepEqual(await receive(result), result);
  }
});

test('rejects missing, mismatched, malformed, or oversized formula payloads', async () => {
  const valid = { ...seriesFixture([]), exactValue: '1/10' };
  const series = seriesFixture();
  const malformed = [
    { ...resultWith(valid), series: undefined },
    { ...resultWith(valid), series: [] },
    { ...resultWith(valid), roots: ['NaN'] },
    { ...resultWith(valid), roots: ['1'.repeat(4097)] },
    ...[
      { ...valid, exactValue: '1/0' },
      { ...valid, exactValue: '<img>' },
      { ...valid, exactValue: '1'.repeat(4097) },
      { kind: 'exact', value: '1/10', reducedDegree: 1 },
      { ...valid, reducedDegree: 3 },
      { kind: 'unavailable', reason: 123, reducedDegree: 2 },
      { ...series, parameters: [] },
      { ...series, parameters: ['1/0'] },
      { ...series, radius: '0' },
      { ...series, radius: '1/2' },
      { ...series, scale: '0' },
      { ...series, steps: series.steps.slice(0, 8) },
      { ...series, steps: series.steps.map(step => ({ ...step, order: 0 })) },
      { ...series, steps: series.steps.map(step => ({ ...step, bound: '-1' })) },
      { ...series, steps: series.steps.map(step => ({ ...step, term: '1/0' })) },
      { ...series, steps: series.steps.map(step => ({ ...step, sum: 'Infinity' })) },
    ].map(resultWith),
  ];
  for (const result of malformed) {
    await assert.rejects(receive(result), /invalid result/i);
  }
});

test('rejects oversized and malformed worked products before rendering', async () => {
  const fixture = seriesFixture();
  const step = fixture.steps[2];
  const product = step.products[0];
  const invalid = [
    { ...step, products: undefined },
    { ...step, products: [] },
    { ...step, products: Array(13).fill(product) },
    { ...step, productCount: 319771 },
    { ...step, productCount: -1 },
    { ...step, productCount: 1.5 },
    ...[
      { ...product, exponents: [1] },
      { ...product, exponents: [2, 0] },
      { ...product, exponents: [2.5] },
      { ...product, exponents: [-2] },
      { ...product, coefficient: '2/3' },
      { ...product, coefficient: '<script>' },
      { ...product, numerator: '9'.repeat(4097) },
      { ...product, denominator: '0' },
      { ...product, value: '1/0' },
    ].map(product => ({ ...step, products: [product] })),
  ];
  for (const entry of invalid) {
    await assert.rejects(receive(resultWith({ ...fixture, steps: fixture.steps.map(step => step.order === 2 ? entry : step) })), /invalid result/i);
  }
});
