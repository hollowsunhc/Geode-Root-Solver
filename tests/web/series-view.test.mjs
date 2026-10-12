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

test('accepts exact, unavailable, and bounded series representations', async () => {
  for (const entry of [
    { kind: 'exact', value: '1/10', reducedDegree: 1 },
    { kind: 'unavailable', reason: 'Series exceeded the display limit.', reducedDegree: 2 },
    {
      kind: 'series', center: '3/2', scale: '-1/12', parameters: ['1/36'], radius: '2', reducedDegree: 2,
      steps: Array.from({ length: 9 }, (_, order) => ({ order, coefficient: '1', term: '-1/12', sum: '1.41666666667', bound: '0.01' })),
    },
  ]) {
    const result = resultWith(entry);
    assert.deepEqual(await receive(result), result);
  }
});

test('rejects missing, mismatched, malformed, or oversized formula payloads', async () => {
  const valid = { kind: 'exact', value: '1/10', reducedDegree: 1 };
  const series = {
    kind: 'series', center: '3/2', scale: '-1/12', parameters: ['1/36'], radius: '2', reducedDegree: 2,
    steps: Array.from({ length: 9 }, (_, order) => ({ order, coefficient: '1', term: '-1/12', sum: '1.41666666667', bound: '0.01' })),
  };
  const malformed = [
    { ...resultWith(valid), series: undefined },
    { ...resultWith(valid), series: [] },
    { ...resultWith(valid), roots: ['NaN'] },
    { ...resultWith(valid), roots: ['1'.repeat(4097)] },
    ...[
      { ...valid, value: '1/0' },
      { ...valid, value: '<img>' },
      { ...valid, value: '1'.repeat(4097) },
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
