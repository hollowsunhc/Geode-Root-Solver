import assert from 'node:assert/strict';
import test from 'node:test';
import { verifyRequest } from '../../web/policy.js';

const request = (coefficients, bits = 128, digits = 20) => ({ coefficients, bits, digits });

test('keeps exact decimal strings while accepting comma and whitespace separators', () => {
  assert.deepEqual(verifyRequest(request('  -1,\t+0.25  .5, 1.e+02\n')), {
    coefficients: ['-1', '+0.25', '.5', '1.e+02'], bits: 128, digits: 20, degree: 3,
  });
  assert.deepEqual(verifyRequest(request('0.12345678901234567890123456789 1')).coefficients,
    ['0.12345678901234567890123456789', '1']);
});

test('removes trailing zero coefficients without losing a nonzero constant', () => {
  assert.deepEqual(verifyRequest(request('-1 1 0.0e10 -0')), {
    coefficients: ['-1', '1'], bits: 128, digits: 20, degree: 1,
  });
  assert.equal(verifyRequest(request('1e-100 0')).degree, 0);
});

test('accepts each precision at its requested digit boundary', () => {
  for (const [bits, digits] of [[128, 30], [256, 60], [512, 120]]) {
    assert.equal(verifyRequest(request('1 1', bits, digits)).digits, digits);
    assert.equal(verifyRequest(request('1 1', bits, 1)).digits, 1);
    assert.throws(() => verifyRequest(request('1 1', bits, digits + 1)), /digits/i);
  }
});

test('accepts exactly the input, coefficient, token, and exponent bounds', () => {
  assert.equal(verifyRequest(request(' '.repeat(2047) + '1')).degree, 0);
  assert.equal(verifyRequest(request(Array(17).fill('1').join(' '))).degree, 16);
  assert.equal(verifyRequest(request('1'.repeat(64))).coefficients[0].length, 64);
  assert.deepEqual(verifyRequest(request('1e-100 1e+100')).coefficients, ['1e-100', '1e+100']);
});

test('rejects oversized input before trailing zero normalization', () => {
  for (const value of [
    ' '.repeat(2048) + '1',
    '1 ' + Array(17).fill('0').join(' '),
    '1'.repeat(65),
    '1e101',
    '1e-101',
    '1e999999999999999999999999999999999999999999999999999999999999',
  ]) {
    assert.throws(() => verifyRequest(request(value)), Error, value);
  }
});

test('rejects empty comma fields and zero polynomials', () => {
  for (const value of ['', ' \n\t ', ',1', '1,', '1,,2', '1, \n,2', '0', '-0 +0.0 0e100']) {
    assert.throws(() => verifyRequest(request(value)), Error, JSON.stringify(value));
  }
});

test('rejects expressions, nondecimal values, injection strings, and embedded nulls', () => {
  for (const value of [
    'NaN', 'Infinity', '-Infinity', '0x10', '0b11', '1_000', '1/2', '2**100000000',
    'x^2-1', '1e', '1e+', '--1', '.', '+', '1.2.3', '1;fetch("https://example.invalid")',
    '<img src=x onerror=alert(1)>', '${globalThis.location}', '1\0 1', '1\u200b2',
  ]) {
    assert.throws(() => verifyRequest(request(value)), Error, JSON.stringify(value));
  }
});

test('requires the documented request and field types without coercion', () => {
  for (const value of [null, undefined, [], '1 1', 4, {},
    { coefficients: ['1', '1'], bits: 128, digits: 20 },
    { coefficients: 1, bits: 128, digits: 20 },
    request('1 1', '128'), request('1 1', 64), request('1 1', 128.5),
    request('1 1', Infinity), request('1 1', NaN), request('1 1', null),
    request('1 1', 128, '20'), request('1 1', 128, 0), request('1 1', 128, -1),
    request('1 1', 128, 1.5), request('1 1', 128, NaN), request('1 1', 128, Infinity),
  ]) {
    assert.throws(() => verifyRequest(value), Error);
  }
});
