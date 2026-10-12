/** @file Browser workload policy shared by the page and its worker. */

/** @brief Hard bounds for a single browser calculation. */
export const LIMITS = Object.freeze({
  inputLength: 2048,
  coefficientCount: 17,
  coefficientLength: 64,
  exponentMagnitude: 100,
  bits: Object.freeze([128, 256, 512]),
  digits: Object.freeze({ 128: 30, 256: 60, 512: 120 }),
  loadTimeoutMs: 20000,
  computeTimeoutMs: 5000,
});

/**
 * @brief Validate a solver request without rounding its coefficient strings.
 * @param request Object containing coefficient text, working bits, and digits.
 * @returns Bounded decimal tokens, precision, and effective polynomial degree.
 * @throws Error when the input is malformed or exceeds a workload limit.
 */
export function verifyRequest(request) {
  if (!request || typeof request !== 'object' || Array.isArray(request)) {
    throw new Error('Provide coefficients, working precision, and requested digits.');
  }
  const { coefficients: text, bits, digits } = request;
  if (typeof text !== 'string') throw new Error('Coefficients must be decimal text.');
  if (text.length > LIMITS.inputLength) {
    throw new Error(`Use at most ${LIMITS.inputLength} input characters.`);
  }
  if (text.includes('\0')) throw new Error('Coefficient text cannot contain null characters.');
  if (!LIMITS.bits.includes(bits)) throw new Error('Working precision must be 128, 256, or 512 bits.');
  if (!Number.isInteger(digits) || digits < 1 || digits > LIMITS.digits[bits]) {
    throw new Error(`Requested digits must be an integer from 1 to ${LIMITS.digits[bits]}.`);
  }

  const fields = text.split(',');
  if (fields.some(field => !field.trim())) {
    throw new Error('Enter decimal coefficients; comma-separated fields cannot be empty.');
  }
  const coefficients = fields.flatMap(field => field.trim().split(/\s+/u));
  if (coefficients.length > LIMITS.coefficientCount) {
    throw new Error(`Use at most ${LIMITS.coefficientCount} coefficients (degree 16).`);
  }
  for (const coefficient of coefficients) {
    if (coefficient.length > LIMITS.coefficientLength) {
      throw new Error(`Each coefficient may contain at most ${LIMITS.coefficientLength} characters.`);
    }
    const decimal = coefficient.match(/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE]([+-]?\d+))?$/);
    if (!decimal) throw new Error('Use decimal numbers only, without expressions or other symbols.');
    if (decimal[1] !== undefined && Math.abs(Number(decimal[1])) > LIMITS.exponentMagnitude) {
      throw new Error(`Decimal exponents must be between -${LIMITS.exponentMagnitude} and ${LIMITS.exponentMagnitude}.`);
    }
  }
  while (coefficients.length && !/[1-9]/.test(coefficients.at(-1).split(/[eE]/)[0])) {
    coefficients.pop();
  }
  if (!coefficients.length) throw new Error('The zero polynomial has no isolated roots; use a nonzero coefficient.');
  return { coefficients, bits, digits, degree: coefficients.length - 1 };
}
