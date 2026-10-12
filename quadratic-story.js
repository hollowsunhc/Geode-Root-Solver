/** @file A bounded, floating-point illustration of Catalan sums and the two quadratic roots. */
import { curvePath } from './plot.js';

/**
 * @brief Sum at most 81 Catalan terms and map the sum to either quadratic root.
 * @param t Quadratic coefficient, between 1/8 and 1/4 inclusive.
 * @param termCount Number of terms, including C_0 = 1.
 * @param upper Reflect about 1/(2t) to follow the larger root.
 * @return Rounded illustration values; tailBound bounds omitted terms before rounding.
 * @throws RangeError When input exceeds the fixed teaching domain or work budget.
 */
export function quadraticSeries(t, termCount, upper = false) {
  if (!Number.isFinite(t) || t < 1 / 8 || t > 1 / 4 ||
      !Number.isInteger(termCount) || termCount < 1 || termCount > 81 || typeof upper !== 'boolean') {
    throw new RangeError('Use 1/8 ≤ t ≤ 1/4 and 1–81 terms.');
  }
  const lower = 2 / (1 + Math.sqrt(1 - 4 * t));
  const roots = [lower, 1 / t - lower];
  const center = upper ? 1 / t : 0;
  const scale = upper ? -1 : 1;
  const terms = [1];
  let term = 1;
  let sum = 1;
  for (let n = 1; n < termCount; n += 1) {
    term *= (2 * (2 * n - 1) / (n + 1)) * t;
    terms.push(term);
    sum += term;
  }
  // At t=1/4 the tail is (2N+1)C_N/4^N; scale it by (4t)^(N+1) for smaller t.
  const tailBound = 4 * t * (2 * termCount - 1) * term;
  return { roots, center, scale, terms, sum, partial: center + scale * sum, tailBound };
}

/** @brief Express the slider's bounded integer ticks as a reduced, exact fraction. */
function parameterFraction(ticks) {
  let divisor = ticks;
  let remainder = 3600;
  while (remainder) [divisor, remainder] = [remainder, divisor % remainder];
  return `${ticks / divisor}/${3600 / divisor}`;
}

/** @brief Attach the teaching controls without changing the polynomial lab's inputs or worker. */
export function initQuadraticStory(panel) {
  const node = id => panel.querySelector(`#quadratic-${id}`);
  const parameter = node('parameter');
  const count = node('terms');
  const root = node('root');

  /** @brief Update the small illustration using only bounded arithmetic and text nodes. */
  function render() {
    const ticks = Number(parameter.value);
    const t = ticks / 3600;
    const termCount = Number(count.value);
    const upper = root.value === 'upper';
    const result = quadraticSeries(t, termCount, upper);
    const target = result.roots[upper ? 1 : 0];
    const boundary = ticks === 900;
    const rationalExample = ticks === 800;
    const order = termCount - 1;
    const orderSubscript = String(order).replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[Number(digit)]);
    const label = `${parameterFraction(ticks)} ≈ ${t.toFixed(6)}`;
    node('parameter-label').textContent = label;
    parameter.setAttribute('aria-valuetext', `t = ${label}`);
    const countLabel = `${termCount} ${termCount === 1 ? 'term' : 'terms'} · N = ${order}`;
    node('terms-label').textContent = countLabel;
    count.setAttribute('aria-valuetext', countLabel);
    node('center').textContent = upper ? String(Number(result.center.toPrecision(9))) : '0';
    node('scale').textContent = String(result.scale).replace('-', '−');
    node('target').textContent = target.toFixed(9);
    node('partial').textContent = result.partial.toFixed(9);
    node('tail').textContent = result.tailBound.toPrecision(3);
    root.options[0].textContent = boundary ? 'Double root 2 / from below' : `Smaller root ≈ ${Number(result.roots[0].toPrecision(6))}`;
    root.options[1].textContent = boundary ? 'Double root 2 / from above' : `Larger root ≈ ${Number(result.roots[1].toPrecision(6))}`;
    node('exact').textContent = boundary ? 'Exact value: 2 (a double root).'
      : rationalExample ? `Exact value: ${upper ? '3' : '3/2'}. The infinite series is another exact representation.`
        : 'The radical expression is exact; the decimal above is rounded.';
    node('substitution').textContent = `R${orderSubscript} ≈ ${Number(result.center.toPrecision(9))} ${upper ? '−' : '+'} ${result.sum.toFixed(9)} ≈ ${result.partial.toFixed(9)}`;
    node('state').textContent = boundary
      ? 'The two roots meet at x = 2: one double root. The series still converges, but even 81 terms leave a visible gap. It approaches from below for a = 0, and from above for a = 4.'
      : 'There are two distinct roots. Move t toward 1/4 to bring them together; keep the term count fixed to see the remainder bound grow.';
    node('normalization').textContent = upper
      ? 'Set x = 1/t − y. Substitution gives y = 1 + ty² again, so the same Catalan sum reaches the other root. Only the center and the sign of the scale change.'
      : 'Set x = y, so a = 0 and δ = 1. The Catalan sum starts at 1 and approaches the smaller root. Select the larger root to reuse the sum with a different center and scale.';
    const preview = result.terms.slice(0, 4).map(value => Number(value.toPrecision(6))).join(' + ');
    node('sum').textContent = `S${orderSubscript} ≈ ${preview}${termCount > 4 ? ' + ⋯' : ''} ≈ ${result.sum.toFixed(9)}`;

    node('curve').setAttribute('d', curvePath([1, -1, t], { xMin: 0, xMax: 8, yMin: -1.2, yMax: 1.5 }, 576, 180));
    node('lower-marker').setAttribute('cx', String(32 + 72 * result.roots[0]));
    node('upper-marker').setAttribute('cx', String(32 + 72 * result.roots[1]));
    node('upper-marker').setAttribute('visibility', boundary ? 'hidden' : 'visible');
    node('partial-marker').setAttribute('transform', `translate(${32 + 72 * result.partial} 137)`);
    node('plot').setAttribute('aria-label', `Graph of t x squared minus x plus 1 for t = ${parameterFraction(ticks)}. ${boundary ? 'A double root at 2.' : `Roots near ${result.roots.map(value => value.toFixed(6)).join(' and ')}.`} The partial sum is ${result.partial.toFixed(6)}.`);
  }

  parameter.addEventListener('input', render);
  count.addEventListener('input', render);
  root.addEventListener('change', render);
  node('reset').addEventListener('click', () => {
    parameter.value = '800';
    count.value = '9';
    root.value = 'lower';
    render();
  });
  render();
  node('controls').disabled = false;
}
