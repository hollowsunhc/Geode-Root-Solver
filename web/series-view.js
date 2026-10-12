/** @file Accessible exact-root formulas and certified series convergence, rendered as text and MathML. */
import { shortRootLabel } from './plot.js';

/** @brief Create text-only content without interpreting solver output as markup. */
function element(tag, text) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  return node;
}

/** @brief Build MathML with text nodes only, including worker-provided exact numbers. */
function math(tag, ...children) {
  const node = document.createElementNS('http://www.w3.org/1998/Math/MathML', tag);
  if (tag === 'math') node.setAttribute('display', 'block');
  if (tag === 'mi' && children[0] === 'δ') node.setAttribute('mathvariant', 'normal');
  node.append(...children);
  return node;
}

/** @brief Label indexed quantities consistently in controls and plain-text calculations. */
function indexed(name, index) {
  return name + String(index).replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[Number(digit)]);
}

/** @brief Keep fractions recognizable while their full strings remain in exact-value details. */
function compact(value) {
  return value.split('/').map(part => part.length > 14 ? `${part.slice(0, 7)}…${part.slice(-3)}` : part).join('/').replace('-', '−');
}

/** @brief Render a variable with a numeric or symbolic subscript. */
function subscript(name, index) {
  return math('msub', math('mi', name), math(typeof index === 'number' ? 'mn' : 'mi', String(index)));
}

/** @brief Render a compact exact rational without converting it to floating point. */
function rational(value) {
  const [numerator, denominator] = compact(value).split('/');
  return denominator ? math('mfrac', math('mn', numerator), math('mn', denominator)) : math('mn', numerator);
}

/** @brief Insert an arithmetic operator between MathML expressions. */
function joinMath(nodes, operator) {
  return nodes.flatMap((node, index) => index ? [math('mo', operator), node] : [node]);
}

/** @brief Own the selected root and partial-sum controls for a validated solver result. */
export class SeriesView {
  #panel;
  #result = null;

  /** @brief Attach the existing form controls once; all formula arithmetic stays in the worker. */
  constructor(panel) {
    this.#panel = panel;
    this.#node('series-root').addEventListener('change', () => this.#selectRoot());
    this.#node('series-order').addEventListener('input', () => this.#showOrders());
    this.#node('series-product').addEventListener('change', () => this.#showProduct());
  }

  /** @brief Find a control within this panel. */
  #node(id) {
    return this.#panel.querySelector(`#${id}`);
  }

  /** @brief Clear every derived value before displaying a new root or an unavailable report. */
  #clearExpansion() {
    for (const id of [
      'series-root-index', 'series-order-label', 'series-parameters', 'series-parameter-preview', 'series-steps', 'series-terms', 'series-arguments',
      'series-full-sum', 'series-full-bound', 'series-sum', 'series-bound', 'series-product',
      'series-product-count', 'series-product-cap', 'series-exponents', 'series-factorial',
      'series-factorial-ratio', 'series-coefficient', 'series-product-symbolic', 'series-substitution',
      'series-product-value', 'series-numerator-exact', 'series-denominator-exact', 'series-coefficient-exact',
      'series-substitution-exact', 'series-product-value-exact', 'series-group-note', 'series-group-formula',
      'series-group-value', 'series-group-exact', 'series-last-order', 'series-partial-formula', 'series-partial-substitution',
    ]) this.#node(id).replaceChildren();
  }

  /** @brief Remove stale data whenever the polynomial or precision changes. */
  clear() {
    this.#result = null;
    this.#node('series-content').hidden = true;
    this.#node('series-selector').hidden = true;
    this.#node('series-empty').hidden = false;
    this.#node('series-empty').textContent = 'Solve a polynomial to inspect each root’s formula and its convergence.';
    this.#clearExpansion();
    this.#node('series-root').replaceChildren();
    this.#node('series-order').value = '3';
    for (const id of ['series-root-value', 'series-exact-value', 'series-unavailable', 'series-reduction-note']) this.#node(id).textContent = '';
    for (const details of this.#panel.querySelectorAll('details')) details.open = false;
  }

  /** @brief Present the first root or an explicit no-real-root state. */
  show(result) {
    this.clear();
    this.#result = result;
    if (!result.roots.length) {
      this.#node('series-empty').textContent = 'No real roots: there is no real-root expansion to inspect for this polynomial.';
      return;
    }
    for (const [index, root] of result.roots.entries()) {
      const option = element('option', `Root ${index + 1} ≈ ${shortRootLabel(root)}`);
      option.value = String(index);
      this.#node('series-root').append(option);
    }
    this.#node('series-empty').hidden = true;
    this.#node('series-selector').hidden = false;
    this.#node('series-content').hidden = false;
    this.#selectRoot();
  }

  /** @brief Display an expansion alongside any exact value, retaining the inspected product on root changes. */
  #selectRoot() {
    const index = Number(this.#node('series-root').value);
    const entry = this.#result.series[index];
    const productIndex = Number(this.#node('series-product').value);
    this.#clearExpansion();
    this.#node('series-root-value').textContent = `Root ${index + 1} ≈ ${shortRootLabel(this.#result.roots[index])}`;
    const reduction = this.#node('series-reduction-note');
    reduction.hidden = entry.reducedDegree === this.#result.degree;
    reduction.textContent = reduction.hidden ? '' : `The formula uses the polynomial with repeated factors removed (degree ${entry.reducedDegree}). It has the same distinct roots.`;
    this.#node('series-exact').hidden = entry.exactValue === undefined;
    this.#node('series-unavailable').hidden = entry.kind !== 'unavailable';
    this.#node('series-expansion').hidden = entry.kind !== 'series';
    this.#node('series-exact-value').textContent = entry.exactValue ?? '';
    this.#node('series-unavailable').textContent = entry.kind === 'unavailable' ? `An expansion is unavailable for this root: ${entry.reason} The decimal root above retains its requested accuracy.` : '';
    if (entry.kind !== 'series') return;

    this.#node('series-root-index').textContent = String(index + 1);
    this.#node('series-linear-note').hidden = entry.reducedDegree !== 1;
    if (entry.reducedDegree === 1) this.#node('series-order').value = '1';
    if (entry.parameters.length) {
      this.#node('series-arguments').append(math('mo', '('), ...joinMath(entry.parameters.map((_, index) => subscript('t', index + 2)), ','), math('mo', ')'));
    }
    const values = [['Center a', 'a', entry.center], ['Scale δ', 'δ', entry.scale], ...entry.parameters.map((value, index) => [indexed('t', index + 2), indexed('t', index + 2), value]), ['Certified radius ρ', 'ρ', entry.radius]];
    for (const [name, symbol, value] of values) {
      const definition = element('dd');
      const code = element('code', value);
      code.tabIndex = 0;
      definition.append(code);
      this.#node('series-parameters').append(element('dt', name), definition);
      if (symbol !== 'ρ') {
        const preview = element('div');
        preview.append(element('dt', symbol), element('dd', compact(value)));
        this.#node('series-parameter-preview').append(preview);
      }
    }
    this.#showOrders(productIndex);
  }

  /** @brief Show the chosen grouped orders with untouched exact terms and full certified strings. */
  #showOrders(productIndex = 0) {
    if (!this.#result) return;
    const entry = this.#result.series[Number(this.#node('series-root').value)];
    if (entry?.kind !== 'series') return;
    const count = Number(this.#node('series-order').value);
    const selected = entry.steps[count - 1];
    const label = `${count} ${count === 1 ? 'order' : 'orders'} · through j = ${selected.order}`;
    this.#node('series-order-label').textContent = label;
    this.#node('series-order').setAttribute('aria-valuetext', label);
    this.#node('series-sum').textContent = `≈ ${shortRootLabel(selected.sum)}`;
    this.#node('series-bound').textContent = shortRootLabel(selected.bound);
    this.#node('series-full-sum').textContent = selected.sum;
    this.#node('series-full-bound').textContent = selected.bound;
    this.#node('series-last-order').textContent = String(selected.order);
    this.#node('series-partial-formula').replaceChildren(math('math', subscript('R', selected.order), math('mo', '='), math('mi', 'a'), math('mo', '+'), math('mi', 'δ'), math('mo', '('), ...joinMath(entry.steps.slice(0, count).map(step => subscript('s', step.order)), '+'), math('mo', ')')));
    this.#node('series-partial-substitution').textContent = `${indexed('R', selected.order)} = (${compact(entry.center)}) + (${compact(entry.scale)}) × (${entry.steps.slice(0, count).map(step => `(${compact(step.coefficient)})`).join(' + ')})`;
    const productCount = selected.productCount;
    this.#node('series-product-count').textContent = productCount ? `Order j = ${selected.order} has ${productCount} ${productCount === 1 ? 'product' : 'products'}. Choose its exponents to follow the calculation.` : `Order j = ${selected.order} has no products. Its grouped value is 0.`;
    this.#node('series-product-cap').hidden = productCount <= selected.products.length;
    this.#node('series-product-cap').textContent = productCount > selected.products.length ? `Showing the first ${selected.products.length} of ${productCount} products. ${indexed('s', selected.order)} and every partial sum include all ${productCount}.` : '';
    this.#node('series-product-picker').hidden = !productCount;
    this.#node('series-product-calculation').hidden = !productCount;
    this.#node('series-product').replaceChildren();
    for (const [index, product] of selected.products.entries()) {
      const label = product.exponents.length ? product.exponents.map((exponent, index) => `${indexed('m', index + 2)} = ${exponent}`).join(', ') : 'Empty tuple (no parameters)';
      const option = element('option', label);
      option.value = String(index);
      this.#node('series-product').append(option);
    }
    this.#node('series-product').value = String(Math.min(productIndex, Math.max(0, selected.products.length - 1)));
    this.#showProduct();
    this.#node('series-group-note').textContent = productCount ? `Add all ${productCount} ${productCount === 1 ? 'product value' : 'product values'} to get ${indexed('s', selected.order)}.${productCount > selected.products.length ? ' The dots below include every remaining product.' : ''}` : 'There are no products at this order. The linear expansion has only s₀ = 1.';
    const groupValues = selected.products.map(product => math('mrow', math('mo', '('), rational(product.value), math('mo', ')')));
    if (productCount > selected.products.length) groupValues.push(math('mo', '⋯'));
    this.#node('series-group-formula').replaceChildren(math('math', subscript('s', selected.order), math('mo', '='), ...joinMath(groupValues.length ? groupValues : [math('mn', '0')], '+')));
    this.#node('series-group-value').textContent = `${indexed('s', selected.order)} = ${compact(selected.coefficient)}`;
    this.#node('series-group-exact').textContent = selected.coefficient;
    this.#node('series-steps').replaceChildren();
    this.#node('series-terms').replaceChildren();
    for (const step of entry.steps.slice(0, count)) {
      const row = element('tr');
      const order = element('th', String(step.order));
      order.scope = 'row';
      row.append(order, element('td', `≈ ${shortRootLabel(step.sum).replace(/(\.\d{8})\d+…?/, '$1…')}`), element('td', `≈ ${shortRootLabel(step.bound)}`));
      this.#node('series-steps').append(row);
      const terms = element('tr');
      const termOrder = element('th', String(step.order));
      termOrder.scope = 'row';
      terms.append(termOrder);
      for (const value of [step.coefficient, step.term]) {
        const cell = element('td');
        cell.append(element('code', value));
        terms.append(cell);
      }
      this.#node('series-terms').append(terms);
    }
  }

  /** @brief Explain one certified product using its exponents and the selected root’s exact parameters. */
  #showProduct() {
    if (!this.#result) return;
    const entry = this.#result.series[Number(this.#node('series-root').value)];
    if (entry?.kind !== 'series') return;
    const step = entry.steps[Number(this.#node('series-order').value) - 1];
    const product = step.products[Number(this.#node('series-product').value)];
    for (const id of ['series-exponents', 'series-factorial', 'series-factorial-ratio', 'series-coefficient', 'series-product-symbolic', 'series-substitution', 'series-product-value', 'series-numerator-exact', 'series-denominator-exact', 'series-coefficient-exact', 'series-substitution-exact', 'series-product-value-exact']) this.#node(id).replaceChildren();
    if (!product) return;
    this.#node('series-exponents').textContent = product.exponents.length ? `${product.exponents.map((value, index) => `${indexed('m', index + 2)} = ${value}`).join(', ')} → j = ${product.exponents.join(' + ')} = ${step.order}` : 'No exponents: the empty product is 1 at order j = 0.';
    const powers = product.exponents.map((exponent, index) => ({ exponent, power: index + 2, value: entry.parameters[index] })).filter(item => item.exponent > 0);
    const numeratorTerms = powers.map(({ exponent, power }) => math('mrow', math('mn', String(power)), math('mo', '×'), math('mn', String(exponent))));
    const denominatorTerms = [math('mn', '1'), ...powers.map(({ exponent, power }) => math('mrow', math('mn', String(power - 1)), math('mo', '×'), math('mn', String(exponent))))];
    const factorials = () => powers.flatMap(({ exponent }) => [math('mo', '×'), math('mn', String(exponent)), math('mo', '!')]);
    const numeratorOrder = powers.reduce((sum, { exponent, power }) => sum + exponent * power, 0);
    const denominatorOrder = 1 + powers.reduce((sum, { exponent, power }) => sum + exponent * (power - 1), 0);
    this.#node('series-factorial').replaceChildren(math('math', subscript('C', 'm'), math('mo', '='), math('mfrac',
      math('mrow', math('mo', '('), ...joinMath(numeratorTerms.length ? numeratorTerms : [math('mn', '0')], '+'), math('mo', ')'), math('mo', '!')),
      math('mrow', math('mo', '('), ...joinMath(denominatorTerms, '+'), math('mo', ')'), math('mo', '!'), ...factorials())),
    math('mo', '='), math('mfrac', math('mrow', math('mn', String(numeratorOrder)), math('mo', '!')), math('mrow', math('mn', String(denominatorOrder)), math('mo', '!'), ...factorials()))));
    this.#node('series-factorial-ratio').textContent = `${compact(product.numerator)} / ${compact(product.denominator)}`;
    this.#node('series-coefficient').textContent = compact(product.coefficient);
    this.#node('series-product-symbolic').replaceChildren(math('math', subscript('C', 'm'), math('mo', '×'), ...joinMath(powers.length ? powers.map(({ exponent, power }) => math('msup', subscript('t', power), math('mn', String(exponent)))) : [math('mn', '1')], '×')));
    this.#node('series-substitution').replaceChildren(math('math', rational(product.coefficient), math('mo', '×'), ...joinMath(powers.length ? powers.map(({ exponent, value }) => math('msup', math('mrow', math('mo', '('), rational(value), math('mo', ')')), math('mn', String(exponent)))) : [math('mn', '1')], '×')));
    this.#node('series-product-value').textContent = compact(product.value);
    this.#node('series-numerator-exact').textContent = product.numerator;
    this.#node('series-denominator-exact').textContent = product.denominator;
    this.#node('series-coefficient-exact').textContent = product.coefficient;
    this.#node('series-substitution-exact').textContent = `${product.coefficient} × ${powers.length ? powers.map(({ exponent, value }) => `(${value})^${exponent}`).join(' × ') : '1'} = ${product.value}`;
    this.#node('series-product-value-exact').textContent = product.value;
  }
}
