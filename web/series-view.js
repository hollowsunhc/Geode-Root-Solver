/** @file Accessible exact-root formulas and certified series convergence, rendered as text and MathML. */
import { shortRootLabel } from './plot.js';

/** @brief Create text-only content without interpreting solver output as markup. */
function element(tag, text) {
  const node = document.createElement(tag);
  if (text !== undefined) node.textContent = text;
  return node;
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
  }

  /** @brief Find a control within this panel. */
  #node(id) {
    return this.#panel.querySelector(`#${id}`);
  }

  /** @brief Remove stale data whenever the polynomial or precision changes. */
  clear() {
    this.#result = null;
    this.#node('series-content').hidden = true;
    this.#node('series-selector').hidden = true;
    this.#node('series-empty').hidden = false;
    this.#node('series-empty').textContent = 'Solve a polynomial to inspect each root’s formula and its convergence.';
    for (const id of ['series-root', 'series-parameters', 'series-steps', 'series-terms', 'series-arguments']) this.#node(id).replaceChildren();
    for (const id of ['series-root-value', 'series-exact-value', 'series-full-sum', 'series-full-bound', 'series-sum', 'series-bound', 'series-unavailable', 'series-reduction-note']) this.#node(id).textContent = '';
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

  /** @brief Select the matching exact value or expansion without recomputing roots. */
  #selectRoot() {
    const index = Number(this.#node('series-root').value);
    const entry = this.#result.series[index];
    this.#node('series-root-value').textContent = `Root ${index + 1} ≈ ${shortRootLabel(this.#result.roots[index])}`;
    const reduction = this.#node('series-reduction-note');
    reduction.hidden = entry.reducedDegree === this.#result.degree;
    reduction.textContent = reduction.hidden ? '' : `The formula uses the polynomial with repeated factors removed (degree ${entry.reducedDegree}). It has the same distinct roots.`;
    this.#node('series-exact').hidden = entry.kind !== 'exact';
    this.#node('series-unavailable').hidden = entry.kind !== 'unavailable';
    this.#node('series-expansion').hidden = entry.kind !== 'series';
    this.#node('series-exact-value').textContent = entry.kind === 'exact' ? entry.value : '';
    this.#node('series-unavailable').textContent = entry.kind === 'unavailable' ? `An expansion is unavailable for this root: ${entry.reason} The decimal root above retains its requested accuracy.` : '';
    for (const id of ['series-parameters', 'series-steps', 'series-terms', 'series-arguments']) this.#node(id).replaceChildren();
    if (entry.kind !== 'series') return;

    this.#node('series-root-index').textContent = String(index + 1);
    const namespace = 'http://www.w3.org/1998/Math/MathML';
    for (let power = 2; power <= entry.reducedDegree; power += 1) {
      if (power > 2) {
        const comma = document.createElementNS(namespace, 'mo');
        comma.textContent = ',';
        this.#node('series-arguments').append(comma);
      }
      const subscript = document.createElementNS(namespace, 'msub');
      const name = document.createElementNS(namespace, 'mi');
      name.textContent = 't';
      const number = document.createElementNS(namespace, 'mn');
      number.textContent = String(power);
      subscript.append(name, number);
      this.#node('series-arguments').append(subscript);
    }
    const values = [['Center a', entry.center], ['Scale δ', entry.scale], ...entry.parameters.map((value, index) => [`t${index + 2}`, value]), ['Certified radius ρ', entry.radius]];
    for (const [name, value] of values) {
      const definition = element('dd');
      const code = element('code', value);
      code.tabIndex = 0;
      definition.append(code);
      this.#node('series-parameters').append(element('dt', name), definition);
    }
    this.#node('series-order').value = '3';
    this.#showOrders();
  }

  /** @brief Show the chosen grouped orders with untouched exact terms and full certified strings. */
  #showOrders() {
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
}
