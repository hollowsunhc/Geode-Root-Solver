/** @file Connect the educational interface to the bounded, worker-based Geode solver. */
import { LIMITS, verifyRequest } from './policy.js';
import { SolverClient } from './client.js';
import { PolynomialPlot, shortRootLabel } from './plot.js';

const examples = {
  crossings: {
    coefficients: '6, -5, -2, 1',
    view: { xMin: -4, xMax: 4, yMin: -12, yMax: 12 },
    explanation: 'A cubic can cross the axis three times. Follow this curve from left to right, then solve to locate each crossing.',
  },
  repeated: {
    coefficients: '2, -3, 0, 1',
    view: { xMin: -3, xMax: 3, yMin: -5, yMax: 9 },
    explanation: 'Watch the curve touch zero and turn around. That is a repeated root: a crossing-only search would miss it. Geode reports each distinct root once.',
  },
  'no-real': {
    coefficients: '1, 0, 1',
    view: { xMin: -4, xMax: 4, yMin: -2, yMax: 10 },
    explanation: 'This parabola stays above zero everywhere. A polynomial can have complex roots without having any real roots for this lab to report.',
  },
  nearby: {
    coefficients: '1.0009765625, -2.0009765625, 1',
    view: { xMin: 0.9995, xMax: 1.0015, yMin: -0.0000004, yMax: 0.0000014 },
    explanation: 'These two roots are only 1/1024 apart. The graph starts zoomed in: zoom out to see how two crossings can seem to become a single touch.',
  },
};

const form = document.querySelector('#solver-form');
const coefficientInput = document.querySelector('#coefficients');
const bitsInput = document.querySelector('#bits');
const digitsInput = document.querySelector('#digits');
const solveButton = document.querySelector('#solve-button');
const cancelButton = document.querySelector('#cancel-button');
const formError = document.querySelector('#form-error');
const status = document.querySelector('#solver-status');
const results = document.querySelector('#results');
const resultHeading = document.querySelector('#result-heading');
const resultSummary = document.querySelector('#result-summary');
const resultBadge = document.querySelector('#result-badge');
const rootList = document.querySelector('#root-list');
const resultNote = document.querySelector('#result-note');
const equation = document.querySelector('#equation');
const client = new SolverClient();
const plot = new PolynomialPlot(document.querySelector('#plot'));
let selectedExample = 'crossings';
let revision = 0;

/** @brief Create a text-only HTML element, keeping user and solver strings out of markup. */
function element(tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
}

/** @brief Read an unrounded coefficient request from the form. */
function currentRequest() {
  return { coefficients: coefficientInput.value, bits: Number(bitsInput.value), digits: Number(digitsInput.value) };
}

/** @brief Set the single live status announcement for solving and input changes. */
function setStatus(message, state = 'ready') {
  status.textContent = message;
  status.dataset.state = state;
}

/** @brief Expose cancellation and announce computation without locking the editable inputs. */
function setBusy(busy) {
  solveButton.disabled = busy;
  solveButton.firstElementChild.textContent = busy ? 'Finding roots…' : 'Find real roots';
  cancelButton.hidden = !busy;
  results.setAttribute('aria-busy', String(busy));
}

/** @brief Clear obsolete output so it can never describe a newly edited polynomial. */
function clearResults() {
  rootList.replaceChildren();
  rootList.hidden = true;
  resultHeading.textContent = 'Find the meeting points';
  resultSummary.textContent = 'Find every distinct real root: the x-values where this polynomial equals zero.';
  resultBadge.textContent = 'READY TO SOLVE';
  resultBadge.classList.remove('is-solved');
  resultNote.hidden = true;
  document.querySelector('#root-legend').lastChild.textContent = ' Roots appear after solving';
  plot.setRoots([]);
}

/** @brief Show the polynomial in descending powers using its original decimal tokens. */
function renderEquation(coefficients) {
  equation.replaceChildren();
  equation.classList.toggle('is-incomplete', !coefficients);
  if (!coefficients) {
    equation.textContent = 'Your polynomial will appear here';
    return;
  }
  equation.append(element('span', 'function-name', 'f(x) = '));
  let first = true;
  for (let power = coefficients.length - 1; power >= 0; power -= 1) {
    const token = coefficients[power];
    if (!/[1-9]/.test(token.split(/[eE]/)[0])) continue;
    const negative = token.startsWith('-');
    const magnitude = token.replace(/^[+-]/, '');
    equation.append(document.createTextNode(first ? (negative ? '−' : '') : (negative ? ' − ' : ' + ')));
    const unit = /^0*1(?:\.0*)?(?:[eE][+-]?0+)?$/.test(magnitude);
    if (power === 0 || !unit) equation.append(document.createTextNode(magnitude));
    if (power > 0) equation.append(element('var', '', 'x'));
    if (power > 1) equation.append(element('sup', '', String(power)));
    first = false;
  }
}

/** @brief Refresh only the lightweight preview; computation always requires an explicit solve. */
function updatePreview() {
  const bits = Number(bitsInput.value);
  const digits = Number(digitsInput.value);
  digitsInput.max = String(LIMITS.digits[bits]);
  const help = document.querySelector('#precision-help');
  help.replaceChildren();
  if (Number.isInteger(digits) && digits >= 1 && digits <= LIMITS.digits[bits]) {
    help.append('Absolute error at most 10', element('sup', '', `−${digits}`), `. Up to ${LIMITS.digits[bits]} digits at ${bits} bits.`);
  } else {
    help.textContent = `Choose 1–${LIMITS.digits[bits]} accuracy digits at ${bits}-bit working precision.`;
  }
  const view = examples[selectedExample]?.view ?? examples.crossings.view;
  try {
    const request = verifyRequest(currentRequest());
    renderEquation(request.coefficients);
    plot.setPolynomial(request.coefficients, view);
  } catch {
    renderEquation(null);
    plot.setPolynomial(null, view);
  }
}

/** @brief Cancel pending work and invalidate old results as soon as any input changes. */
function invalidate(message = 'Inputs changed. Find the roots when you’re ready.') {
  revision += 1;
  client.cancel();
  setBusy(false);
  clearResults();
  formError.hidden = true;
  coefficientInput.removeAttribute('aria-invalid');
  digitsInput.removeAttribute('aria-invalid');
  updatePreview();
  setStatus(message);
}

/** @brief Copy the original MPFR string, with a selectable fallback if clipboard access fails. */
async function copyRoot(root, button, details, code, resultRevision) {
  try {
    await navigator.clipboard.writeText(root);
    if (revision !== resultRevision) return;
    button.textContent = 'Copied';
    setStatus('All printed digits copied.');
  } catch {
    if (revision !== resultRevision) return;
    details.open = true;
    const selection = window.getSelection();
    const range = document.createRange();
    range.selectNodeContents(code);
    selection.removeAllRanges();
    selection.addRange(range);
    code.focus();
    setStatus('The full value is selected. Use your browser’s copy command.');
  }
}

/** @brief Render exact returned strings alongside explicitly approximate short previews. */
function showResults(result, resultRevision) {
  const count = result.roots.length;
  resultHeading.textContent = count ? `${count} distinct real ${count === 1 ? 'root' : 'roots'}` : 'No real roots';
  resultBadge.textContent = `DEGREE ${result.degree} · ${result.bits} BITS`;
  resultBadge.classList.add('is-solved');
  resultSummary.replaceChildren();
  if (count) {
    resultSummary.append('Each returned value is within 10', element('sup', '', `−${result.digits}`), ' of a distinct real root. Short previews below are approximate.');
  } else {
    resultSummary.textContent = 'This polynomial never equals zero on the real number line. Complex roots are outside this lab’s scope.';
  }
  rootList.replaceChildren();
  for (const [index, root] of result.roots.entries()) {
    const row = element('li', 'root-row');
    const label = element('span', 'root-index', 'x');
    label.append(element('sub', '', String(index + 1)));
    const content = element('div', 'root-content');
    const topline = element('div', 'root-topline');
    const preview = element('span', 'root-approx', `≈ ${shortRootLabel(root)}`);
    const copy = element('button', 'copy-root', 'Copy');
    copy.type = 'button';
    copy.setAttribute('aria-label', `Copy all printed digits of root ${index + 1}`);
    const details = element('details', 'root-full');
    const summary = element('summary', '', 'All printed digits');
    const code = element('code', '', root);
    code.tabIndex = 0;
    details.append(summary, code);
    copy.addEventListener('click', () => copyRoot(root, copy, details, code, resultRevision));
    topline.append(preview, copy);
    content.append(topline, details);
    row.append(label, content);
    rootList.append(row);
  }
  rootList.hidden = count === 0;
  resultNote.hidden = false;
  resultNote.textContent = count ? 'The full MPFR values preserve every printed digit. Extra printed digits are not all certified: the error bound guarantees only the requested absolute accuracy, for the stored binary coefficients.' : 'Geode checks the whole real number line, including regions outside the graph’s current view.';
  plot.setRoots(result.roots);
  document.querySelector('#root-legend').lastChild.textContent = count ? ' Computed roots · approximate positions' : ' No real roots';
  setStatus(count ? `${count} distinct real ${count === 1 ? 'root' : 'roots'} found. Results are ready.` : 'Complete. No real roots found.', 'complete');
}

form.addEventListener('submit', async event => {
  event.preventDefault();
  const solveRevision = ++revision;
  client.cancel();
  clearResults();
  formError.hidden = true;
  const request = currentRequest();
  try {
    verifyRequest(request);
  } catch (error) {
    formError.textContent = error.message;
    formError.hidden = false;
    const invalidInput = !Number.isInteger(request.digits) || request.digits < 1 || request.digits > LIMITS.digits[request.bits] ? digitsInput : coefficientInput;
    invalidInput.setAttribute('aria-invalid', 'true');
    invalidInput.focus();
    setStatus('Check the highlighted input, then try again.', 'error');
    return;
  }
  setBusy(true);
  resultBadge.textContent = 'COMPUTING';
  try {
    const result = await client.solve(request, state => {
      if (revision !== solveRevision) return;
      setStatus(state === 'loading' ? 'Loading the local math engine…' : 'Isolating and refining the real roots…', state);
    });
    if (revision !== solveRevision) return;
    showResults(result, solveRevision);
  } catch (error) {
    if (revision !== solveRevision) return;
    if (error.name === 'AbortError') {
      setStatus('Calculation cancelled. You can try again.');
    } else {
      formError.textContent = error.message;
      formError.hidden = false;
      setStatus('The calculation stopped. See the message above.', 'error');
    }
    resultBadge.textContent = 'READY TO TRY AGAIN';
  } finally {
    if (revision === solveRevision) setBusy(false);
  }
});

cancelButton.addEventListener('click', () => invalidate('Calculation cancelled. You can edit the inputs or try again.'));
coefficientInput.addEventListener('input', () => {
  selectedExample = null;
  for (const button of document.querySelectorAll('[data-example]')) {
    button.classList.remove('is-selected');
    button.setAttribute('aria-pressed', 'false');
  }
  document.querySelector('#example-explanation').textContent = 'Your own experiment. Look for crossings and touches, then let the solver check the whole real number line—including beyond this graph.';
  invalidate();
});
bitsInput.addEventListener('change', () => {
  if (Number(digitsInput.value) > LIMITS.digits[Number(bitsInput.value)]) digitsInput.value = String(LIMITS.digits[Number(bitsInput.value)]);
  invalidate();
});
digitsInput.addEventListener('input', () => invalidate());
for (const button of document.querySelectorAll('[data-example]')) {
  button.addEventListener('click', () => {
    selectedExample = button.dataset.example;
    coefficientInput.value = examples[selectedExample].coefficients;
    for (const other of document.querySelectorAll('[data-example]')) {
      other.classList.toggle('is-selected', other === button);
      other.setAttribute('aria-pressed', String(other === button));
    }
    document.querySelector('#example-explanation').textContent = examples[selectedExample].explanation;
    invalidate('Example selected. Find its real roots when you’re ready.');
  });
}
coefficientInput.maxLength = LIMITS.inputLength;
updatePreview();
