/** @file Cancellable browser interface to the isolated solver. */
import { LIMITS, verifyRequest } from './policy.js';

/** @brief Bound numeric text before parsing or rendering any worker output. */
function validNumber(value, pattern) {
  return typeof value === 'string' && value.length <= 4096 && pattern.test(value);
}

/** @brief Accept finite decimal strings without converting their precision to Number. */
function validDecimal(value) {
  return validNumber(value, /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/);
}

/** @brief Accept integer fractions with a strictly positive denominator. */
function validRational(value) {
  return validNumber(value, /^[+-]?\d+(?:\/[1-9]\d*)?$/);
}

/** @brief Compare a positive rational with one using bounded digit strings only. */
function radiusAboveOne(value) {
  if (!validRational(value) || value.startsWith('-')) return false;
  const [numerator, denominator = '1'] = value.replace(/^\+?0*/, '').split('/');
  return numerator.length > denominator.length || (numerator.length === denominator.length && numerator > denominator);
}

/** @brief Validate one exact value or certified expansion before exposing it to the view. */
function validSeries(entry, degree) {
  if (!entry || !Number.isInteger(entry.reducedDegree) || entry.reducedDegree < 1 || entry.reducedDegree > degree) return false;
  if (entry.kind === 'exact') return validRational(entry.value);
  if (entry.kind === 'unavailable') return typeof entry.reason === 'string' && entry.reason.length > 0 && entry.reason.length <= 4096;
  return entry.kind === 'series' && entry.reducedDegree >= 2 && validRational(entry.center) &&
    validRational(entry.scale) && !/^[+-]?0+(?:\/|$)/.test(entry.scale) && radiusAboveOne(entry.radius) &&
    Array.isArray(entry.parameters) && entry.parameters.length === entry.reducedDegree - 1 && entry.parameters.every(validRational) &&
    Array.isArray(entry.steps) && entry.steps.length === 9 && entry.steps.every((step, order) => step &&
      step.order === order && validRational(step.coefficient) && validRational(step.term) &&
      validDecimal(step.sum) && validDecimal(step.bound) && !step.bound.startsWith('-'));
}

/** @brief Own at most one worker and terminate it on every terminal path. */
export class SolverClient {
  #cancel = null;

  /**
   * @brief Solve a bounded request, cancelling any earlier calculation.
   * @param request Coefficient text, working bits, and requested digits.
   * @param onStatus Callback receiving loading and computing status changes.
   * @returns Promise of root strings, exact formula data, and calculation metadata.
   * @throws Error through the returned Promise on invalid input or worker failure.
   */
  async solve(request, onStatus = () => {}) {
    this.cancel();
    const normalized = verifyRequest(request);
    return new Promise((resolve, reject) => {
      let worker;
      let timer;
      let computing = false;
      let settled = false;
      const finish = (error, result) => {
        if (settled) return;
        settled = true;
        clearTimeout(timer);
        if (worker) {
          worker.onmessage = worker.onerror = worker.onmessageerror = null;
          worker.terminate();
        }
        if (this.#cancel === cancel) this.#cancel = null;
        if (error) reject(error);
        else resolve(result);
      };
      const cancel = () => finish(new DOMException('Calculation cancelled.', 'AbortError'));
      this.#cancel = cancel;

      try {
        worker = new Worker(new URL('./solver.worker.js', import.meta.url), { type: 'module' });
        timer = setTimeout(() => finish(new Error('The solver did not load within 20 seconds. Please try again.')), LIMITS.loadTimeoutMs);
        worker.onmessage = ({ data }) => {
          if (settled) return;
          try {
            if (data?.type === 'ready' && !computing) {
              computing = true;
              clearTimeout(timer);
              timer = setTimeout(() => finish(new Error('The calculation exceeded 5 seconds. Try a lower degree or lower precision.')), LIMITS.computeTimeoutMs);
              onStatus('computing');
              if (!settled) worker.postMessage({ coefficients: normalized.coefficients.join(' '), bits: normalized.bits, digits: normalized.digits });
            } else if (data?.type === 'result' && computing) {
              const result = data.result;
              if (!result || !Array.isArray(result.roots) || result.roots.length > normalized.degree ||
                  !result.roots.every(validDecimal) || result.degree !== normalized.degree ||
                  !Array.isArray(result.series) || result.series.length !== result.roots.length ||
                  !result.series.every(entry => validSeries(entry, normalized.degree)) ||
                  result.bits !== normalized.bits || result.digits !== normalized.digits) {
                throw new Error('The solver returned an invalid result.');
              }
              finish(null, result);
            } else if (data?.type === 'error' && typeof data.message === 'string') {
              finish(new Error(data.message));
            } else {
              throw new Error('The solver returned an unexpected response.');
            }
          } catch (error) {
            finish(error);
          }
        };
        worker.onerror = event => {
          event.preventDefault();
          finish(new Error('The solver could not load or stopped unexpectedly. Please try again.'));
        };
        worker.onmessageerror = () => finish(new Error('The solver response could not be read.'));
        onStatus('loading');
      } catch (error) {
        finish(error);
      }
    });
  }

  /** @brief Terminate the current worker and reject its Promise with AbortError. */
  cancel() {
    this.#cancel?.();
  }
}
