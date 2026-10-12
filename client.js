/** @file Cancellable browser interface to the isolated solver. */
import { LIMITS, verifyRequest } from './policy.js';

/** @brief Own at most one worker and terminate it on every terminal path. */
export class SolverClient {
  #cancel = null;

  /**
   * @brief Solve a bounded request, cancelling any earlier calculation.
   * @param request Coefficient text, working bits, and requested digits.
   * @param onStatus Callback receiving loading and computing status changes.
   * @returns Promise of MPFR root strings and the calculation metadata.
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
                  !result.roots.every(root => typeof root === 'string') || result.degree !== normalized.degree ||
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
