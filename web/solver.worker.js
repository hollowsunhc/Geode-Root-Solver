/** @file One-request worker wrapping the compiled MPFR solver. */
import createGeode from './solver.mjs';
import { verifyRequest } from './policy.js';

try {
  const solver = await createGeode();
  self.onmessage = ({ data }) => {
    self.onmessage = null;
    try {
      const request = verifyRequest(data);
      const result = JSON.parse(solver.ccall('geode_solve', 'string', ['string', 'number', 'number'],
        [request.coefficients.join(' '), request.bits, request.digits]));
      if (typeof result.error === 'string') throw new Error(result.error);
      self.postMessage({ type: 'result', result });
    } catch (error) {
      self.postMessage({ type: 'error', message: error instanceof Error ? error.message : 'The calculation failed.' });
    } finally {
      self.close();
    }
  };
  self.postMessage({ type: 'ready' });
} catch (error) {
  self.postMessage({ type: 'error', message: error instanceof Error ? error.message : 'The solver could not load.' });
  self.close();
}
