/** @file Behavioral checks for browser workload limits and isolation. */
import { test, expect } from '@playwright/test';

const normalRequest = { coefficients: '-1 0 1', bits: 128, digits: 20 };

test.beforeEach(async ({ page }) => {
  await page.goto('./');
  await page.evaluate(() => {
    const NativeWorker = globalThis.Worker;
    globalThis.workerAudit = { created: 0, active: 0, peak: 0, terminated: 0 };
    // Observe actual workers and their termination without replacing their work.
    globalThis.Worker = class extends NativeWorker {
      constructor(...args) {
        super(...args);
        this.counted = true;
        workerAudit.created += 1;
        workerAudit.active += 1;
        workerAudit.peak = Math.max(workerAudit.peak, workerAudit.active);
      }
      terminate() {
        if (this.counted) {
          this.counted = false;
          workerAudit.active -= 1;
          workerAudit.terminated += 1;
        }
        return super.terminate();
      }
    };
  });
});

test('rejects hostile requests before allocating a worker', async ({ page }) => {
  const outcome = await page.evaluate(async () => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    const client = new SolverClient();
    const errors = [];
    for (const coefficients of ['1'.repeat(2049), '1 0\0 1', '1,,2', '0', '1e101',
      Array(18).fill('1').join(' '), '<img src=x onerror=alert(1)>', 'fetch("https://example.invalid")']) {
      try {
        await client.solve({ coefficients, bits: 128, digits: 20 });
        errors.push(null);
      } catch (error) {
        errors.push(error.message);
      }
    }
    return { errors, audit: workerAudit };
  });
  expect(outcome.errors).toHaveLength(8);
  expect(outcome.errors.every(message => typeof message === 'string' && message.length > 0)).toBe(true);
  expect(outcome.audit.created).toBe(0);
});

test('computes precise strings locally and terminates the completed worker', async ({ page }) => {
  const requests = [];
  page.on('request', request => requests.push({ url: request.url(), method: request.method(), body: request.postData() }));
  const outcome = await page.evaluate(async request => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    const statuses = [];
    const result = await new SolverClient().solve(request, status => statuses.push(status));
    return { result, statuses, audit: workerAudit };
  }, normalRequest);
  expect(outcome.statuses).toEqual(['loading', 'computing']);
  expect(outcome.result).toMatchObject({ degree: 2, bits: 128, digits: 20 });
  expect(outcome.result.roots.every(root => typeof root === 'string')).toBe(true);
  expect(outcome.result.roots).toHaveLength(2);
  expect(Number(outcome.result.roots[0])).toBeCloseTo(-1, 10);
  expect(Number(outcome.result.roots[1])).toBeCloseTo(1, 10);
  expect(outcome.audit).toEqual({ created: 1, active: 0, peak: 1, terminated: 1 });
  for (const request of requests) {
    expect(new URL(request.url).origin).toBe(new URL(page.url()).origin);
    expect(new URL(request.url).pathname.split('/').at(-1)).toMatch(/^(client\.js|policy\.js|solver\.worker\.js|solver\.mjs|solver\.wasm)$/);
    expect(new URL(request.url).search).toBe('');
    expect(request.method).toBe('GET');
    expect(request.body).toBeNull();
  }
});

test('worker independently rejects invalid requests that bypass the client', async ({ page }) => {
  const outcomes = await page.evaluate(async () => {
    const results = [];
    for (const request of [
      { coefficients: '1'.repeat(2049), bits: 128, digits: 20 },
      { coefficients: '1\0 1', bits: 128, digits: 20 },
      { coefficients: '1,,2', bits: 128, digits: 20 },
      { coefficients: Array(18).fill('1').join(' '), bits: 128, digits: 20 },
      { coefficients: '0', bits: 128, digits: 20 },
      { coefficients: '1e101 1', bits: 128, digits: 20 },
      { coefficients: '1 1', bits: '128', digits: 20 },
      { coefficients: '1 1', bits: 128, digits: 31 },
    ]) {
      results.push(await new Promise((resolve, reject) => {
        const worker = new Worker(new URL('solver.worker.js', document.baseURI), { type: 'module' });
        const timer = setTimeout(() => { worker.terminate(); reject(new Error('Worker did not answer.')); }, 20000);
        worker.onerror = () => { clearTimeout(timer); worker.terminate(); reject(new Error('Worker failed to load.')); };
        worker.onmessage = ({ data }) => {
          if (data.type === 'ready') worker.postMessage(request);
          else { clearTimeout(timer); worker.terminate(); resolve(data); }
        };
      }));
    }
    return results;
  });
  expect(outcomes).toHaveLength(8);
  for (const outcome of outcomes) {
    expect(outcome.type).toBe('error');
    expect(outcome.message).toEqual(expect.any(String));
    expect(outcome.message.length).toBeGreaterThan(0);
  }
});

test('a worker accepts only one request even when more are already queued', async ({ page }) => {
  const messages = await page.evaluate(request => new Promise((resolve, reject) => {
    const worker = new Worker(new URL('solver.worker.js', document.baseURI), { type: 'module' });
    const results = [];
    const timer = setTimeout(() => { worker.terminate(); reject(new Error('Worker did not answer.')); }, 20000);
    worker.onerror = () => { clearTimeout(timer); worker.terminate(); reject(new Error('Worker failed to load.')); };
    worker.onmessage = ({ data }) => {
      if (data.type === 'ready') {
        worker.postMessage(request);
        worker.postMessage({ coefficients: '1 0 1', bits: 128, digits: 20 });
      } else {
        results.push(data);
        if (results.length === 1) setTimeout(() => { clearTimeout(timer); worker.terminate(); resolve(results); }, 150);
      }
    };
  }), normalRequest);
  expect(messages).toHaveLength(1);
  expect(messages[0].type).toBe('result');
  expect(messages[0].result.roots).toHaveLength(2);
});

test('a burst supersedes old requests with at most one live worker', async ({ page }) => {
  const outcome = await page.evaluate(async request => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    const client = new SolverClient();
    const pending = Array.from({ length: 8 }, () => client.solve(request));
    const results = await Promise.allSettled(pending);
    return { results: results.map(result => ({ status: result.status, name: result.reason?.name, roots: result.value?.roots })), audit: workerAudit };
  }, normalRequest);
  expect(outcome.results.slice(0, -1)).toEqual(Array(7).fill({ status: 'rejected', name: 'AbortError', roots: undefined }));
  expect(outcome.results.at(-1).status).toBe('fulfilled');
  expect(outcome.results.at(-1).roots).toHaveLength(2);
  expect(outcome.audit).toEqual({ created: 8, active: 0, peak: 1, terminated: 8 });
});

test('cancels loading immediately and can solve again', async ({ page }) => {
  const outcome = await page.evaluate(async request => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    const client = new SolverClient();
    const pending = client.solve(request).then(() => 'unexpected success', error => error.name);
    client.cancel();
    client.cancel();
    const cancelled = await pending;
    const result = await client.solve(request);
    return { cancelled, result, audit: workerAudit };
  }, normalRequest);
  expect(outcome.cancelled).toBe('AbortError');
  expect(outcome.result.roots).toHaveLength(2);
  expect(outcome.audit).toEqual({ created: 2, active: 0, peak: 1, terminated: 2 });
});

/** @brief Start a real worker that receives its request, then consumes CPU forever. */
async function startHungWorker(page) {
  await page.route('**/solver.worker.js', route => route.fulfill({
    contentType: 'text/javascript',
    body: `self.onmessage = () => { console.log('security-worker-started'); while (true) {} }; self.postMessage({ type: 'ready' });`,
  }), { times: 1 });
  const started = page.waitForEvent('console', { predicate: message => message.text() === 'security-worker-started' });
  await page.evaluate(async request => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    globalThis.testClient = new SolverClient();
    globalThis.heartbeat = 0;
    globalThis.successes = 0;
    const interval = setInterval(() => heartbeat += 1, 10);
    const start = performance.now();
    globalThis.pendingSolve = testClient.solve(request).then(result => {
      successes += 1;
      return { result };
    }, error => ({ name: error.name, message: error.message })).then(outcome => {
      clearInterval(interval);
      return { ...outcome, elapsed: performance.now() - start };
    });
  }, normalRequest);
  await started;
}

test('cancels a running CPU-bound worker without freezing the page or delivering stale results', async ({ page }) => {
  await startHungWorker(page);
  await expect.poll(() => page.evaluate(() => heartbeat)).toBeGreaterThan(5);
  const outcome = await page.evaluate(async request => {
    testClient.cancel();
    const cancelled = await pendingSolve;
    const result = await testClient.solve(request);
    return { cancelled, result, successes, audit: workerAudit };
  }, normalRequest);
  expect(outcome.cancelled.name).toBe('AbortError');
  expect(outcome.successes).toBe(0);
  expect(outcome.result.roots).toHaveLength(2);
  expect(outcome.audit).toEqual({ created: 2, active: 0, peak: 1, terminated: 2 });
});

test('terminates CPU-bound work at the compute deadline and recovers', async ({ page }) => {
  await startHungWorker(page);
  const outcome = await page.evaluate(async request => {
    const timeout = await pendingSolve;
    const ticks = heartbeat;
    const result = await testClient.solve(request);
    return { timeout, ticks, result, successes, audit: workerAudit };
  }, normalRequest);
  expect(outcome.timeout.message).toMatch(/5 seconds|timed out/i);
  expect(outcome.timeout.elapsed).toBeGreaterThanOrEqual(4900);
  expect(outcome.timeout.elapsed).toBeLessThan(10000);
  expect(outcome.ticks).toBeGreaterThan(20);
  expect(outcome.successes).toBe(0);
  expect(outcome.result.roots).toHaveLength(2);
  expect(outcome.audit).toEqual({ created: 2, active: 0, peak: 1, terminated: 2 });
});

test('terminates a worker that never becomes ready at the loading deadline', async ({ page }) => {
  await page.route('**/solver.worker.js', route => route.fulfill({ contentType: 'text/javascript', body: 'self.onmessage = () => {};' }), { times: 1 });
  const outcome = await page.evaluate(async request => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    const client = new SolverClient();
    const start = performance.now();
    const error = await client.solve(request).then(() => null, error => error.message);
    const elapsed = performance.now() - start;
    const result = await client.solve(request);
    return { error, elapsed, result, audit: workerAudit };
  }, normalRequest);
  expect(outcome.error).toMatch(/20 seconds|timed out/i);
  expect(outcome.elapsed).toBeGreaterThanOrEqual(19900);
  expect(outcome.elapsed).toBeLessThan(26000);
  expect(outcome.result.roots).toHaveLength(2);
  expect(outcome.audit).toEqual({ created: 2, active: 0, peak: 1, terminated: 2 });
});

test('terminates a failed worker load and recovers without an automatic retry', async ({ page }) => {
  await page.route('**/solver.worker.js', route => route.abort('failed'), { times: 1 });
  const outcome = await page.evaluate(async request => {
    const { SolverClient } = await import(new URL('client.js', document.baseURI));
    const client = new SolverClient();
    const error = await client.solve(request).then(() => null, error => error.message);
    const afterFailure = { ...workerAudit };
    const result = await client.solve(request);
    return { error, afterFailure, result, audit: workerAudit };
  }, normalRequest);
  expect(outcome.error).toEqual(expect.any(String));
  expect(outcome.afterFailure).toEqual({ created: 1, active: 0, peak: 1, terminated: 1 });
  expect(outcome.result.roots).toHaveLength(2);
  expect(outcome.audit).toEqual({ created: 2, active: 0, peak: 1, terminated: 2 });
});

test('CSP denies injected scripts, dynamic JavaScript, foreign connections, and blob workers', async ({ page }) => {
  // DevTools evaluation bypasses eval restrictions, so run the probe as a page script.
  const probe = async () => {
    const violations = [];
    document.addEventListener('securitypolicyviolation', event => violations.push(event.effectiveDirective));
    globalThis.injectedScriptRan = false;
    const script = document.createElement('script');
    script.textContent = 'globalThis.injectedScriptRan = true;';
    document.head.append(script);
    let evaluated = false;
    try { evaluated = Function('return true')(); } catch { /* Expected CSP denial. */ }
    let connected = false;
    try { await fetch('https://example.invalid/geode-security-check'); connected = true; } catch { /* Expected CSP denial. */ }
    const blob = URL.createObjectURL(new Blob(['self.postMessage("unsafe");'], { type: 'text/javascript' }));
    const blobWorkerStarted = await new Promise(resolve => {
      try {
        const worker = new Worker(blob);
        worker.onmessage = () => { worker.terminate(); resolve(true); };
        worker.onerror = event => { event.preventDefault(); worker.terminate(); resolve(false); };
      } catch { resolve(false); }
    });
    URL.revokeObjectURL(blob);
    await new Promise(resolve => setTimeout(resolve, 100));
    return { injected: injectedScriptRan, evaluated, connected, blobWorkerStarted, violations };
  };
  await page.route('**/__csp_probe.js', route => route.fulfill({
    contentType: 'text/javascript', body: `globalThis.cspProbe = (${probe.toString()})();`,
  }));
  await page.evaluate(url => new Promise((resolve, reject) => {
    const script = document.createElement('script');
    script.src = url;
    script.onload = () => resolve();
    script.onerror = () => reject(new Error('CSP probe did not load.'));
    document.head.append(script);
  }), new URL('__csp_probe.js', page.url()).href);
  const outcome = await page.evaluate(() => globalThis.cspProbe);
  expect(outcome.injected).toBe(false);
  expect(outcome.evaluated).toBe(false);
  expect(outcome.connected).toBe(false);
  expect(outcome.blobWorkerStarted).toBe(false);
  expect(outcome.violations).toContain('script-src-elem');
  expect(outcome.violations).toContain('script-src');
  expect(outcome.violations).toContain('connect-src');
  expect(outcome.violations).toContain('worker-src');
});
