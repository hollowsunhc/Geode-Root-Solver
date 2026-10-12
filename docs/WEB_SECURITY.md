# Browser security and workload controls

The public lab is a static site. Solving happens in a disposable Web Worker on the
visitor's device. There is no calculation endpoint, account, cookie, analytics
request, or stored input. The page downloads its own HTML, scripts, stylesheet,
and WebAssembly module; it does not send coefficients to a server.

## Enforced limits

| Resource | Limit |
| --- | --- |
| Coefficient text | 2,048 characters |
| Coefficients | 17, including trailing zeros; degree at most 16 |
| One coefficient | 64 characters; decimal notation only |
| Explicit decimal exponent | -100 to 100 |
| Working precision | 128, 256, or 512 bits |
| Requested error | `10^-d`, with `d` from 1 to 30, 60, or 120 respectively |
| Module loading | 20 seconds |
| One calculation | 5 seconds after the module is ready |
| Concurrent calculations | One worker per solver client |
| WebAssembly linear memory | 128 MiB maximum, including its 2 MiB stack |
| Analytical expansion | Orders 0 through 8 per root; 7 precision levels, at most 2 center checks each |
| Worked products | At most 12 per order; factorial arguments at most 128 |
| Exact series fraction | At most 4,096 printed characters |

The page and worker use the same request validator. The C++ entry point independently
checks the numeric grammar and bounds before constructing exact rational values. The page
never evaluates input as JavaScript, HTML, or a symbolic expression.

Series generation runs in the same worker and shares its five-second deadline
and memory ceiling. Exact fractions are checked during coefficient arithmetic,
before conversion to decimal text. A controlled fraction-size or convergence
failure omits only that root's formula; a watchdog timeout or allocation failure
stops the complete request. The client checks formula counts, numeric formats,
string lengths, and the number of steps before rendering. It also checks product
counts, tuple lengths, nonnegative integral exponents, total order, and positive
integer factorial values. Formula values are
inserted as text, never executed as code or interpreted as HTML.

Worked-product enumeration stops as soon as twelve products have been emitted;
it does not generate an unbounded list and then truncate it. Full grouped
coefficients use polynomial convolution independently of this display cap.
Identified exact roots follow the same workload limits and retain an expansion;
linear equations show zero higher-order terms without artificial work.

Cancel, a superseding request, timeout, successful completion, and worker errors
all terminate the worker. This releases its Wasm instance and prevents a stalled
synchronous solve from blocking the interface. There are no automatic retries.
Browser timers can be delayed when a tab is suspended; five seconds is a watchdog
deadline, not an operating-system CPU quota.

The static page supplies a restrictive meta Content Security Policy. Scripts,
styles, workers, and connections are confined to the site's own origin. Inline
scripts and JavaScript `eval` are prohibited; WebAssembly compilation is allowed
explicitly. User-provided data is inserted through text nodes.

## Run the verifier

After building the web output:

```sh
npm ci
npx playwright install chromium
npm run verify:security
npm run test:browser
```

The security command tests malformed numbers, expression and HTML payloads,
oversized inputs, precision/exponent bounds, and bypasses that call the compiled
entry point directly. It checks the actual Wasm memory ceiling. Browser tests
exercise real worker lifetimes, cancellation, timeout, recovery, CSP, and network
behavior. A controlled worker fixture that deliberately stalls makes watchdog
tests deterministic; normal solving uses the real compiled module.

`npm test` also checks root correctness, including an independent rational check
of the 120-digit square-root result and exact decimal output for very large roots.
It also checks preservation of decimal repeated roots, hyper-Catalan series
coefficients against known values, and independently verifies partial-sum error
bounds with exact integer arithmetic. Walkthrough tests check explicit factorial
quotients, zero and exact roots, and a quartic order with 45 products, confirming
that displaying only 12 does not truncate the calculated coefficient.
Graph tests exercise large and degenerate axis ranges so floating-point rounding
cannot turn tick generation into an unbounded loop. `npm audit` checks development dependencies.
Passing these checks is evidence for the tested controls, not a promise that every
possible security flaw has been eliminated.

## Scope and remaining limits

- These controls protect ordinary visitors from accidental excessive work. A
  visitor who modifies JavaScript, directly calls Wasm, or opens many tabs can
  consume their own resources. There is no backend compute budget to exhaust.
- The memory ceiling covers Wasm linear memory, not all browser or JavaScript
  allocations. Large valid inputs can still reach a time or memory limit.
- GitHub Pages controls hosting bandwidth, response headers, and platform abuse
  protection. Client-side validation is not a server-side rate limiter.
- A meta CSP does not support `frame-ancestors` or reporting, and does not set CSP
  response headers on worker files. This site has no privileged actions to protect
  from clickjacking. A future server, authentication feature, or external input
  source would require a new threat review and server-enforced controls.
- Decimal input coefficients are exact rationals. Computed roots and partial sums
  are finite-precision approximations with separate absolute error bounds.
  Graphs use ordinary floating-point arithmetic and are explanatory, not proof
  of root accuracy. The finite series preview need not achieve the root solver's
  requested accuracy, and some valid inputs have no formula within the series budget.

## Dependency and source distribution

The build pins Emscripten 4.0.23, GMP 6.3.0, and MPFR 4.2.2. It verifies source
archive hashes before extracting dependencies. The deployed site includes license
notices and `source/geode-browser-source.tar.gz`, containing application sources,
the exact GMP/MPFR source archives, and the build script. Users can modify and
rebuild the application and libraries. No secrets or local build caches belong
in the publication branch.
