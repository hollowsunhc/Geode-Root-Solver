/** @file Approximate, interactive polynomial plotting in browser floating-point arithmetic. */

/** @brief Produce a short root preview without a floating-point conversion. */
export function shortRootLabel(root) {
  const match = root.match(/^([+-]?)(\d+)(?:\.(\d*))?(?:e([+-]?\d+))?$/i);
  if (!match) return root.slice(0, 16);
  const sign = match[1] === '-' ? '−' : '';
  const combined = match[2] + (match[3] ?? '');
  const leadingZeros = combined.match(/^0*/)[0].length;
  const digits = combined.slice(leadingZeros).replace(/0+$/, '');
  if (!digits) return `${sign}0`;
  const point = match[2].length + Number(match[4] ?? 0) - leadingZeros;
  let label;
  if (point >= -3 && point <= 6) {
    label = point <= 0 ? `0.${'0'.repeat(-point)}${digits}`
      : point >= digits.length ? digits + '0'.repeat(point - digits.length)
        : `${digits.slice(0, point)}.${digits.slice(point)}`;
    return sign + (label.length > 16 ? `${label.slice(0, 16)}…` : label);
  }
  label = digits[0] + (digits.length > 1 ? `.${digits.slice(1, 13)}${digits.length > 13 ? '…' : ''}` : '');
  return `${sign}${label} × 10^${String(point - 1).replace('-', '−')}`;
}

/** @brief Evaluate ascending coefficients with Horner's method for visualization only. */
export function evaluatePolynomial(coefficients, x) {
  let value = 0;
  for (let i = coefficients.length - 1; i >= 0; i -= 1) value = value * x + coefficients[i];
  return value;
}

/** @brief Choose evenly spaced, human-readable ticks for a finite axis range. */
export function axisTicks(min, max, count = 6) {
  if (!Number.isFinite(min) || !Number.isFinite(max) || min >= max) return [];
  count = Number.isFinite(count) ? Math.max(2, Math.min(12, count)) : 6;
  const roughStep = (max - min) / count;
  const magnitude = 10 ** Math.floor(Math.log10(roughStep));
  const step = ([1, 2, 5, 10].find(value => value * magnitude >= roughStep) ?? 10) * magnitude;
  if (!Number.isFinite(step) || step <= 0) return [];
  const first = Math.ceil(min / step);
  if (!Number.isFinite(first)) return [];
  const result = [];
  // A bounded offset still advances when a large coordinate plus one rounds back to itself.
  for (let offset = 0; offset <= count + 2; offset += 1) {
    const value = (first + offset) * step;
    if (value > max) break;
    if (!Number.isFinite(value) || value < min || result.at(-1)?.value === value) continue;
    const scientific = value !== 0 && (Math.abs(value) < 0.001 || Math.abs(value) >= 100000);
    const label = scientific ? value.toExponential(1) : value.toFixed(Math.max(0, Math.min(12, -Math.floor(Math.log10(step)))));
    const formatted = label.replace('-', '−');
    if (result.at(-1)?.label !== formatted) result.push({ value, label: formatted });
  }
  return result;
}

/** @brief Zoom around a world-coordinate anchor while retaining a usable numerical scale. */
export function zoomView(view, factor, x = (view.xMin + view.xMax) / 2, y = (view.yMin + view.yMax) / 2) {
  const xSpan = (view.xMax - view.xMin) * factor;
  const ySpan = (view.yMax - view.yMin) * factor;
  if (!Number.isFinite(xSpan + ySpan) || Math.min(xSpan, ySpan) < 1e-12 || Math.max(xSpan, ySpan) > 1e12) return { ...view };
  const next = {
    xMin: x + (view.xMin - x) * factor,
    xMax: x + (view.xMax - x) * factor,
    yMin: y + (view.yMin - y) * factor,
    yMax: y + (view.yMax - y) * factor,
  };
  return Object.values(next).every(Number.isFinite) && next.xMin < next.xMax && next.yMin < next.yMax ? next : { ...view };
}

/** @brief Sample a bounded SVG curve; overflow breaks a segment instead of entering the path. */
export function curvePath(coefficients, view, width, height) {
  if (width <= 0 || height <= 0) return '';
  const samples = Math.ceil(width);
  const points = [];
  let connected = false;
  for (let i = 0; i <= samples; i += 1) {
    const x = view.xMin + (i / samples) * (view.xMax - view.xMin);
    const value = evaluatePolynomial(coefficients, x);
    if (!Number.isFinite(value)) {
      connected = false;
      continue;
    }
    const y = Math.max(-2 * height, Math.min(3 * height, ((view.yMax - value) / (view.yMax - view.yMin)) * height));
    points.push(`${connected ? 'L' : 'M'}${((i / samples) * width).toFixed(2)},${y.toFixed(2)}`);
    connected = true;
  }
  return points.join('');
}

/** @brief Create an SVG element without interpreting any text as markup. */
function svgElement(name, attributes, text) {
  const element = document.createElementNS('http://www.w3.org/2000/svg', name);
  for (const [key, value] of Object.entries(attributes)) element.setAttribute(key, String(value));
  if (text !== undefined) element.textContent = text;
  return element;
}

/** @brief Show a short, explicitly approximate coordinate in the plot readout. */
function coordinate(value) {
  return Number.isFinite(value) ? Number(value.toPrecision(6)).toString().replace('-', '−') : 'outside graph range';
}

/** @brief A local SVG plot with mouse, touch, and keyboard navigation. */
export class PolynomialPlot {
  /** @brief Bind controls and resize behavior to the plot's existing semantic markup. */
  constructor(svg) {
    this.svg = svg;
    this.coefficients = null;
    this.roots = [];
    this.defaultView = { xMin: -4, xMax: 4, yMin: -12, yMax: 12 };
    this.view = { ...this.defaultView };
    this.readout = document.querySelector('#plot-readout');
    this.range = document.querySelector('#plot-range');
    this.message = document.querySelector('#plot-message');
    this.drag = null;
    this.resizeObserver = new ResizeObserver(() => this.render());
    this.resizeObserver.observe(svg);

    document.querySelector('#zoom-in').addEventListener('click', () => this.zoom(0.7));
    document.querySelector('#zoom-out').addEventListener('click', () => this.zoom(1 / 0.7));
    document.querySelector('#reset-view').addEventListener('click', () => this.reset());
    svg.addEventListener('keydown', event => this.onKey(event));
    svg.addEventListener('pointerdown', event => {
      if (event.button !== 0) return;
      this.onPointer(event);
      this.drag = { x: event.clientX, y: event.clientY, view: { ...this.view } };
      svg.setPointerCapture(event.pointerId);
      svg.classList.add('is-dragging');
    });
    svg.addEventListener('pointermove', event => this.onPointer(event));
    const finishDrag = () => {
      if (this.drag) this.announceRange();
      this.drag = null;
      svg.classList.remove('is-dragging');
    };
    svg.addEventListener('pointerup', finishDrag);
    svg.addEventListener('pointercancel', finishDrag);
    svg.addEventListener('lostpointercapture', finishDrag);
    svg.addEventListener('pointerleave', () => this.clearHover());
    svg.addEventListener('wheel', event => {
      if (document.activeElement !== svg && !event.ctrlKey) return;
      event.preventDefault();
      const point = this.pointAt(event);
      this.zoom(event.deltaY > 0 ? 1.15 : 1 / 1.15, point.x, point.y);
    }, { passive: false });
  }

  /** @brief Replace the approximate curve and clear markers for previous computed roots. */
  setPolynomial(coefficients, view = this.defaultView) {
    this.coefficients = coefficients ? coefficients.map(Number) : null;
    this.roots = [];
    this.defaultView = { ...view };
    this.view = { ...view };
    this.clearHover();
    this.render();
  }

  /** @brief Convert root strings solely for approximate graph markers, never numerical output. */
  setRoots(roots) {
    this.roots = roots.map(Number).filter(Number.isFinite);
    this.render();
  }

  /** @brief Restore the selected example's original viewport. */
  reset() {
    this.view = { ...this.defaultView };
    this.render();
    this.announceRange();
  }

  /** @brief Change graph scale around an optional anchor. */
  zoom(factor, x, y) {
    this.view = zoomView(this.view, factor, x, y);
    this.render();
    this.announceRange();
  }

  /** @brief Give assistive technology a concise viewport update after navigation. */
  announceRange() {
    this.range.textContent = `Approximate horizontal range ${coordinate(this.view.xMin)} to ${coordinate(this.view.xMax)}. Vertical range ${coordinate(this.view.yMin)} to ${coordinate(this.view.yMax)}.`;
  }

  /** @brief Translate a pointer position into graph coordinates. */
  pointAt(event) {
    const rectangle = this.svg.getBoundingClientRect();
    const px = event.clientX - rectangle.left - this.left;
    const py = event.clientY - rectangle.top - this.top;
    return {
      x: this.view.xMin + (px / this.width) * (this.view.xMax - this.view.xMin),
      y: this.view.yMax - (py / this.height) * (this.view.yMax - this.view.yMin),
      px,
      py,
    };
  }

  /** @brief Pan while dragging, otherwise inspect one approximate curve sample. */
  onPointer(event) {
    if (this.drag) {
      const original = this.drag.view;
      const dx = ((event.clientX - this.drag.x) / this.width) * (original.xMax - original.xMin);
      const dy = ((event.clientY - this.drag.y) / this.height) * (original.yMax - original.yMin);
      this.view = { xMin: original.xMin - dx, xMax: original.xMax - dx, yMin: original.yMin + dy, yMax: original.yMax + dy };
      this.render();
      return;
    }
    const point = this.pointAt(event);
    if (!this.coefficients || point.px < 0 || point.px > this.width || point.py < 0 || point.py > this.height) {
      this.clearHover();
      return;
    }
    const value = evaluatePolynomial(this.coefficients, point.x);
    this.readout.textContent = `x ≈ ${coordinate(point.x)}    f(x) ≈ ${coordinate(value)}`;
    const hover = this.svg.querySelector('#plot-hover');
    hover.replaceChildren(svgElement('line', { x1: point.px + this.left, x2: point.px + this.left, y1: this.top, y2: this.top + this.height, class: 'hover-line' }));
    if (value >= this.view.yMin && value <= this.view.yMax) {
      const y = this.top + ((this.view.yMax - value) / (this.view.yMax - this.view.yMin)) * this.height;
      hover.append(svgElement('circle', { cx: point.px + this.left, cy: y, r: 5, class: 'hover-dot' }));
    }
  }

  /** @brief Clear the cursor without retaining a sample from an old polynomial. */
  clearHover() {
    this.svg.querySelector('#plot-hover').replaceChildren();
    this.readout.textContent = 'Hover or tap to inspect the curve';
  }

  /** @brief Provide the same viewport operations without a pointer. */
  onKey(event) {
    if (['+', '='].includes(event.key)) this.zoom(0.7);
    else if (event.key === '-') this.zoom(1 / 0.7);
    else if (event.key === '0' || event.key === 'Home') this.reset();
    else if (['ArrowLeft', 'ArrowRight', 'ArrowUp', 'ArrowDown'].includes(event.key)) {
      const dx = (this.view.xMax - this.view.xMin) * 0.2 * (event.key === 'ArrowLeft' ? -1 : event.key === 'ArrowRight' ? 1 : 0);
      const dy = (this.view.yMax - this.view.yMin) * 0.2 * (event.key === 'ArrowDown' ? -1 : event.key === 'ArrowUp' ? 1 : 0);
      this.view = { xMin: this.view.xMin + dx, xMax: this.view.xMax + dx, yMin: this.view.yMin + dy, yMax: this.view.yMax + dy };
      this.render();
      this.announceRange();
    } else return;
    event.preventDefault();
  }

  /** @brief Draw the viewport grid, sampled curve, and approximate root markers. */
  render() {
    const fullWidth = Math.max(240, this.svg.clientWidth);
    const fullHeight = Math.max(180, this.svg.clientHeight);
    this.left = fullWidth < 400 ? 48 : 60;
    this.top = 18;
    this.width = fullWidth - this.left - 22;
    this.height = fullHeight - this.top - 34;
    this.svg.setAttribute('viewBox', `0 0 ${fullWidth} ${fullHeight}`);
    const rectangle = this.svg.querySelector('#plot-clip rect');
    for (const [attribute, value] of Object.entries({ x: this.left, y: this.top, width: this.width, height: this.height })) rectangle.setAttribute(attribute, value);
    const xPixel = x => this.left + ((x - this.view.xMin) / (this.view.xMax - this.view.xMin)) * this.width;
    const yPixel = y => this.top + ((this.view.yMax - y) / (this.view.yMax - this.view.yMin)) * this.height;
    const grid = this.svg.querySelector('#plot-grid');
    const axes = this.svg.querySelector('#plot-axes');
    const labels = this.svg.querySelector('#plot-labels');
    grid.replaceChildren();
    axes.replaceChildren();
    labels.replaceChildren();
    for (const tick of axisTicks(this.view.xMin, this.view.xMax, fullWidth < 400 ? 4 : 7)) {
      const x = xPixel(tick.value);
      grid.append(svgElement('line', { x1: x, x2: x, y1: this.top, y2: this.top + this.height }));
      labels.append(svgElement('text', { x, y: fullHeight - 10, 'text-anchor': 'middle' }, tick.label));
    }
    for (const tick of axisTicks(this.view.yMin, this.view.yMax, 5)) {
      const y = yPixel(tick.value);
      grid.append(svgElement('line', { x1: this.left, x2: this.left + this.width, y1: y, y2: y }));
      labels.append(svgElement('text', { x: this.left - 10, y: y + 4, 'text-anchor': 'end' }, tick.label));
    }
    if (this.view.xMin <= 0 && this.view.xMax >= 0) axes.append(svgElement('line', { x1: xPixel(0), x2: xPixel(0), y1: this.top, y2: this.top + this.height }));
    if (this.view.yMin <= 0 && this.view.yMax >= 0) axes.append(svgElement('line', { x1: this.left, x2: this.left + this.width, y1: yPixel(0), y2: yPixel(0), class: 'zero-axis' }));
    labels.append(svgElement('text', { x: fullWidth - 8, y: fullHeight - 10, class: 'axis-name' }, 'x'));
    labels.append(svgElement('text', { x: this.left - 10, y: 12, 'text-anchor': 'end', class: 'axis-name' }, 'f(x)'));
    const curve = this.svg.querySelector('#plot-curve');
    const path = this.coefficients ? curvePath(this.coefficients, this.view, this.width, this.height) : '';
    curve.setAttribute('d', path);
    curve.setAttribute('transform', `translate(${this.left} ${this.top})`);
    this.message.hidden = Boolean(path);
    this.message.textContent = this.coefficients ? 'This view exceeds the approximate graph’s range.' : 'Complete the coefficients to see your curve.';
    const markers = this.svg.querySelector('#plot-roots');
    markers.replaceChildren();
    for (const [index, root] of this.roots.entries()) {
      if (root < this.view.xMin || root > this.view.xMax || this.view.yMin > 0 || this.view.yMax < 0) continue;
      const marker = svgElement('circle', { cx: xPixel(root), cy: yPixel(0), r: 6, class: 'root-marker' });
      marker.append(svgElement('title', {}, `Approximate position of root ${index + 1}`));
      markers.append(marker);
    }
    this.svg.querySelector('#plot-hover').replaceChildren();
  }
}
