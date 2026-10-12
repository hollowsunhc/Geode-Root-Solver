/** @file User-facing checks for the real educational browser solver. */
import { test, expect } from '@playwright/test';

/** @brief Verify the mixed cubic product against its exact root-specific parameters. */
async function expectMixedProduct(page) {
  const parameters = (await page.locator('#series-parameters dd code').allTextContents()).slice(2, 4);
  const rational = value => {
    const [numerator, denominator = '1'] = value.split('/');
    return [BigInt(numerator), BigInt(denominator)];
  };
  const [[n2, d2], [n3, d3]] = parameters.map(rational);
  const value = await page.locator('#series-product-value-exact').textContent();
  const [numerator, denominator] = rational(value);
  expect(numerator * d2 * d3).toBe(5n * n2 * n3 * denominator);
  await expect(page.locator('#series-coefficient')).toHaveText('5');
  await expect(page.locator('#series-factorial-ratio')).toHaveText('120 / 24');
  await expect(page.locator('#series-factorial')).toContainText('2×1+3×1');
  await expect(page.locator('#series-substitution-exact')).toHaveText(`5 × (${parameters[0]})^1 × (${parameters[1]})^1 = ${value}`);
  return { parameters, value };
}

test('a visitor can select a root, inspect its series, and invalidate stale results', async ({ page }) => {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto('./');
  await expect(page.locator('#plot')).toBeVisible();
  await expect(page.locator('.root-row')).toHaveCount(0);
  await expect(page.locator('#series-content')).toBeHidden();
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(3);
  expect((await page.locator('.root-full code').allTextContents()).map(Number)).toEqual([-2, 1, 3]);
  await expect(page.locator('#series-root option')).toHaveCount(3);
  await expect(page.locator('#series-expansion')).toBeVisible();
  await page.locator('#series-root').selectOption('1');
  await expect(page.locator('#series-root-index')).toHaveText('2');
  const selectedPreview = (await page.locator('#series-root-value').textContent()).split('≈ ')[1].replace('…', '');
  expect(Number(selectedPreview)).toBeCloseTo(1, 12);
  await page.locator('#series-order').focus();
  await page.keyboard.press('Home');
  await expect(page.locator('#series-steps tr')).toHaveCount(1);
  const firstSum = await page.locator('#series-full-sum').textContent();
  const firstBound = Number(await page.locator('#series-full-bound').textContent());
  await page.keyboard.press('End');
  await expect(page.locator('#series-steps tr')).toHaveCount(9);
  await expect(page.locator('#series-order-label')).toContainText('j = 8');
  expect(await page.locator('#series-full-sum').textContent()).not.toBe(firstSum);
  expect(Number(await page.locator('#series-full-bound').textContent())).toBeLessThan(firstBound);
  await page.locator('#series-panel').screenshot({ path: 'out/screenshots/series-desktop.png' });
  await page.locator('#coefficients').fill('-2, 0, 1');
  await expect(page.locator('.root-row')).toHaveCount(0);
  await expect(page.locator('#series-content')).toBeHidden();
  await expect(page.locator('#series-root option')).toHaveCount(0);
  await expect(page.locator('#series-steps tr')).toHaveCount(0);
  await expect(page.locator('#series-product option')).toHaveCount(0);
  await expect(page.locator('#series-coefficient')).toBeEmpty();
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(2);
  expect((await page.locator('.root-full code').allTextContents()).map(Number)[1]).toBeCloseTo(Math.SQRT2, 12);
  expect(errors).toEqual([]);
});

test('the cubic guide evaluates the mixed coefficient and substitutes the selected root', async ({ page }) => {
  await page.goto('./');
  await page.locator('#solve-button').click();
  await expect(page.locator('#series-expansion')).toBeVisible();
  await expect(page.locator('#series-product option')).toHaveCount(3);
  await expect(page.locator('#series-product option')).toHaveText(['m₂ = 2, m₃ = 0', 'm₂ = 1, m₃ = 1', 'm₂ = 0, m₃ = 2']);
  await page.locator('#series-product').focus();
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Enter');
  await expect(page.locator('#series-product')).toHaveValue('1');
  const first = await expectMixedProduct(page);
  await expect(page.locator('#series-exponents')).toContainText('j = 1 + 1 = 2');
  await expect(page.locator('#series-group-note')).toContainText('all 3 product values');
  await expect(page.locator('#series-partial-formula')).toContainText('R2=');
  await expect(page.locator('#series-partial-formula')).toContainText('s0+s1+s2');
  await expect(page.locator('#series-partial-formula mi[mathvariant="normal"]')).toHaveText('δ');
  await page.locator('#series-root').selectOption('1');
  await expect(page.locator('#series-product')).toHaveValue('1');
  const second = await expectMixedProduct(page);
  expect(second.parameters).not.toEqual(first.parameters);
  expect(second.value).not.toBe(first.value);
  await page.locator('#series-panel').screenshot({ path: 'out/screenshots/series-worked-cubic.png' });
});

test('guided examples teach repeated and absent real roots', async ({ page }) => {
  await page.goto('./');
  await page.locator('[data-example="repeated"]').click();
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(2);
  await page.locator('[data-example="no-real"]').click();
  await page.locator('#solve-button').click();
  await expect(page.locator('#results')).toContainText(/no real roots/i);
  await expect(page.locator('.root-row')).toHaveCount(0);
  await expect(page.locator('#series-content')).toBeHidden();
  await expect(page.locator('#series-empty')).toContainText('No real roots');
});

test('exact decimal coefficients keep a repeated rational root exact', async ({ page }) => {
  await page.goto('./');
  await page.locator('#coefficients').fill('0.01, -0.2, 1');
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(1);
  await expect(page.locator('#series-exact-value')).toHaveText('1/10');
  await expect(page.locator('#series-reduction-note')).toContainText('repeated factors removed');
  await expect(page.locator('#series-expansion')).toBeVisible();
  await expect(page.locator('#series-linear-note')).toContainText('S = 1 + 0 + 0 + ⋯');
  await expect(page.locator('#series-order')).toHaveValue('1');
  await expect(page.locator('#series-product option')).toHaveText(['Empty tuple (no parameters)']);
  await expect(page.locator('#series-coefficient')).toHaveText('1');
  await expect(page.locator('#series-product-value')).toHaveText('1');
  await page.locator('#series-order').focus();
  await page.keyboard.press('End');
  await expect(page.locator('#series-steps tr')).toHaveCount(9);
  await expect(page.locator('#series-product option')).toHaveCount(0);
  await expect(page.locator('#series-product-calculation')).toBeHidden();
  await expect(page.locator('#series-product-count')).toContainText('no products');
  await expect(page.locator('#series-group-value')).toHaveText('s₈ = 0');
  await expect(page.locator('#series-coefficient')).toBeEmpty();
  await expect(page.locator('#series-terms tr').nth(1).locator('td')).toHaveText(['0', '0']);
  expect(Number(await page.locator('#series-full-sum').textContent())).toBe(0.1);
});

test('a nonlinear zero root keeps a nontrivial expansion', async ({ page }) => {
  await page.goto('./');
  await page.locator('#coefficients').fill('0, -2, 0, 1');
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(3);
  await page.locator('#series-root').selectOption('1');
  await expect(page.locator('#series-exact-value')).toHaveText('0');
  await expect(page.locator('#series-expansion')).toBeVisible();
  await expect(page.locator('#series-linear-note')).toBeHidden();
  await expect(page.locator('#series-product option')).toHaveCount(3);
  await page.locator('#series-order').focus();
  await page.keyboard.press('Home');
  expect(Number(await page.locator('#series-full-sum').textContent())).not.toBe(0);
  await page.keyboard.press('End');
  await expect(page.locator('#series-steps tr')).toHaveCount(9);
});

test('a capped product list still explains the full grouped sum', async ({ page }) => {
  await page.goto('./');
  await page.locator('#coefficients').fill('-1, 0, 0, 0, 1');
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(2);
  await page.locator('#series-order').focus();
  await page.keyboard.press('End');
  await expect(page.locator('#series-product option')).toHaveCount(12);
  await expect(page.locator('#series-product-cap')).toHaveText('Showing the first 12 of 45 products. s₈ and every partial sum include all 45.');
  await expect(page.locator('#series-group-note')).toContainText('all 45 product values');
  await expect(page.locator('#series-group-formula')).toContainText('⋯');
});

test('high precision remains readable on a narrow screen', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('./');
  await page.locator('#coefficients').fill('-2 0 1');
  await page.locator('#bits').selectOption('512');
  await page.locator('#digits').fill('120');
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(2);
  const roots = await page.locator('.root-full code').allTextContents();
  expect(roots[1].length).toBeGreaterThan(120);
  await page.locator('.root-full summary').first().click();
  await expect(page.locator('.root-full code').first()).toBeVisible();
  await page.locator('#series-order').focus();
  await page.keyboard.press('End');
  for (const id of ['series-parameters-details', 'series-product-details', 'series-group-details', 'series-sum-details', 'series-terms-details', 'series-math-details']) {
    await page.locator(`#${id} > summary`).click();
  }
  await expect(page.locator('#series-full-sum')).toBeVisible();
  await expect(page.locator('#series-terms tr')).toHaveCount(9);
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.setViewportSize({ width: 320, height: 720 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.screenshot({ path: 'out/screenshots/mobile.png', fullPage: true });
  for (const id of ['series-parameters-details', 'series-product-details', 'series-group-details', 'series-sum-details', 'series-terms-details', 'series-math-details']) {
    await page.locator(`#${id} > summary`).click();
  }
  await page.locator('#series-panel').screenshot({ path: 'out/screenshots/series-mobile.png' });
});

test('graph tools remain keyboard accessible', async ({ page }) => {
  await page.goto('./');
  await page.locator('#zoom-in').focus();
  await expect(page.locator('#zoom-in')).toBeFocused();
  await page.keyboard.press('Enter');
  await page.locator('#reset-view').focus();
  await page.keyboard.press('Enter');
  await expect(page.locator('#plot')).toBeVisible();
  await page.screenshot({ path: 'out/screenshots/desktop.png', fullPage: true });
});

test('the quadratic story adds Catalan terms and reaches the other root by reflection', async ({ page }) => {
  await page.goto('./');
  await expect(page.locator('#quadratic-parameter-label')).toContainText('2/9');
  await expect(page.locator('#quadratic-exact')).toContainText('3/2');
  await page.locator('#quadratic-terms').focus();
  await page.keyboard.press('Home');
  await expect(page.locator('#quadratic-partial')).toHaveText('1.000000000');
  await page.keyboard.press('ArrowRight');
  await page.keyboard.press('ArrowRight');
  expect(Number(await page.locator('#quadratic-partial').textContent())).toBeCloseTo(107 / 81, 8);
  await page.locator('#quadratic-root').selectOption('upper');
  await expect(page.locator('#quadratic-center')).toHaveText('4.5');
  await expect(page.locator('#quadratic-scale')).toHaveText('−1');
  await expect(page.locator('#quadratic-exact')).toContainText('3');
  expect(Number(await page.locator('#quadratic-target').textContent())).toBe(3);
  expect(Number(await page.locator('#quadratic-partial').textContent())).toBeCloseTo(515 / 162, 8);
  const before = Number(await page.locator('#quadratic-tail').textContent());
  await page.locator('#quadratic-terms').focus();
  await page.keyboard.press('End');
  expect(Number(await page.locator('#quadratic-tail').textContent())).toBeLessThan(before);
  expect(Number(await page.locator('#quadratic-partial').textContent())).toBeCloseTo(3, 6);
  await page.locator('#quadratic-reset').click();
  await expect(page.locator('#quadratic-root')).toHaveValue('lower');
  await expect(page.locator('#quadratic-terms')).toHaveValue('9');
  await expect(page.locator('#quadratic-exact')).toContainText('3/2');
  await page.locator('#quadratic-story').screenshot({ path: 'out/screenshots/quadratic-desktop.png' });
});

test('the quadratic story shows a convergent boundary and stays usable on a small screen', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('./');
  await page.locator('#quadratic-parameter').focus();
  await page.keyboard.press('End');
  await expect(page.locator('#quadratic-parameter-label')).toContainText('1/4');
  await expect(page.locator('#quadratic-state')).toContainText('double root');
  await expect(page.locator('#quadratic-state')).toContainText('still converges');
  expect(Number(await page.locator('#quadratic-target').textContent())).toBe(2);
  await page.locator('#quadratic-terms').focus();
  await page.keyboard.press('Home');
  await page.keyboard.press('ArrowRight');
  await page.keyboard.press('ArrowRight');
  expect(Number(await page.locator('#quadratic-partial').textContent())).toBe(1.375);
  expect(Number(await page.locator('#quadratic-tail').textContent())).toBe(0.625);
  await page.locator('#quadratic-root').selectOption('upper');
  expect(Number(await page.locator('#quadratic-partial').textContent())).toBe(2.625);
  await page.locator('#quadratic-bound-details > summary').click();
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.setViewportSize({ width: 320, height: 720 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.locator('#quadratic-story').screenshot({ path: 'out/screenshots/quadratic-mobile.png' });
  await page.locator('#quadratic-parameter').focus();
  await page.keyboard.press('ArrowLeft');
  await expect(page.locator('#quadratic-state')).toContainText('two distinct roots');
  await page.locator('#coefficients').fill('-2 0 1');
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(2);
  await expect(page.locator('#quadratic-parameter')).toHaveValue('899');
});
