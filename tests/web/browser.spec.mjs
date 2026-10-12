/** @file User-facing checks for the real educational browser solver. */
import { test, expect } from '@playwright/test';

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
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(2);
  expect((await page.locator('.root-full code').allTextContents()).map(Number)[1]).toBeCloseTo(Math.SQRT2, 12);
  expect(errors).toEqual([]);
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
  await expect(page.locator('#series-expansion')).toBeHidden();
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
  for (const id of ['series-parameters-details', 'series-sum-details', 'series-terms-details', 'series-math-details']) {
    await page.locator(`#${id} > summary`).click();
  }
  await expect(page.locator('#series-full-sum')).toBeVisible();
  await expect(page.locator('#series-terms tr')).toHaveCount(9);
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.setViewportSize({ width: 320, height: 720 });
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.screenshot({ path: 'out/screenshots/mobile.png', fullPage: true });
  for (const id of ['series-parameters-details', 'series-sum-details', 'series-terms-details', 'series-math-details']) {
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
