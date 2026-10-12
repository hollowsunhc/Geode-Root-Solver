/** @file User-facing checks for the real educational browser solver. */
import { test, expect } from '@playwright/test';

test('a visitor can solve, inspect precise roots, and invalidate stale results', async ({ page }) => {
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto('./');
  await expect(page.locator('#plot')).toBeVisible();
  await expect(page.locator('.root-row')).toHaveCount(0);
  await page.locator('#solve-button').click();
  await expect(page.locator('.root-row')).toHaveCount(3);
  expect((await page.locator('.root-full code').allTextContents()).map(Number)).toEqual([-2, 1, 3]);
  await page.locator('#coefficients').fill('-2, 0, 1');
  await expect(page.locator('.root-row')).toHaveCount(0);
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
  expect(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth)).toBe(true);
  await page.screenshot({ path: 'out/screenshots/mobile.png', fullPage: true });
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
