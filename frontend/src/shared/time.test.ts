import { DateTime } from 'luxon';
import { describe, expect, it } from 'vitest';

import { formatObserved, isQuickRange, quickRange, rangeQuery } from './time';

describe('quickRange', () => {
  it('includes the end day in the span', () => {
    const end = DateTime.utc(2026, 11, 1);
    const range = quickRange(30, end);
    expect(range.end.toISODate()).toBe('2026-11-01');
    expect(range.start.toISODate()).toBe('2026-10-03');
    expect(isQuickRange(range, 30)).toBe(false);
    expect(isQuickRange(quickRange(30), 30)).toBe(true);
  });
});

describe('rangeQuery', () => {
  it('emits ISO dates for the API', () => {
    const range = {
      start: DateTime.utc(2025, 6, 1),
      end: DateTime.utc(2025, 6, 30),
    };
    expect(rangeQuery(range)).toEqual({ start: '2025-06-01', end: '2025-06-30' });
  });
});

describe('formatObserved', () => {
  it('formats UTC timestamps', () => {
    expect(formatObserved('2025-06-25T18:13:04.330Z')).toBe('2025-06-25 18:13');
  });
});
