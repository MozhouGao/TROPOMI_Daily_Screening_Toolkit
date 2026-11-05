import { describe, expect, it } from 'vitest';

import { DateTime } from 'luxon';

import { plumesUrl } from './plumes';

describe('plumesUrl', () => {
  it('puts bbox and dates on the catalog path', () => {
    const url = plumesUrl(
      'carbonmapper',
      { west: -102.8, south: 31.6, east: -101.6, north: 32.4 },
      { start: DateTime.utc(2025, 6, 1), end: DateTime.utc(2025, 6, 30) },
    );
    expect(url).toBe(
      '/api/plumes/carbonmapper?west=-102.8&south=31.6&east=-101.6&north=32.4&start=2025-06-01&end=2025-06-30',
    );
  });
});
