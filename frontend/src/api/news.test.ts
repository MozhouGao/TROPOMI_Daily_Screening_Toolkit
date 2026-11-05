import { describe, expect, it } from 'vitest';

import { newsUrl } from './news';

describe('newsUrl', () => {
  it('puts the bbox on the news path without a time window', () => {
    expect(newsUrl({ west: -102.8, south: 31.6, east: -101.6, north: 32.4 })).toBe(
      '/api/news?west=-102.8&south=31.6&east=-101.6&north=32.4',
    );
  });
});
