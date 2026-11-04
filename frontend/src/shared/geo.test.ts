import { describe, expect, it } from 'vitest';

import { bboxFromCorners, bboxQuery, bboxToPolygon, formatDegrees } from './geo';

describe('bboxFromCorners', () => {
  it('orders west south east north regardless of click order', () => {
    expect(bboxFromCorners({ lng: 12, lat: 23 }, { lng: 10, lat: 20 })).toEqual({
      west: 10,
      south: 20,
      east: 12,
      north: 23,
    });
  });
});

describe('bboxQuery', () => {
  it('serializes the four edges as strings', () => {
    expect(bboxQuery({ west: -102.8, south: 31.6, east: -101.6, north: 32.4 })).toEqual({
      west: '-102.8',
      south: '31.6',
      east: '-101.6',
      north: '32.4',
    });
  });
});

describe('bboxToPolygon', () => {
  it('closes a rectangle ring', () => {
    const polygon = bboxToPolygon({ west: 10, south: 20, east: 12, north: 23 });
    expect(polygon.geometry.coordinates[0]).toEqual([
      [10, 20],
      [12, 20],
      [12, 23],
      [10, 23],
      [10, 20],
    ]);
  });
});

describe('formatDegrees', () => {
  it('shows two decimal places', () => {
    expect(formatDegrees(51.0447)).toBe('51.04');
  });
});
