import { describe, expect, it } from 'vitest';

import type { Plume } from '@/api/types';

import { formatEmission, mergeNewestFirst, plumeKey, toFeatureCollection } from './plumes';

const plume = (overrides: Partial<Plume> = {}): Plume => ({
  id: 'p1',
  source: 'carbonmapper',
  observedAt: '2025-06-25T18:13:04Z',
  lon: -101.8,
  lat: 32.2,
  emissionKgHr: 554.5,
  uncertaintyKgHr: 105.1,
  instrument: 'tan',
  sector: '1B2',
  imageUrl: 'https://example.org/plume.png',
  attribution: 'Source: Carbon Mapper',
  details: {},
  ...overrides,
});

describe('plumeKey', () => {
  it('namespaces ids by source so layers cannot collide', () => {
    expect(plumeKey({ source: 'mars', id: 'p1' })).toBe('mars:p1');
  });
});

describe('formatEmission', () => {
  it('includes rounded uncertainty when present', () => {
    expect(formatEmission(plume())).toBe('555 ± 105');
  });

  it('is blank when emission is missing', () => {
    expect(formatEmission(plume({ emissionKgHr: null, uncertaintyKgHr: null }))).toBe('');
  });
});

describe('mergeNewestFirst', () => {
  it('flattens layer groups by observation time', () => {
    const older = plume({ id: 'old', observedAt: '2025-01-01T00:00:00Z' });
    const newer = plume({ id: 'new', source: 'mars', observedAt: '2025-08-01T00:00:00Z' });
    expect(mergeNewestFirst([[older], [newer]]).map((item) => item.id)).toEqual(['new', 'old']);
  });
});

describe('toFeatureCollection', () => {
  it('puts lon/lat on a point and the selection key on properties', () => {
    const collection = toFeatureCollection([plume()]);
    expect(collection.features[0]?.geometry.coordinates).toEqual([-101.8, 32.2]);
    expect(collection.features[0]?.properties).toEqual({ key: 'carbonmapper:p1', emission: 554.5 });
  });
});
