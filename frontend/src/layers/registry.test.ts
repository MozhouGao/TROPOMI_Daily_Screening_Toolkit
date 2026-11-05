import { describe, expect, it } from 'vitest';

import { CATALOG_LAYERS, LAYER_REGISTRY, layerById } from './registry';

describe('LAYER_REGISTRY', () => {
  it('registers each source once', () => {
    const ids = LAYER_REGISTRY.map((layer) => layer.id);
    expect(ids).toEqual([
      'carbonmapper',
      'mars',
      'ghgsat',
      'cusworth',
      'ghgrp_us',
      'ghgrp_ca',
      'nesser',
      'canada_nir',
      'methanesat',
      'tropomi',
    ]);
    expect(new Set(ids).size).toBe(ids.length);
  });

  it('keeps screening layers out of the catalog fetch list', () => {
    expect(CATALOG_LAYERS.map((layer) => layer.id)).toEqual([
      'carbonmapper',
      'mars',
      'ghgsat',
      'cusworth',
      'ghgrp_us',
      'ghgrp_ca',
      'nesser',
      'canada_nir',
      'methanesat',
    ]);
  });

  it('looks up a layer by id', () => {
    expect(layerById('mars').labelKey).toBe('layers.mars');
  });
});
