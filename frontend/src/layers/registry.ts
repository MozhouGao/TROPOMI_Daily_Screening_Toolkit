import type { PlumeSource } from '@/api/types';
import type { TranslationKey } from '@/i18n';

export type LayerKind = 'catalog' | 'screening' | 'facility' | 'inventory';

export type LayerDefinition = {
  id: PlumeSource;
  kind: LayerKind;
  labelKey: TranslationKey;
  freshnessKey: TranslationKey;
  disclaimerKey: TranslationKey;
  color: string;
  defaultEnabled: boolean;
};

export const LAYER_REGISTRY: readonly LayerDefinition[] = [
  {
    id: 'carbonmapper',
    kind: 'catalog',
    labelKey: 'layers.carbonmapper',
    freshnessKey: 'layers.carbonmapper.freshness',
    disclaimerKey: 'detail.independent',
    color: '#d7263d',
    defaultEnabled: true,
  },
  {
    id: 'mars',
    kind: 'catalog',
    labelKey: 'layers.mars',
    freshnessKey: 'layers.mars.freshness',
    disclaimerKey: 'detail.independent',
    color: '#1b6ca8',
    defaultEnabled: true,
  },
  {
    id: 'ghgsat',
    kind: 'catalog',
    labelKey: 'layers.ghgsat',
    freshnessKey: 'layers.ghgsat.freshness',
    disclaimerKey: 'detail.independent',
    color: '#7c3aed',
    defaultEnabled: true,
  },
  {
    id: 'cusworth',
    kind: 'catalog',
    labelKey: 'layers.cusworth',
    freshnessKey: 'layers.cusworth.freshness',
    disclaimerKey: 'detail.campaign',
    color: '#be185d',
    defaultEnabled: false,
  },
  {
    id: 'ghgrp_us',
    kind: 'facility',
    labelKey: 'layers.ghgrp_us',
    freshnessKey: 'layers.ghgrp_us.freshness',
    disclaimerKey: 'detail.facility',
    color: '#b45309',
    defaultEnabled: true,
  },
  {
    id: 'ghgrp_ca',
    kind: 'facility',
    labelKey: 'layers.ghgrp_ca',
    freshnessKey: 'layers.ghgrp_ca.freshness',
    disclaimerKey: 'detail.facility',
    color: '#c2410c',
    defaultEnabled: true,
  },
  {
    id: 'nesser',
    kind: 'facility',
    labelKey: 'layers.nesser',
    freshnessKey: 'layers.nesser.freshness',
    disclaimerKey: 'detail.inversion',
    color: '#0f766e',
    defaultEnabled: false,
  },
  {
    id: 'canada_nir',
    kind: 'inventory',
    labelKey: 'layers.canada_nir',
    freshnessKey: 'layers.canada_nir.freshness',
    disclaimerKey: 'detail.grid',
    color: '#4b5563',
    defaultEnabled: false,
  },
  {
    id: 'methanesat',
    kind: 'catalog',
    labelKey: 'layers.methanesat',
    freshnessKey: 'layers.methanesat.freshness',
    disclaimerKey: 'detail.methanesat',
    color: '#0e7490',
    defaultEnabled: false,
  },
  {
    id: 'tropomi',
    kind: 'screening',
    labelKey: 'layers.tropomi',
    freshnessKey: 'layers.tropomi.freshness',
    disclaimerKey: 'tropomi.attribution',
    color: '#415f4a',
    defaultEnabled: true,
  },
];

export const CATALOG_LAYERS = LAYER_REGISTRY.filter((layer) => layer.kind !== 'screening');

export const layerById = (id: PlumeSource) =>
  LAYER_REGISTRY.find((layer) => layer.id === id) ?? LAYER_REGISTRY[0]!;

export const defaultEnabledLayers = () =>
  LAYER_REGISTRY.filter((layer) => layer.defaultEnabled).map((layer) => layer.id);

export const countKeyForKind = (kind: LayerKind): TranslationKey => {
  if (kind === 'facility') {
    return 'layers.countFacilities';
  }
  if (kind === 'inventory') {
    return 'layers.countCells';
  }
  return 'layers.count';
};
