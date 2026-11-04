import type { FeatureCollection, Point } from 'geojson';
import { DateTime } from 'luxon';

import type { Plume } from '@/api/types';

export type PlumeProperties = { key: string; emission: number };

export const plumeKey = (plume: Pick<Plume, 'source' | 'id'>) => `${plume.source}:${plume.id}`;

export const toFeatureCollection = (
  plumes: Plume[],
): FeatureCollection<Point, PlumeProperties> => ({
  type: 'FeatureCollection',
  features: plumes.map((plume) => ({
    type: 'Feature',
    geometry: { type: 'Point', coordinates: [plume.lon, plume.lat] },
    properties: { key: plumeKey(plume), emission: plume.emissionKgHr ?? 0 },
  })),
});

const observedMillis = (plume: Plume) => DateTime.fromISO(plume.observedAt).toMillis();

export const mergeNewestFirst = (groups: Plume[][]) =>
  groups.flat().sort((a, b) => observedMillis(b) - observedMillis(a));

export const formatEmission = (plume: Pick<Plume, 'emissionKgHr' | 'uncertaintyKgHr'>) => {
  if (plume.emissionKgHr === null) {
    return '';
  }
  const emission = Math.round(plume.emissionKgHr).toLocaleString();
  if (plume.uncertaintyKgHr === null) {
    return emission;
  }
  return `${emission} ± ${Math.round(plume.uncertaintyKgHr).toLocaleString()}`;
};
