import { useQueries } from '@tanstack/react-query';

import type { PlumeCollection, PlumeSource } from '@/api/types';
import { CATALOG_LAYERS } from '@/layers/registry';
import { bboxQuery, type BBox } from '@/shared/geo';
import { rangeQuery, type DateRange } from '@/shared/time';

const FIVE_MINUTES = 5 * 60 * 1000;

export const plumesUrl = (source: PlumeSource, bbox: BBox, range: DateRange) =>
  `/api/plumes/${source}?${new URLSearchParams({ ...bboxQuery(bbox), ...rangeQuery(range) })}`;

const readErrorMessage = async (response: Response) => {
  const body = (await response.json().catch(() => null)) as { detail?: unknown } | null;
  return typeof body?.detail === 'string' ? body.detail : `HTTP ${response.status}`;
};

export const fetchPlumes = async (
  source: PlumeSource,
  bbox: BBox,
  range: DateRange,
): Promise<PlumeCollection> => {
  const response = await fetch(plumesUrl(source, bbox, range));
  if (response.ok) {
    return (await response.json()) as PlumeCollection;
  }
  throw new Error(await readErrorMessage(response));
};

export type LayerQuery = {
  source: PlumeSource;
  enabled: boolean;
  data: PlumeCollection | undefined;
  isLoading: boolean;
  error: Error | null;
};

const usePlumeLayers = (
  enabledLayers: PlumeSource[],
  bbox: BBox | null,
  range: DateRange,
): LayerQuery[] => {
  const results = useQueries({
    queries: CATALOG_LAYERS.map((layer) => ({
      queryKey: ['plumes', layer.id, bbox, rangeQuery(range)],
      queryFn: () => fetchPlumes(layer.id, bbox as BBox, range),
      enabled: bbox !== null && enabledLayers.includes(layer.id),
      staleTime: FIVE_MINUTES,
      retry: 1,
    })),
  });

  return CATALOG_LAYERS.map((layer, index) => {
    const result = results[index]!;
    return {
      source: layer.id,
      enabled: enabledLayers.includes(layer.id),
      data: result.data,
      isLoading: result.isFetching,
      error: result.error,
    };
  });
};

export default usePlumeLayers;
