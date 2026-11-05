import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';

import type { CredentialStatus, TropomiScreenResult } from '@/api/types';
import type { BBox } from '@/shared/geo';
import { rangeQuery, type DateRange } from '@/shared/time';

const readErrorMessage = async (response: Response) => {
  const body = (await response.json().catch(() => null)) as { detail?: unknown } | null;
  return typeof body?.detail === 'string' ? body.detail : `HTTP ${response.status}`;
};

const readJson = async <T>(response: Response): Promise<T> => {
  if (response.ok) {
    return (await response.json()) as T;
  }
  throw new Error(await readErrorMessage(response));
};

export const fetchCredentials = () =>
  fetch('/api/tropomi/credentials').then((response) => readJson<CredentialStatus>(response));

export const saveCredentials = (payload: {
  accessKeyId?: string;
  secretAccessKey?: string;
  fileText?: string;
}) =>
  fetch('/api/tropomi/credentials', {
    method: 'PUT',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload),
  }).then((response) => readJson<CredentialStatus>(response));

export const downloadTropomi = (range: DateRange) =>
  fetch('/api/tropomi/download', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(rangeQuery(range)),
  }).then((response) => readJson<{ message: string }>(response));

export const screenTropomi = (
  bbox: BBox,
  range: DateRange,
  thresholdDelta: number,
  minPixels: number,
) =>
  fetch('/api/tropomi/screen', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      west: bbox.west,
      south: bbox.south,
      east: bbox.east,
      north: bbox.north,
      ...rangeQuery(range),
      thresholdDelta,
      minPixels,
    }),
  }).then((response) => readJson<TropomiScreenResult>(response));

export const useCredentials = () =>
  useQuery({
    queryKey: ['tropomi', 'credentials'],
    queryFn: fetchCredentials,
  });

export const useSaveCredentials = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: saveCredentials,
    onSuccess: (status) => {
      queryClient.setQueryData(['tropomi', 'credentials'], status);
    },
  });
};

export default screenTropomi;
