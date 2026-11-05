import { useQuery } from '@tanstack/react-query';

import type { NewsCollection } from '@/api/types';
import { bboxQuery, type BBox } from '@/shared/geo';

const FIVE_MINUTES = 5 * 60 * 1000;

export const newsUrl = (bbox: BBox) => `/api/news?${new URLSearchParams(bboxQuery(bbox))}`;

const readErrorMessage = async (response: Response) => {
  const body = (await response.json().catch(() => null)) as { detail?: unknown } | null;
  return typeof body?.detail === 'string' ? body.detail : `HTTP ${response.status}`;
};

export const fetchNews = async (bbox: BBox): Promise<NewsCollection> => {
  const response = await fetch(newsUrl(bbox));
  if (response.ok) {
    return (await response.json()) as NewsCollection;
  }
  throw new Error(await readErrorMessage(response));
};

const useNews = (bbox: BBox | null) =>
  useQuery({
    queryKey: ['news', bbox],
    queryFn: () => fetchNews(bbox as BBox),
    enabled: bbox !== null,
    staleTime: FIVE_MINUTES,
    retry: 1,
  });

export default useNews;
