import { screen, waitFor } from '@testing-library/react';
import { afterEach, describe, expect, it, vi } from 'vitest';

import { bboxAtom } from '@/state/atoms';
import { createTestStore, renderWithProviders } from '@/test/render';

import NewsPanel from './index';

const collection = {
  items: [
    {
      id: 'permian-1',
      title: 'Permian methane leak found',
      url: 'https://example.org/permian-leak',
      source: 'Reuters',
      publishedAt: '2026-11-04T12:00:00Z',
      snippet: 'Operators in the Permian Basin reported a methane leak.',
      basins: ['Permian'],
      attribution: 'Headlines via Google News RSS. Link through to the publisher.',
    },
  ],
  totalCount: 1,
  truncated: false,
  attribution: 'Headlines via Google News RSS. Link through to the publisher.',
  licenseUrl: 'https://news.google.com',
};

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('NewsPanel', () => {
  it('asks for a region before loading headlines', () => {
    renderWithProviders(<NewsPanel />);
    expect(screen.getByText('Draw a region to load methane news.')).toBeInTheDocument();
  });

  it('lists headlines with basin tags and an external link', async () => {
    const store = createTestStore();
    store.set(bboxAtom, { west: -102.8, south: 31.6, east: -101.6, north: 32.4 });
    vi.stubGlobal(
      'fetch',
      vi.fn().mockResolvedValue({
        ok: true,
        json: async () => collection,
      }),
    );
    renderWithProviders(<NewsPanel />, store);
    expect(await screen.findByRole('link', { name: 'Permian methane leak found' })).toHaveAttribute(
      'href',
      'https://example.org/permian-leak',
    );
    const link = screen.getByRole('link', { name: 'Permian methane leak found' });
    expect(link).toHaveAttribute('target', '_blank');
    expect(link).toHaveAttribute('rel', 'noreferrer');
    expect(screen.getByText('Permian')).toBeInTheDocument();
    expect(screen.getByText(/Reuters/)).toBeInTheDocument();
    await waitFor(() => expect(fetch).toHaveBeenCalled());
  });
});
