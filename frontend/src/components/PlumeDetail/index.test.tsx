import { screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';

import type { Plume } from '@/api/types';
import { renderWithProviders } from '@/test/render';

import PlumeDetail from './index';

const plume: Plume = {
  id: 'tan20250625t181304c33s4001-A',
  source: 'carbonmapper',
  observedAt: '2025-06-25T18:13:04Z',
  lon: -101.828,
  lat: 32.264,
  emissionKgHr: 554.5,
  uncertaintyKgHr: 105.1,
  instrument: 'tan',
  sector: '1B2',
  imageUrl: 'https://catalog.carbonmapper.org/x/plume.png',
  attribution: 'Source: Carbon Mapper',
  details: { Platform: 'Planet Tanager 1' },
};

describe('PlumeDetail', () => {
  it('prompts for a selection', () => {
    renderWithProviders(<PlumeDetail plume={null} licenseUrl={undefined} />);
    expect(screen.getByText('Select a plume on the map or in the list.')).toBeInTheDocument();
  });

  it('shows the image, emission, and attribution', () => {
    renderWithProviders(<PlumeDetail plume={plume} licenseUrl="https://carbonmapper.org/terms" />);
    expect(screen.getByRole('img', { name: 'Plume image from Carbon Mapper' })).toHaveAttribute(
      'src',
      plume.imageUrl,
    );
    expect(screen.getByText('555 ± 105 kg/hr')).toBeInTheDocument();
    expect(screen.getByText('Source: Carbon Mapper')).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'Terms' })).toHaveAttribute(
      'href',
      'https://carbonmapper.org/terms',
    );
  });
});
