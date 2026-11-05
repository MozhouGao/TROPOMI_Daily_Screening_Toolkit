import { screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';

import type { LayerQuery } from '@/api/plumes';
import { renderWithProviders } from '@/test/render';

import LayerPanel from './index';

const catalogQueries: LayerQuery[] = [
  {
    source: 'carbonmapper',
    enabled: true,
    isLoading: false,
    error: null,
    data: {
      source: 'carbonmapper',
      items: [],
      totalCount: 237,
      truncated: true,
      attribution: 'Source: Carbon Mapper',
      licenseUrl: 'https://carbonmapper.org/terms',
    },
  },
  {
    source: 'mars',
    enabled: true,
    isLoading: false,
    error: null,
    data: {
      source: 'mars',
      items: [],
      totalCount: 1,
      truncated: false,
      attribution: 'Source: UNEP IMEO Eye on Methane (CC BY-NC-SA 4.0)',
      licenseUrl: 'https://methanedata.unep.org/dict-mars-plumes',
    },
  },
];

describe('LayerPanel', () => {
  it('lists every registered layer with attribution', () => {
    renderWithProviders(<LayerPanel queries={catalogQueries} />);
    expect(screen.getByText('Carbon Mapper')).toBeInTheDocument();
    expect(screen.getByText('UNEP IMEO MARS')).toBeInTheDocument();
    expect(screen.getByText('TROPOMI screening')).toBeInTheDocument();
    expect(screen.getByText('GHGSat landfills')).toBeInTheDocument();
    expect(screen.getByText('US GHGRP facilities')).toBeInTheDocument();
    expect(screen.getByText('Canada GHGRP facilities')).toBeInTheDocument();
    expect(screen.getByText('EDF MethaneSAT')).toBeInTheDocument();
    expect(
      screen.getByText('Showing 0 of 237. Narrow the region or time window.'),
    ).toBeInTheDocument();
    expect(screen.getByText('Source: Carbon Mapper')).toBeInTheDocument();
  });
});
