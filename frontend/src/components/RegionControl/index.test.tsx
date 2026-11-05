import { fireEvent, screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';

import { bboxAtom } from '@/state/atoms';
import { createTestStore, renderWithProviders } from '@/test/render';

import RegionControl from './index';

describe('RegionControl', () => {
  it('asks the user to draw a rectangle before any request', () => {
    renderWithProviders(<RegionControl />);
    expect(screen.getByText('No region yet. Draw a rectangle to load plumes.')).toBeInTheDocument();
  });

  it('shows the selected extent', () => {
    const store = createTestStore();
    store.set(bboxAtom, { west: 10, south: 20, east: 12, north: 23 });
    renderWithProviders(<RegionControl />, store);
    expect(screen.getByText('lon [10.00, 12.00], lat [20.00, 23.00]')).toBeInTheDocument();
  });

  it('clears the region', () => {
    const store = createTestStore();
    store.set(bboxAtom, { west: 10, south: 20, east: 12, north: 23 });
    renderWithProviders(<RegionControl />, store);
    fireEvent.click(screen.getByRole('button', { name: 'Clear' }));
    expect(store.get(bboxAtom)).toBeNull();
    expect(screen.getByText('No region yet. Draw a rectangle to load plumes.')).toBeInTheDocument();
  });
});
