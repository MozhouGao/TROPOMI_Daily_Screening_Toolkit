import { fireEvent, screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';

import { createTestStore, renderWithProviders } from '@/test/render';

import TropomiPanel from './index';

describe('TropomiPanel', () => {
  it('asks for a region before screening', () => {
    renderWithProviders(<TropomiPanel />, createTestStore());
    fireEvent.click(screen.getByRole('button', { name: 'Start screening' }));
    expect(screen.getByText('Please define the region.')).toBeInTheDocument();
  });
});
