import { atom } from 'jotai';

import type { PlumeSource, TropomiScreenResult } from '@/api/types';
import { defaultEnabledLayers } from '@/layers/registry';
import type { BBox } from '@/shared/geo';
import { DEFAULT_RANGE_DAYS, quickRange } from '@/shared/time';

export const bboxAtom = atom<BBox | null>(null);

export const dateRangeAtom = atom(quickRange(DEFAULT_RANGE_DAYS));

export const enabledLayersAtom = atom<PlumeSource[]>(defaultEnabledLayers());

export const selectedPlumeKeyAtom = atom<string | null>(null);

export const drawModeAtom = atom(false);

export const tropomiResultAtom = atom<TropomiScreenResult | null>(null);
