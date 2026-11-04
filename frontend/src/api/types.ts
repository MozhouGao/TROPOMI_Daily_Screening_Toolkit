export type PlumeSource =
  | 'carbonmapper'
  | 'mars'
  | 'tropomi'
  | 'ghgsat'
  | 'cusworth'
  | 'ghgrp_us'
  | 'ghgrp_ca'
  | 'nesser'
  | 'canada_nir';

export type Plume = {
  id: string;
  source: PlumeSource;
  observedAt: string;
  lon: number;
  lat: number;
  emissionKgHr: number | null;
  uncertaintyKgHr: number | null;
  instrument: string | null;
  sector: string | null;
  imageUrl: string | null;
  attribution: string;
  details: Record<string, string>;
};

export type PlumeCollection = {
  source: PlumeSource;
  items: Plume[];
  totalCount: number;
  truncated: boolean;
  attribution: string;
  licenseUrl: string;
};

export type GridOverlay = {
  imageUrl: string;
  west: number;
  south: number;
  east: number;
  north: number;
};

export type TropomiScreenResult = PlumeCollection & {
  overlay: GridOverlay | null;
  message: string;
  figurePath: string | null;
  screenedDays: number;
  plumeDays: number;
};

export type CredentialStatus = {
  configured: boolean;
  accessKeyId: string;
  fileText: string;
  example: string;
};
