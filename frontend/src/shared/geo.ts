import type { Feature, Polygon } from 'geojson';

export type LngLat = { lng: number; lat: number };

export type BBox = { west: number; south: number; east: number; north: number };

export const bboxFromCorners = (first: LngLat, second: LngLat): BBox => ({
  west: Math.min(first.lng, second.lng),
  east: Math.max(first.lng, second.lng),
  south: Math.min(first.lat, second.lat),
  north: Math.max(first.lat, second.lat),
});

export const bboxToPolygon = (bbox: BBox): Feature<Polygon> => ({
  type: 'Feature',
  properties: {},
  geometry: {
    type: 'Polygon',
    coordinates: [
      [
        [bbox.west, bbox.south],
        [bbox.east, bbox.south],
        [bbox.east, bbox.north],
        [bbox.west, bbox.north],
        [bbox.west, bbox.south],
      ],
    ],
  },
});

export const formatDegrees = (value: number) => value.toFixed(2);

export const bboxQuery = (bbox: BBox) => ({
  west: String(bbox.west),
  south: String(bbox.south),
  east: String(bbox.east),
  north: String(bbox.north),
});
