import { useAtom, useAtomValue } from 'jotai';
import type { ExpressionSpecification, StyleSpecification } from 'maplibre-gl';
import 'maplibre-gl/dist/maplibre-gl.css';
import { useEffect, useRef } from 'react';
import Map, {
  Layer,
  NavigationControl,
  Source,
  type MapLayerMouseEvent,
  type MapRef,
} from 'react-map-gl/maplibre';

import type { LayerQuery } from '@/api/plumes';
import type { GridOverlay, Plume } from '@/api/types';
import { LAYER_REGISTRY, layerById } from '@/layers/registry';
import { bboxToPolygon } from '@/shared/geo';
import { plumeKey, toFeatureCollection } from '@/shared/plumes';
import { bboxAtom, selectedPlumeKeyAtom } from '@/state/atoms';

import styles from './index.module.scss';
import useRectangleDraw from './useRectangleDraw';

const OSM_STYLE: StyleSpecification = {
  version: 8,
  sources: {
    osm: {
      type: 'raster',
      tiles: ['https://tile.openstreetmap.org/{z}/{x}/{y}.png'],
      tileSize: 256,
      attribution: '© OpenStreetMap contributors',
    },
  },
  layers: [{ id: 'osm', type: 'raster', source: 'osm' }],
};

const INITIAL_VIEW = { longitude: -114.07, latitude: 51.04, zoom: 4 };
const SELECTED_ZOOM = 10;
const plumeLayerId = (source: string) => `plumes-${source}`;

type MapViewProps = {
  queries: LayerQuery[];
  selectedPlume: Plume | null;
  tropomiOverlay: GridOverlay | null;
};

const circleRadius: ExpressionSpecification = [
  'interpolate',
  ['linear'],
  ['get', 'emission'],
  0,
  4,
  1000,
  7,
  10000,
  12,
];

const MapView = ({ queries, selectedPlume, tropomiOverlay }: MapViewProps) => {
  const mapRef = useRef<MapRef>(null);

  const bbox = useAtomValue(bboxAtom);

  const [selectedKey, setSelectedKey] = useAtom(selectedPlumeKeyAtom);

  const { drawMode, preview, handleDrawClick, handleDrawMove } = useRectangleDraw();

  useEffect(() => {
    function flyToSelectedPlume() {
      const map = mapRef.current;
      if (map && selectedPlume) {
        const zoom = Math.max(map.getZoom(), SELECTED_ZOOM);
        map.flyTo({ center: [selectedPlume.lon, selectedPlume.lat], zoom, duration: 800 });
      }
    }
    flyToSelectedPlume();
  }, [selectedPlume]);

  const handleClick = (event: MapLayerMouseEvent) => {
    if (handleDrawClick(event)) {
      return;
    }
    const key = event.features?.[0]?.properties?.key;
    setSelectedKey(typeof key === 'string' ? key : null);
  };

  const visibleQueries = queries.filter((query) => query.enabled && query.data);
  const region = preview ?? bbox;

  return (
    <div className={styles.map}>
      <Map
        ref={mapRef}
        initialViewState={INITIAL_VIEW}
        mapStyle={OSM_STYLE}
        cursor={drawMode ? 'crosshair' : 'auto'}
        interactiveLayerIds={LAYER_REGISTRY.map((layer) => plumeLayerId(layer.id))}
        onClick={handleClick}
        onMouseMove={handleDrawMove}
      >
        <NavigationControl position="top-left" />
        {tropomiOverlay && (
          <Source
            id="tropomi-grid"
            type="image"
            url={tropomiOverlay.imageUrl}
            coordinates={[
              [tropomiOverlay.west, tropomiOverlay.north],
              [tropomiOverlay.east, tropomiOverlay.north],
              [tropomiOverlay.east, tropomiOverlay.south],
              [tropomiOverlay.west, tropomiOverlay.south],
            ]}
          >
            <Layer id="tropomi-grid-raster" type="raster" paint={{ 'raster-opacity': 0.65 }} />
          </Source>
        )}
        {region && (
          <Source id="region" type="geojson" data={bboxToPolygon(region)}>
            <Layer
              id="region-fill"
              type="fill"
              paint={{ 'fill-color': '#f59e0b', 'fill-opacity': 0.08 }}
            />
            <Layer
              id="region-line"
              type="line"
              paint={{ 'line-color': '#f59e0b', 'line-width': 2 }}
            />
          </Source>
        )}
        {visibleQueries.map((query) => (
          <Source
            key={query.source}
            id={plumeLayerId(query.source)}
            type="geojson"
            data={toFeatureCollection(query.data?.items ?? [])}
          >
            <Layer
              id={plumeLayerId(query.source)}
              type="circle"
              paint={{
                'circle-color': layerById(query.source).color,
                'circle-radius': circleRadius,
                'circle-opacity': 0.75,
                'circle-stroke-color': '#ffffff',
                'circle-stroke-width': 1,
              }}
            />
            <Layer
              id={`${plumeLayerId(query.source)}-selected`}
              type="circle"
              filter={['==', ['get', 'key'], selectedKey ?? '']}
              paint={{
                'circle-color': 'transparent',
                'circle-radius': 14,
                'circle-stroke-color': '#111827',
                'circle-stroke-width': 3,
              }}
            />
          </Source>
        ))}
      </Map>
      {selectedPlume && <div className={styles.badge}>{plumeKey(selectedPlume)}</div>}
    </div>
  );
};

export default MapView;
