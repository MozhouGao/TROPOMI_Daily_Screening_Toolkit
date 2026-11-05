import { useAtom, useSetAtom } from 'jotai';
import { useEffect, useState } from 'react';
import type { MapLayerMouseEvent } from 'react-map-gl/maplibre';

import { bboxFromCorners, type BBox, type LngLat } from '@/shared/geo';
import { bboxAtom, drawModeAtom, selectedPlumeKeyAtom } from '@/state/atoms';

const useRectangleDraw = () => {
  const [drawMode, setDrawMode] = useAtom(drawModeAtom);

  const setBbox = useSetAtom(bboxAtom);

  const setSelectedKey = useSetAtom(selectedPlumeKeyAtom);

  const [firstCorner, setFirstCorner] = useState<LngLat | null>(null);

  const [preview, setPreview] = useState<BBox | null>(null);

  useEffect(() => {
    function resetWhenDrawingStops() {
      if (drawMode === false) {
        setFirstCorner(null);
        setPreview(null);
      }
    }
    resetWhenDrawingStops();
  }, [drawMode]);

  const handleDrawClick = (event: MapLayerMouseEvent): boolean => {
    if (drawMode === false) {
      return false;
    }
    if (firstCorner === null) {
      setFirstCorner(event.lngLat);
      return true;
    }
    setBbox(bboxFromCorners(firstCorner, event.lngLat));
    setSelectedKey(null);
    setDrawMode(false);
    return true;
  };

  const handleDrawMove = (event: MapLayerMouseEvent) => {
    if (firstCorner) {
      setPreview(bboxFromCorners(firstCorner, event.lngLat));
    }
  };

  return { drawMode, preview, handleDrawClick, handleDrawMove };
};

export default useRectangleDraw;
