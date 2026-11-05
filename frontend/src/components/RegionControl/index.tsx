import { Button, Space, Typography } from 'antd';
import { useAtom, useSetAtom } from 'jotai';

import { t } from '@/i18n';
import { formatDegrees } from '@/shared/geo';
import { bboxAtom, drawModeAtom, selectedPlumeKeyAtom } from '@/state/atoms';

import styles from './index.module.scss';

const RegionControl = () => {
  const [bbox, setBbox] = useAtom(bboxAtom);

  const [drawMode, setDrawMode] = useAtom(drawModeAtom);

  const setSelectedKey = useSetAtom(selectedPlumeKeyAtom);

  const clearRegion = () => {
    setBbox(null);
    setSelectedKey(null);
    setDrawMode(false);
  };

  return (
    <section className={styles.panel} aria-labelledby="region-title">
      <Typography.Title level={5} id="region-title" className={styles.title}>
        {t('region.title')}
      </Typography.Title>
      <Space wrap>
        <Button type={drawMode ? 'default' : 'primary'} onClick={() => setDrawMode(!drawMode)}>
          {drawMode ? t('region.cancel') : t('region.draw')}
        </Button>
        <Button onClick={clearRegion} disabled={bbox === null}>
          {t('region.clear')}
        </Button>
      </Space>
      <Typography.Paragraph className={styles.extent}>
        {drawMode && t('region.drawing')}
        {drawMode === false && bbox === null && t('region.none')}
        {drawMode === false &&
          bbox &&
          t('region.extent', {
            west: formatDegrees(bbox.west),
            east: formatDegrees(bbox.east),
            south: formatDegrees(bbox.south),
            north: formatDegrees(bbox.north),
          })}
      </Typography.Paragraph>
    </section>
  );
};

export default RegionControl;
