import { Button, InputNumber, Typography } from 'antd';
import { useAtom, useAtomValue } from 'jotai';
import { useState } from 'react';

import { downloadTropomi, screenTropomi } from '@/api/tropomi';
import { t } from '@/i18n';
import Loading from '@/shared/components/ui/Loading';
import { bboxAtom, dateRangeAtom, tropomiResultAtom } from '@/state/atoms';

import styles from './index.module.scss';

const TropomiPanel = () => {
  const bbox = useAtomValue(bboxAtom);

  const range = useAtomValue(dateRangeAtom);

  const [, setResult] = useAtom(tropomiResultAtom);

  const [thresholdDelta, setThresholdDelta] = useState(15);

  const [minPixels, setMinPixels] = useState(1);

  const [busy, setBusy] = useState<'download' | 'screen' | null>(null);

  const [message, setMessage] = useState('');

  const runDownload = async () => {
    setBusy('download');
    try {
      const body = await downloadTropomi(range);
      setMessage(body.message);
    } catch (error) {
      setMessage(error instanceof Error ? error.message : String(error));
    } finally {
      setBusy(null);
    }
  };

  const runScreen = async () => {
    if (bbox === null) {
      setMessage(t('tropomi.needRegion'));
      return;
    }
    setBusy('screen');
    try {
      const body = await screenTropomi(bbox, range, thresholdDelta, minPixels);
      setResult(body);
      setMessage(body.message);
    } catch (error) {
      setMessage(error instanceof Error ? error.message : String(error));
    } finally {
      setBusy(null);
    }
  };

  return (
    <section className={styles.panel} aria-labelledby="tropomi-title">
      <Typography.Title level={5} id="tropomi-title" className={styles.title}>
        {t('tropomi.title')}
      </Typography.Title>
      <Button onClick={() => void runDownload()} disabled={busy !== null} className={styles.field}>
        {t('tropomi.download')}
      </Button>
      <label htmlFor="tropomi-threshold">{t('tropomi.threshold')}</label>
      <InputNumber
        id="tropomi-threshold"
        className={styles.field}
        min={1}
        value={thresholdDelta}
        onChange={(value) => setThresholdDelta(value ?? 15)}
      />
      <label htmlFor="tropomi-min-pixels">{t('tropomi.minPixels')}</label>
      <InputNumber
        id="tropomi-min-pixels"
        className={styles.field}
        min={1}
        value={minPixels}
        onChange={(value) => setMinPixels(value ?? 1)}
      />
      <Button
        type="primary"
        onClick={() => void runScreen()}
        disabled={busy !== null}
        className={styles.field}
      >
        {t('tropomi.screen')}
      </Button>
      {busy ? <Loading /> : null}
      {message ? (
        <Typography.Paragraph className={styles.message}>{message}</Typography.Paragraph>
      ) : null}
    </section>
  );
};

export default TropomiPanel;
