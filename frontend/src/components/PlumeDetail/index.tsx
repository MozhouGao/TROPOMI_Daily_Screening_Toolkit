import { Descriptions, Typography } from 'antd';

import type { Plume } from '@/api/types';
import { t } from '@/i18n';
import { layerById } from '@/layers/registry';
import { formatDegrees } from '@/shared/geo';
import { formatEmission } from '@/shared/plumes';
import { formatObserved } from '@/shared/time';

import styles from './index.module.scss';

type PlumeDetailProps = { plume: Plume | null; licenseUrl: string | undefined };

const PlumeDetail = ({ plume, licenseUrl }: PlumeDetailProps) => {
  if (plume === null) {
    return (
      <section className={styles.panel} aria-labelledby="detail-title">
        <Typography.Title level={5} id="detail-title" className={styles.title}>
          {t('detail.title')}
        </Typography.Title>
        <Typography.Text type="secondary">{t('detail.empty')}</Typography.Text>
      </section>
    );
  }

  const layer = layerById(plume.source);
  const layerLabel = t(layer.labelKey);
  const emission = formatEmission(plume);
  const showImage = layer.kind === 'catalog' || layer.kind === 'screening';

  return (
    <section className={styles.panel} aria-labelledby="detail-title">
      <Typography.Title level={5} id="detail-title" className={styles.title}>
        {t('detail.title')}: {layerLabel}
      </Typography.Title>
      {showImage &&
        (plume.imageUrl ? (
          <img
            className={styles.image}
            src={plume.imageUrl}
            alt={t('detail.imageAlt', { layer: layerLabel })}
          />
        ) : (
          <Typography.Text type="secondary" className={styles.noImage}>
            {t('detail.noImage')}
          </Typography.Text>
        ))}
      <Descriptions size="small" column={1} className={styles.fields}>
        <Descriptions.Item label={t('detail.observed')}>
          {formatObserved(plume.observedAt)} UTC
        </Descriptions.Item>
        <Descriptions.Item label={t('detail.location')}>
          {formatDegrees(plume.lat)}, {formatDegrees(plume.lon)}
        </Descriptions.Item>
        <Descriptions.Item label={t('detail.emission')}>
          {emission ? `${emission} kg/hr` : t('common.unknown')}
        </Descriptions.Item>
        <Descriptions.Item label={t('detail.instrument')}>
          {plume.instrument ?? t('common.unknown')}
        </Descriptions.Item>
        <Descriptions.Item label={t('detail.sector')}>
          {plume.sector ?? t('common.unknown')}
        </Descriptions.Item>
        {Object.entries(plume.details).map(([label, value]) => (
          <Descriptions.Item key={label} label={label}>
            {value}
          </Descriptions.Item>
        ))}
        <Descriptions.Item label={t('detail.id')}>
          <Typography.Text code copyable>
            {plume.id}
          </Typography.Text>
        </Descriptions.Item>
      </Descriptions>
      <Typography.Paragraph className={styles.attribution}>
        <strong>{plume.attribution}</strong>
        {licenseUrl && (
          <>
            {' '}
            <a href={licenseUrl} target="_blank" rel="noreferrer">
              {t('layers.license')}
            </a>
          </>
        )}
        <br />
        {t(layer.disclaimerKey)}
      </Typography.Paragraph>
    </section>
  );
};

export default PlumeDetail;
