import { Checkbox, Typography } from 'antd';
import { useAtom } from 'jotai';

import type { LayerQuery } from '@/api/plumes';
import type { PlumeSource } from '@/api/types';
import { t } from '@/i18n';
import { countKeyForKind, LAYER_REGISTRY, type LayerDefinition } from '@/layers/registry';
import { Loading } from '@/shared/components/ui';
import { enabledLayersAtom } from '@/state/atoms';

import styles from './index.module.scss';

type LayerPanelProps = { queries: LayerQuery[] };

type LayerRowProps = {
  layer: LayerDefinition;
  query: LayerQuery | undefined;
  onToggle: (id: PlumeSource, enabled: boolean) => void;
};

const LayerStatus = ({
  query,
  layer,
}: {
  query: LayerQuery | undefined;
  layer: LayerDefinition;
}) => {
  if (query?.isLoading) {
    return <Loading label={t('layers.loading')} />;
  }
  if (query?.error) {
    return (
      <Typography.Text type="danger" className={styles.status}>
        {t('layers.error', { message: query.error.message })}
      </Typography.Text>
    );
  }
  if (query?.data) {
    const { items, totalCount, truncated } = query.data;
    return (
      <Typography.Text className={styles.status}>
        {truncated
          ? t('layers.truncated', { shown: items.length, total: totalCount })
          : t(countKeyForKind(layer.kind), { count: items.length })}
      </Typography.Text>
    );
  }
  return null;
};

const LayerRow = ({ layer, query, onToggle }: LayerRowProps) => (
  <li className={styles.row}>
    <Checkbox
      checked={query?.enabled ?? layer.defaultEnabled}
      onChange={(event) => onToggle(layer.id, event.target.checked)}
    >
      <span className={`${styles.swatch} ${styles[`swatch_${layer.id}`] ?? ''}`} />
      <span className={styles.label}>{t(layer.labelKey)}</span>
    </Checkbox>
    <div className={styles.meta}>
      <LayerStatus query={query} layer={layer} />
      <Typography.Text type="secondary" className={styles.note}>
        {t(layer.freshnessKey)}
      </Typography.Text>
      {query?.data && (
        <Typography.Text type="secondary" className={styles.note}>
          {query.data.attribution}{' '}
          <a href={query.data.licenseUrl} target="_blank" rel="noreferrer">
            {t('layers.license')}
          </a>
        </Typography.Text>
      )}
    </div>
  </li>
);

const LayerPanel = ({ queries }: LayerPanelProps) => {
  const [enabledLayers, setEnabledLayers] = useAtom(enabledLayersAtom);

  const toggleLayer = (id: PlumeSource, enabled: boolean) =>
    setEnabledLayers(
      enabled ? [...enabledLayers, id] : enabledLayers.filter((layerId) => layerId !== id),
    );

  return (
    <section className={styles.panel} aria-labelledby="layers-title">
      <Typography.Title level={5} id="layers-title" className={styles.title}>
        {t('layers.title')}
      </Typography.Title>
      <ul className={styles.list}>
        {LAYER_REGISTRY.map((layer) => (
          <LayerRow
            key={layer.id}
            layer={layer}
            query={queries.find((query) => query.source === layer.id)}
            onToggle={toggleLayer}
          />
        ))}
      </ul>
    </section>
  );
};

export default LayerPanel;
