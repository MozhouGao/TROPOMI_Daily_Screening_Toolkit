import { Typography } from 'antd';
import { useAtomValue } from 'jotai';

import usePlumeLayers, { type LayerQuery } from '@/api/plumes';
import CredentialManager from '@/components/CredentialManager';
import LayerPanel from '@/components/LayerPanel';
import MapView from '@/components/MapView';
import PlumeDetail from '@/components/PlumeDetail';
import PlumeList from '@/components/PlumeList';
import RegionControl from '@/components/RegionControl';
import TimeWindow from '@/components/TimeWindow';
import TropomiPanel from '@/components/TropomiPanel';
import { t } from '@/i18n';
import { mergeNewestFirst, plumeKey } from '@/shared/plumes';
import {
  bboxAtom,
  dateRangeAtom,
  enabledLayersAtom,
  selectedPlumeKeyAtom,
  tropomiResultAtom,
} from '@/state/atoms';

import styles from './index.module.scss';

const App = () => {
  const bbox = useAtomValue(bboxAtom);

  const range = useAtomValue(dateRangeAtom);

  const enabledLayers = useAtomValue(enabledLayersAtom);

  const selectedKey = useAtomValue(selectedPlumeKeyAtom);

  const tropomiResult = useAtomValue(tropomiResultAtom);

  const catalogQueries = usePlumeLayers(enabledLayers, bbox, range);

  const tropomiEnabled = enabledLayers.includes('tropomi');

  const tropomiQuery: LayerQuery = {
    source: 'tropomi',
    enabled: tropomiEnabled,
    data: tropomiEnabled && tropomiResult ? tropomiResult : undefined,
    isLoading: false,
    error: null,
  };

  const queries = [...catalogQueries, tropomiQuery];
  const plumes = mergeNewestFirst(
    queries.filter((query) => query.enabled).map((query) => query.data?.items ?? []),
  );
  const selectedPlume = plumes.find((plume) => plumeKey(plume) === selectedKey) ?? null;
  const selectedLicense = queries.find((query) => query.source === selectedPlume?.source)?.data
    ?.licenseUrl;
  const tropomiOverlay = tropomiEnabled ? (tropomiResult?.overlay ?? null) : null;

  return (
    <div className={styles.shell}>
      <header className={styles.header}>
        <Typography.Title level={4} className={styles.brand}>
          {t('app.title')}
        </Typography.Title>
        <Typography.Text type="secondary">{t('app.version')}</Typography.Text>
        <div className={styles.headerActions}>
          <CredentialManager />
        </div>
      </header>
      <aside className={styles.left}>
        <RegionControl />
        <TimeWindow />
        <TropomiPanel />
        <LayerPanel queries={queries} />
      </aside>
      <main className={styles.center}>
        <MapView queries={queries} selectedPlume={selectedPlume} tropomiOverlay={tropomiOverlay} />
      </main>
      <aside className={styles.right}>
        <PlumeList plumes={plumes} />
        <PlumeDetail plume={selectedPlume} licenseUrl={selectedLicense} />
      </aside>
    </div>
  );
};

export default App;
