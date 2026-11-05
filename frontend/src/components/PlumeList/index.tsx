import { Table, Typography } from 'antd';
import type { ColumnsType } from 'antd/es/table';
import { useAtom } from 'jotai';

import type { Plume } from '@/api/types';
import { t } from '@/i18n';
import { layerById } from '@/layers/registry';
import { formatEmission, plumeKey } from '@/shared/plumes';
import { formatObserved } from '@/shared/time';
import { selectedPlumeKeyAtom } from '@/state/atoms';

import styles from './index.module.scss';

type PlumeListProps = { plumes: Plume[] };

const columns: ColumnsType<Plume> = [
  {
    title: t('list.column.time'),
    dataIndex: 'observedAt',
    render: (value: string) => formatObserved(value),
    width: 140,
  },
  {
    title: t('list.column.layer'),
    dataIndex: 'source',
    render: (_: string, plume: Plume) => t(layerById(plume.source).labelKey),
  },
  {
    title: t('list.column.instrument'),
    dataIndex: 'instrument',
    render: (value: string | null) => value ?? '',
  },
  {
    title: t('list.column.emission'),
    key: 'emission',
    align: 'right',
    render: (_: unknown, plume: Plume) => formatEmission(plume),
  },
];

const PlumeList = ({ plumes }: PlumeListProps) => {
  const [selectedKey, setSelectedKey] = useAtom(selectedPlumeKeyAtom);

  return (
    <section className={styles.panel} aria-labelledby="list-title">
      <Typography.Title level={5} id="list-title" className={styles.title}>
        {t('list.title')} ({plumes.length})
      </Typography.Title>
      <Table<Plume>
        size="small"
        rowKey={plumeKey}
        columns={columns}
        dataSource={plumes}
        pagination={{ pageSize: 8, size: 'small', hideOnSinglePage: true }}
        locale={{ emptyText: t('list.empty') }}
        rowClassName={(plume) =>
          (plumeKey(plume) === selectedKey ? styles.selected : styles.row) ?? ''
        }
        onRow={(plume) => ({ onClick: () => setSelectedKey(plumeKey(plume)) })}
      />
    </section>
  );
};

export default PlumeList;
