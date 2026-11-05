import { LoadingOutlined } from '@ant-design/icons';

import { t } from '@/i18n';

import styles from './index.module.scss';

type LoadingProps = { label?: string };

const Loading = ({ label = t('common.loading') }: LoadingProps) => (
  <span className={styles.loading} role="status" aria-live="polite">
    <LoadingOutlined spin />
    <span>{label}</span>
  </span>
);

export default Loading;
