import { Space, Tag, Typography } from 'antd';
import { useAtom } from 'jotai';
import type { DateTime } from 'luxon';

import { t } from '@/i18n';
import { DatePicker } from '@/shared/components/ui';
import { QUICK_RANGE_DAYS, isQuickRange, quickRange } from '@/shared/time';
import { dateRangeAtom } from '@/state/atoms';

import styles from './index.module.scss';

const { RangePicker } = DatePicker;

const TimeWindow = () => {
  const [range, setRange] = useAtom(dateRangeAtom);

  const applyPicker = (values: [DateTime | null, DateTime | null] | null) => {
    const [start, end] = values ?? [null, null];
    if (start && end) {
      setRange({ start: start.startOf('day'), end: end.startOf('day') });
    }
  };

  return (
    <section className={styles.panel} aria-labelledby="time-title">
      <Typography.Title level={5} id="time-title" className={styles.title}>
        {t('time.title')}
      </Typography.Title>
      <RangePicker
        value={[range.start, range.end]}
        onChange={applyPicker}
        allowClear={false}
        className={styles.picker}
      />
      <Space wrap className={styles.chips}>
        {QUICK_RANGE_DAYS.map((days) => (
          <Tag.CheckableTag
            key={days}
            checked={isQuickRange(range, days)}
            onChange={() => setRange(quickRange(days))}
          >
            {t('time.last', { days })}
          </Tag.CheckableTag>
        ))}
      </Space>
    </section>
  );
};

export default TimeWindow;
