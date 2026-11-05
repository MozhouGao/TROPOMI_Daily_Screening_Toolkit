import { DateTime } from 'luxon';

export type DateRange = { start: DateTime; end: DateTime };

export const QUICK_RANGE_DAYS = [30, 90, 180, 365] as const;

export const DEFAULT_RANGE_DAYS = 180;

const today = () => DateTime.utc().startOf('day');

export const quickRange = (days: number, end: DateTime = today()): DateRange => ({
  start: end.minus({ days: days - 1 }),
  end,
});

export const isQuickRange = (range: DateRange, days: number) => {
  const expected = quickRange(days);
  return range.start.hasSame(expected.start, 'day') && range.end.hasSame(expected.end, 'day');
};

export const rangeQuery = (range: DateRange) => ({
  start: range.start.toISODate() ?? '',
  end: range.end.toISODate() ?? '',
});

export const formatObserved = (iso: string) =>
  DateTime.fromISO(iso, { zone: 'utc' }).toFormat('yyyy-LL-dd HH:mm');
