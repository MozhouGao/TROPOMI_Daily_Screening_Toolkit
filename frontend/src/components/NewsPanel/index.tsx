import { Typography } from 'antd';
import { useAtomValue } from 'jotai';

import useNews from '@/api/news';
import { t } from '@/i18n';
import Loading from '@/shared/components/ui/Loading';
import { formatObserved } from '@/shared/time';
import { bboxAtom } from '@/state/atoms';

import styles from './index.module.scss';

const NewsPanel = () => {
  const bbox = useAtomValue(bboxAtom);

  const query = useNews(bbox);

  return (
    <section className={styles.panel} aria-labelledby="news-title">
      <Typography.Title level={5} id="news-title" className={styles.title}>
        {t('news.title')}
      </Typography.Title>
      {bbox === null && <Typography.Text type="secondary">{t('news.empty')}</Typography.Text>}
      {bbox !== null && query.isFetching && query.data === undefined && (
        <Loading label={t('common.loading')} />
      )}
      {bbox !== null && query.error && (
        <Typography.Text type="danger">
          {t('news.error', { message: query.error.message })}
        </Typography.Text>
      )}
      {bbox !== null && query.data && query.data.items.length === 0 && (
        <Typography.Text type="secondary">{t('news.none')}</Typography.Text>
      )}
      {query.data && query.data.items.length > 0 && (
        <>
          <ul className={styles.list}>
            {query.data.items.map((article) => (
              <li key={article.id} className={styles.item}>
                <a className={styles.headline} href={article.url} target="_blank" rel="noreferrer">
                  {article.title}
                </a>
                <span className={styles.meta}>
                  {article.source}
                  {article.publishedAt ? ` · ${formatObserved(article.publishedAt)}` : ''}
                </span>
                {article.basins.length > 0 && (
                  <span className={styles.tags}>
                    {article.basins.map((basin) => (
                      <span key={basin} className={styles.tag}>
                        {basin}
                      </span>
                    ))}
                  </span>
                )}
              </li>
            ))}
          </ul>
          <Typography.Paragraph type="secondary" className={styles.attribution}>
            {t('news.attribution')}
          </Typography.Paragraph>
        </>
      )}
    </section>
  );
};

export default NewsPanel;
