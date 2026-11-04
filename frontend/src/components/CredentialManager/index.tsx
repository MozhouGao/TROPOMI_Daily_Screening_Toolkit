import { Button, Input, Modal, Typography } from 'antd';
import { useEffect, useState } from 'react';

import { useCredentials, useSaveCredentials } from '@/api/tropomi';
import { t } from '@/i18n';
import Loading from '@/shared/components/ui/Loading';

import styles from './index.module.scss';

const CREDENTIALS_URL = 'https://eodata-s3keysmanager.dataspace.copernicus.eu/';

const CredentialManager = () => {
  const [open, setOpen] = useState(false);

  const [accessKeyId, setAccessKeyId] = useState('');

  const [secretAccessKey, setSecretAccessKey] = useState('');

  const [fileText, setFileText] = useState('');

  const credentials = useCredentials();

  const save = useSaveCredentials();

  useEffect(() => {
    function hydrateFromServer() {
      const status = credentials.data;
      if (status) {
        setAccessKeyId(status.accessKeyId);
        setFileText(status.fileText);
      }
    }
    hydrateFromServer();
  }, [credentials.data]);

  const configured = credentials.data?.configured ?? false;

  return (
    <>
      <div>
        <Button size="small" onClick={() => setOpen(true)}>
          {t('credentials.open')}
        </Button>
        <div className={styles.banner}>
          {configured ? t('credentials.configured') : t('credentials.missing')}
        </div>
      </div>
      <Modal
        title={t('credentials.title')}
        open={open}
        onCancel={() => setOpen(false)}
        footer={null}
      >
        <Typography.Paragraph className={styles.help}>
          {t('credentials.help')}{' '}
          <a href={CREDENTIALS_URL} target="_blank" rel="noreferrer">
            {CREDENTIALS_URL}
          </a>
        </Typography.Paragraph>
        {credentials.isLoading ? <Loading /> : null}
        <label className={styles.field} htmlFor="cred-access">
          {t('credentials.access')}
        </label>
        <Input
          id="cred-access"
          value={accessKeyId}
          autoComplete="off"
          className={styles.field}
          onChange={(event) => setAccessKeyId(event.target.value)}
        />
        <label className={styles.field} htmlFor="cred-secret">
          {t('credentials.secret')}
        </label>
        <Input.Password
          id="cred-secret"
          value={secretAccessKey}
          autoComplete="off"
          className={styles.field}
          onChange={(event) => setSecretAccessKey(event.target.value)}
        />
        <label className={styles.field} htmlFor="cred-file">
          {t('credentials.file')}
        </label>
        <Input.TextArea
          id="cred-file"
          value={fileText}
          rows={5}
          className={`${styles.field} ${styles.file}`}
          onChange={(event) => setFileText(event.target.value)}
        />
        <Button type="primary" onClick={() => save.mutate({ accessKeyId, secretAccessKey })}>
          {t('credentials.save')}
        </Button>{' '}
        <Button onClick={() => save.mutate({ fileText })}>{t('credentials.saveFile')}</Button>{' '}
        <Button onClick={() => setOpen(false)}>{t('credentials.close')}</Button>
        {save.data || save.error ? (
          <div className={styles.status}>
            {save.error ? save.error.message : t('credentials.configured')}
          </div>
        ) : null}
      </Modal>
    </>
  );
};

export default CredentialManager;
