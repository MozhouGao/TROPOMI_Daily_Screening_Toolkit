import en from './locales/en.json';

export type TranslationKey = keyof typeof en;

type Variables = Record<string, string | number>;

const interpolate = (template: string, variables: Variables) =>
  template.replace(/\{(\w+)\}/g, (_, name: string) => String(variables[name] ?? `{${name}}`));

export const t = (key: TranslationKey, variables: Variables = {}) =>
  interpolate(en[key], variables);
