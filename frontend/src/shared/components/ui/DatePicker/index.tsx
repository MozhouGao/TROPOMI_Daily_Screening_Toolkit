import generatePicker from 'antd/es/date-picker/generatePicker';
import type { DateTime } from 'luxon';
import luxonGenerateConfig from 'rc-picker/es/generate/luxon';

const DatePicker = generatePicker<DateTime>(luxonGenerateConfig);

export default DatePicker;
