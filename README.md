# Methane Data Screening Toolkit - Version 2.1

Browse methane plumes from several public sources on one map, and run the
existing TROPOMI patch screening over a drawn region and date range.

Layers:

- **TROPOMI screening** — Sentinel-5P XCH4 ultra-emitter patches (this toolkit's algorithm).
- **Carbon Mapper** — facility-scale CH4 plumes from the public catalog API.
- **UNEP IMEO MARS** — the public Eye on Methane plume export (CC BY-NC-SA 4.0).
- **GHGSat landfills** — Dogniaux et al. 2021–2022 landfill plumes (Zenodo 16641834).
- **US / Canada GHGRP** — latest self-reported facility CH4 (EPA Envirofacts, ECCC Data Mart).
- **Cusworth landfills 2024** — US landfill aircraft campaign with GHGRP IDs (default off).
- **Nesser landfills 2019** — TROPOMI inversion posteriors (default off).
- **Canada NIR 2018 grid** — Scarpelli 0.1° inventory cells ≥ 1 kg/hr (default off).
- **EDF MethaneSAT** — layer listing and terms only (default off). Raw L3/L4 stay on Earth Engine / Google Cloud.

Headlines: `GET /api/news` proxies Google News RSS for methane emissions/leaks in basins that overlap the drawn rectangle. The right sidebar lists publisher links only (not detections). XML is cached under `NEWS_data/` for about 30 minutes.

Carbon Mapper and MARS detections are independent of TROPOMI patch screening.
They are facility-scale observations with the provider's stated uncertainty.
Both live catalogs are licensed for **non-commercial** use. See
[Carbon Mapper terms](https://carbonmapper.org/terms) and the
[MARS data dictionary](https://methanedata.unep.org/dict-mars-plumes).

Archive layers (GHGSat, GHGRP, Cusworth, Nesser, Canada NIR) ignore the time
window because they are annual or campaign vintages. US GHGRP CH4 is converted
from reported tCO2e with AR4 GWP 25. Cache directories (`GHGSat_data/`,
`GHGRP_data/`, `Cusworth_data/`, `Nesser_data/`, `CanadaNIR_data/`,
`NEWS_data/`) stay local and are gitignored.

## Installation

Python 3.9 or newer, Node.js 20 or newer, and pnpm. Cartopy is easiest from conda-forge.

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
cd frontend
pnpm install
```

If `cartopy` fails to build with pip:

```bash
conda create -n tropomi python=3.10
conda activate tropomi
conda install -c conda-forge cartopy
pip install -r requirements.txt
```

## Run locally

From the repository root, start the API:

```bash
source .venv/bin/activate
uvicorn backend.main:app --reload --port 8000
```

In a second terminal:

```bash
cd frontend
pnpm dev
```

Open http://localhost:5173. The Vite dev server proxies `/api` to port 8000.

### Copernicus credentials (TROPOMI download)

Use **Copernicus credentials** in the header, or create `AWS_Keys.txt` in the
project folder:

```
access_key_id: YOUR_ACCESS_KEY
secret_access_key: YOUR_SECRET_KEY
```

Generate keys at https://eodata-s3keysmanager.dataspace.copernicus.eu/.
`AWS_Keys.txt` stays local and is gitignored.

## Quick guide

1. Draw a rectangle on the map.
2. Choose a time window (date range or a last-N-days chip).
3. Carbon Mapper and MARS load for that region and window. Toggle layers in the left panel.
4. Select a plume on the map or in the list to see time, instrument, emission,
   attribution, and (when Carbon Mapper publishes one) the plume image.
5. For TROPOMI: save Copernicus keys, click **Download Level-2**, set the
   enhancement threshold and minimum pixel count, then **Start screening**.
   The last day's XCH4 grid is drawn on the map; hotspot rows join the plume list.

MARS data is published 30 or more days after detection. The zip is cached under
`MARS_data/` and refreshed at most once a week. Carbon Mapper images are loaded
from the provider CDN and are not stored locally.

## Tests

```bash
python -m unittest discover -s tests -t .
python -m unittest discover -s backend/tests -t .
cd frontend && pnpm test && pnpm lint && pnpm build
```

## Internal default settings (TROPOMI)

The datasets are the Level-2 user products from the TROPOMI instrument onboard
the Sentinel-5P satellite. Atmospheric methane is the column averaged dry air
mole mixing ratio of methane (XCH4) in ppb. For each file,
`methane_mixing_ratio_bias_corrected` is remapped on a 0.05° by 0.05° grid.
Only quality-assured points (`qa_value` greater than 0.5) are used.

The plume screening algorithm originally proposed by Lauvaux et al. (2022) is
simplified to highlight regions with high probability of detecting methane
ultra emitters. Screening is patch-based. Each patch is 0.5° by 0.5° (roughly
55 km by 55 km at low latitudes). This toolkit is not designed to pinpoint
facility- or component-level sources, or to screen small leaks (generally
<25 tons/hour).

### Background XCH4

Background is computed from the skewness of the patch,
`((XCH4,mean - XCH4,median) / XCH4,std)`. Strongly skewed patches
(skewness > 0.3) use the median; otherwise a robust estimator with a typical
clipping parameter of 2.5 is used (Lauvaux et al., 2022).

### User defined parameters

Two parameters control how closely an anomaly must resemble a plume: mean
XCH4 enhancement over positive pixels (ΔXCH4,thr, default 15 ppb) and the
number of valid positive pixels (n, default 1). Larger values are more
conservative.

Screening maps and CSVs are still written to `./assets`.

## Feedback

zhenyu.xing2@ucalgary.ca and mozhou.gao@ucalgary.ca.

## References

Lauvaux, T., Giron, C., Mazzolini, M., d’Aspremont, A., Duren, R., Cusworth, D.,
Shindell, D. and Ciais, P., 2022. Global assessment of oil and gas methane
ultra-emitters. Science, 375(6580), pp.557-561.
