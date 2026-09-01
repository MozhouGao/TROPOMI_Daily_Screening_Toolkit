import os
import re
from datetime import date, datetime, timedelta
from math import sqrt

import numpy as np
import pandas as pd
from matplotlib.path import Path

GRID_RESOLUTION = 0.05
CH4_FILL_VALUE = 9.96921e36
DEFAULT_PRODUCT_TYPES = ("OFFL", "NRTI")


def iter_screening_dates(start_date, end_date=None):
    """Yield each calendar date from start through end, inclusive."""
    start = _as_date(start_date)
    end = start if end_date is None else _as_date(end_date)
    if end < start:
        start, end = end, start
    current = start
    while current <= end:
        yield current
        current += timedelta(days=1)


def _as_date(value):
    if isinstance(value, datetime):
        return value.date()
    if isinstance(value, date):
        return value
    return datetime.strptime(str(value)[:10], "%Y-%m-%d").date()


def bbox_from_coordinates(coordinates):
    """Return (minlat, maxlat, minlon, maxlon) from a GeoJSON ring."""
    lons = [float(point[0]) for point in coordinates]
    lats = [float(point[1]) for point in coordinates]
    if not lons or not lats:
        raise ValueError("Polygon has no coordinates.")
    return min(lats), max(lats), min(lons), max(lons)


def tropomi_sensing_start_date(filename):
    """Return YYYYMMDD of the first sensing timestamp in an S5P filename."""
    match = re.search(r"_(\d{8})T\d{6}_", os.path.basename(filename))
    return match.group(1) if match else None


def is_ch4_l2_nc(filename, product_types=None):
    name = os.path.basename(filename)
    if not name.endswith(".nc") or "L2__CH4" not in name:
        return False
    if product_types is None:
        return name.startswith("S5P_")
    return any(f"S5P_{ptype}_" in name for ptype in product_types)


def local_download_path(target_dir, filename):
    return os.path.join(target_dir, os.path.basename(filename))


def read_aws_keys(keys_path=None):
    keys_path = keys_path or os.path.join(os.getcwd(), "AWS_Keys.txt")
    if not os.path.isfile(keys_path):
        raise FileNotFoundError(
            "AWS_Keys.txt not found. Create it in the project folder with "
            "access_key_id and secret_access_key from "
            "https://eodata-s3keysmanager.dataspace.copernicus.eu/"
        )
    access_key = ""
    secret_key = ""
    with open(keys_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if "access_key_id" in line:
                access_key = line.split(":")[-1].strip()
            if "secret_access_key" in line:
                secret_key = line.split(":")[-1].strip()
    if not access_key or not secret_key:
        raise ValueError(
            "AWS_Keys.txt is missing access_key_id or secret_access_key."
        )
    return access_key, secret_key


def _select_download_objects(objects):
    candidates = [obj for obj in objects if is_ch4_l2_nc(obj.key)]
    offl = [obj for obj in candidates if is_ch4_l2_nc(obj.key, ("OFFL",))]
    if offl:
        return offl, "OFFL"
    nrti = [obj for obj in candidates if is_ch4_l2_nc(obj.key, ("NRTI",))]
    if nrti:
        return nrti, "NRTI"
    return [], None


def download_TROPOMI_CH4_L2_data(start_date, end_date=None):
    import boto3

    access_key, secret_key = read_aws_keys()
    s3 = boto3.resource(
        "s3",
        endpoint_url="https://eodata.dataspace.copernicus.eu",
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        region_name="default",
    )
    bucket = s3.Bucket("eodata")

    prefixes = []
    for screening_date in iter_screening_dates(start_date, end_date):
        prefixes.append(
            f"Sentinel-5P/TROPOMI/L2__CH4___/{screening_date:%Y/%m/%d}/"
        )

    objects = []
    for prefix in prefixes:
        objects.extend(bucket.objects.filter(Prefix=prefix))

    selected, product_kind = _select_download_objects(objects)
    if not selected:
        raise FileNotFoundError(
            f"Could not find TROPOMI CH4 files between {start_date} and {end_date}."
        )

    target = os.path.abspath(os.path.join(os.getcwd(), "TROPOMI_data"))
    os.makedirs(target, exist_ok=True)
    existing = {
        name for name in os.listdir(target)
        if os.path.isfile(os.path.join(target, name))
    }

    count_download = 0
    for obj in selected:
        filename = os.path.basename(obj.key)
        if not filename or filename in existing:
            continue
        bucket.download_file(obj.key, local_download_path(target, filename))
        count_download += 1

    already_present = len(selected) - count_download
    if count_download == 0:
        download_message = (
            f"All {len(selected)} {product_kind} files are already downloaded."
        )
    else:
        download_message = (
            f"{count_download} {product_kind} files downloaded"
            f" ({already_present} already present)."
        )
    if product_kind == "NRTI":
        download_message += " OFFL was not available yet; used NRTI instead."
    return download_message


def inpolygon(xq, yq, xv, yv):
    vertices = np.vstack((xv, yv)).T
    path = Path(vertices)
    test_points = np.hstack([xq.reshape(xq.size, -1), yq.reshape(yq.size, -1)])
    _in = path.contains_points(test_points)
    _in_on = path.contains_points(test_points, radius=-1e-10)
    _on = _in ^ _in_on
    return _in_on, _on


def nanargmax(a):
    if np.all(np.isnan(a)):
        raise ValueError("Cannot find a maximum in an all-NaN array.")
    finite = np.where(np.isnan(a), -np.inf, a)
    return np.unravel_index(np.argmax(finite, axis=None), a.shape)


def find_nc_filenames(path_to_dir, sample_date, suffix=".nc"):
    date_str = sample_date.strftime("%Y%m%d")
    file_list = []
    for filename in os.listdir(path_to_dir):
        if filename.endswith(suffix) and tropomi_sensing_start_date(filename) == date_str:
            file_list.append(filename)
    return file_list


def find_nearest(array, value):
    idx = np.unravel_index((np.abs(array - value)).argmin(), array.shape)
    return idx


def quantification_mass_balance(enhancements, winds, pressures, lons, lats):
    if len(enhancements) == 0:
        return np.array([])
    lons, lats = np.array(lons), np.array(lats)
    enhancements = np.array(enhancements)
    winds = np.array(winds) * 3.6
    pressures = np.array(pressures)

    L = 0.05 * 111.32 * np.cos(lats * np.pi / 180) * 111.31 * 0.05
    Ls = np.array([sqrt(l) for l in L])
    return enhancements * 5.345 * (pressures / 1013) * Ls * winds * 2


def _grid_indices(lats, lons, grid_latitudes, grid_longitudes):
    lat_idx = np.rint((lats - grid_latitudes[0]) / GRID_RESOLUTION).astype(int)
    lon_idx = np.rint((lons - grid_longitudes[0]) / GRID_RESOLUTION).astype(int)
    lat_idx = np.clip(lat_idx, 0, grid_latitudes.size - 1)
    lon_idx = np.clip(lon_idx, 0, grid_longitudes.size - 1)
    return lat_idx, lon_idx


def _list_orbit_files(data_path, date_str):
    if not os.path.isdir(data_path):
        raise FileNotFoundError(
            f"TROPOMI_data folder not found at {data_path}. Download files first."
        )
    matches = []
    for filename in os.listdir(data_path):
        full_path = os.path.join(data_path, filename)
        if not os.path.isfile(full_path) or not is_ch4_l2_nc(filename):
            continue
        if tropomi_sensing_start_date(filename) == date_str:
            matches.append(full_path)
    offl = [path for path in matches if is_ch4_l2_nc(path, ("OFFL",))]
    return offl or matches


def _as_ndarray(variable):
    values = variable[:]
    if hasattr(values, "filled"):
        return np.asarray(values.filled(np.nan))
    return np.asarray(values)


def Load_CH4(minlat, maxlat, minlon, maxlon, sample_date, qa_pass=0.5):
    import netCDF4 as nc

    date_str = sample_date.strftime("%Y%m%d")
    daily_folder = "Daily_Global_TROPOMI_Concentration_Maps"
    os.makedirs(daily_folder, exist_ok=True)
    daily_file_path = os.path.join(daily_folder, f"{date_str}.nc")
    merge_file = not os.path.exists(daily_file_path)

    if merge_file:
        grid_longitudes = np.arange(-180, 180, GRID_RESOLUTION)
        grid_latitudes = np.arange(-90, 90, GRID_RESOLUTION)
        X, Y = np.meshgrid(grid_longitudes, grid_latitudes)
        fliped_Y = np.flipud(Y)
        global_ch4 = np.full(X.shape, np.nan, dtype=np.float32)
        global_winds = np.full(X.shape, np.nan, dtype=np.float32)
        global_qa = np.full(X.shape, np.nan, dtype=np.float32)
        global_pressures = np.full(X.shape, np.nan, dtype=np.float32)

        data_path = os.path.abspath(os.path.join(os.getcwd(), "TROPOMI_data"))
        file_list = _list_orbit_files(data_path, date_str)
        if not file_list:
            raise FileNotFoundError(
                f"No TROPOMI CH4 files for {date_str} in {data_path}."
            )

        all_ch4, all_qa, all_winds, all_pressures, all_lats, all_lons = (
            [], [], [], [], [], []
        )
        for file in file_list:
            data = nc.Dataset(file)
            xch4 = data.groups["PRODUCT"]["methane_mixing_ratio_bias_corrected"][:].data[0, :, :]
            valid_mask = np.isfinite(xch4) & (xch4 != CH4_FILL_VALUE) & (xch4 > 0)
            all_ch4.extend(xch4[valid_mask].tolist())
            latitudes = data.groups["PRODUCT"]["latitude"][:].data[0, :, :]
            longitudes = data.groups["PRODUCT"]["longitude"][:].data[0, :, :]
            all_lats.extend(latitudes[valid_mask].tolist())
            all_lons.extend(longitudes[valid_mask].tolist())
            qa_value = data.groups["PRODUCT"]["qa_value"][:].data[0, :, :]
            all_qa.extend(qa_value[valid_mask].tolist())
            u = _as_ndarray(
                data.groups["PRODUCT"].groups["SUPPORT_DATA"]["INPUT_DATA"]["eastward_wind"][0, :, :]
            )
            v = _as_ndarray(
                data.groups["PRODUCT"].groups["SUPPORT_DATA"]["INPUT_DATA"]["northward_wind"][0, :, :]
            )
            surface_pressure = _as_ndarray(
                data.groups["PRODUCT"].groups["SUPPORT_DATA"]["INPUT_DATA"]["surface_pressure"][0, :, :]
            )
            data.close()
            all_winds.extend(np.sqrt(u[valid_mask] ** 2 + v[valid_mask] ** 2).tolist())
            all_pressures.extend((surface_pressure[valid_mask] / 100).tolist())

        xdf = pd.DataFrame(
            {
                "xch4": all_ch4,
                "qa_value": all_qa,
                "wind_speed": all_winds,
                "pressure": all_pressures,
                "lon": np.round(all_lons, decimals=2),
                "lat": np.round(all_lats, decimals=2),
            }
        )
        xdf = xdf[xdf.qa_value >= qa_pass]
        if xdf.empty:
            raise FileNotFoundError(
                f"No quality-assured CH4 observations for {date_str}."
            )
        grouped_xdf = xdf.groupby(["lon", "lat"], as_index=False).mean()
        lat_idx, lon_idx = _grid_indices(
            grouped_xdf["lat"].to_numpy(),
            grouped_xdf["lon"].to_numpy(),
            grid_latitudes,
            grid_longitudes,
        )
        global_ch4[lat_idx, lon_idx] = grouped_xdf["xch4"].to_numpy()
        global_winds[lat_idx, lon_idx] = grouped_xdf["wind_speed"].to_numpy()
        global_qa[lat_idx, lon_idx] = grouped_xdf["qa_value"].to_numpy()
        global_pressures[lat_idx, lon_idx] = grouped_xdf["pressure"].to_numpy()

        flipped_ch4 = np.flipud(global_ch4)
        flipped_qa = np.flipud(global_qa)
        flipped_winds = np.flipud(global_winds)
        flipped_pressures = np.flipud(global_pressures)

        with nc.Dataset(daily_file_path, "w", format="NETCDF4") as ncfile:
            lat_size, lon_size = flipped_ch4.shape
            ncfile.createDimension("lat", lat_size)
            ncfile.createDimension("lon", lon_size)
            latitudes = ncfile.createVariable("latitude", "f4", ("lat", "lon"))
            longitudes = ncfile.createVariable("longitude", "f4", ("lat", "lon"))
            ch4 = ncfile.createVariable("CH4", "f4", ("lat", "lon"))
            wind = ncfile.createVariable("Wind", "f4", ("lat", "lon"))
            pressure = ncfile.createVariable("Pressure", "f4", ("lat", "lon"))
            qa_value = ncfile.createVariable("QA_value", "f4", ("lat", "lon"))
            latitudes[:, :] = fliped_Y
            longitudes[:, :] = X
            ch4[:, :] = flipped_ch4
            wind[:, :] = flipped_winds
            pressure[:, :] = flipped_pressures
            qa_value[:, :] = flipped_qa
    else:
        daily_data = nc.Dataset(daily_file_path)
        flipped_ch4 = daily_data.variables["CH4"][:].data
        flipped_winds = daily_data.variables["Wind"][:].data
        flipped_pressures = daily_data.variables["Pressure"][:].data
        X = daily_data.variables["longitude"][:].data
        fliped_Y = daily_data.variables["latitude"][:].data
        daily_data.close()

    north_row, west_col = find_indices(minlon, maxlat, X, fliped_Y)
    south_row, east_col = find_indices(maxlon, minlat, X, fliped_Y)
    row0, row1 = sorted((north_row, south_row))
    col0, col1 = sorted((west_col, east_col))

    fch4 = flipped_ch4[row0:row1 + 1, col0:col1 + 1]
    fwind = flipped_winds[row0:row1 + 1, col0:col1 + 1]
    fpressure = flipped_pressures[row0:row1 + 1, col0:col1 + 1]
    grid_lons = X[row0:row1 + 1, col0:col1 + 1]
    grid_lats = fliped_Y[row0:row1 + 1, col0:col1 + 1]
    fch4 = np.where((fch4 == 0) | ~np.isfinite(fch4), np.nan, fch4)
    return grid_lons, grid_lats, fch4, fwind, fpressure


def screening_plumes(ch4_obs, wind, pressure, grid_lons, grid_lats, threshold_delta, min_pixelcount):
    detected_plumes = []
    detected_plumes_lons = []
    detected_plumes_lats = []
    detected_plume_wind = []
    detected_plume_pressure = []

    if ch4_obs is None or np.size(ch4_obs) == 0:
        return (
            detected_plumes,
            detected_plume_wind,
            detected_plume_pressure,
            detected_plumes_lons,
            detected_plumes_lats,
        )

    x, y = np.shape(ch4_obs)
    i = 5
    while i < x - 6:
        j = 5
        while j < y - 6:
            patch = ch4_obs[i - 5:i + 6, j - 5:j + 6]
            if np.nansum(patch) > 0:
                mean_patch = np.nanmean(patch)
                median_patch = np.nanmedian(patch)
                std_patch = np.nanstd(patch, ddof=0)
                if std_patch > 0:
                    c = (mean_patch - median_patch) / std_patch
                    if c > 0.3:
                        bgd_patch = median_patch
                    else:
                        bgd_patch = (2.5 * median_patch) - (1.5 * mean_patch)
                    ano_patch = patch - bgd_patch - (3 * std_patch)
                    delta_patch = patch - bgd_patch
                    ind_patch = np.any(ano_patch > 0)
                    mean_delta = (
                        np.nanmean(delta_patch[delta_patch > 0])
                        if ind_patch
                        else np.nan
                    )
                    count_pixel = int(np.sum(delta_patch > 0))
                    if (
                        ind_patch
                        and mean_delta > threshold_delta
                        and count_pixel >= min_pixelcount
                    ):
                        detected_plumes.append(delta_patch)
                        detected_plume_wind.append(wind[i - 5:i + 6, j - 5:j + 6])
                        detected_plume_pressure.append(pressure[i - 5:i + 6, j - 5:j + 6])
                        detected_plumes_lons.append(grid_lons[i - 5:i + 6, j - 5:j + 6])
                        detected_plumes_lats.append(grid_lats[i - 5:i + 6, j - 5:j + 6])
            j += 1
        i += 1

    return (
        detected_plumes,
        detected_plume_wind,
        detected_plume_pressure,
        detected_plumes_lons,
        detected_plumes_lats,
    )


def unique_hotspots(detected_plumes, detected_plume_wind, detected_plume_pressure,
                    detected_plumes_lons, detected_plumes_lats):
    from shapely.geometry.polygon import Polygon

    polygons = []
    max_enhance = []
    max_lons = []
    max_lats = []
    max_winds = []
    max_pressures = []
    seen = set()

    for plume_lons, plume_lats, plume, plume_wind, plume_pressure in zip(
        detected_plumes_lons,
        detected_plumes_lats,
        detected_plumes,
        detected_plume_wind,
        detected_plume_pressure,
    ):
        ulon = np.max(plume_lons)
        llon = np.min(plume_lons)
        ulat = np.max(plume_lats)
        llat = np.min(plume_lats)
        pgon = Polygon((
            (llon, llat),
            (llon, ulat),
            (ulon, ulat),
            (ulon, llat),
            (llon, llat),
        ))
        i, j = nanargmax(plume)
        max_lon = float(plume_lons[i, j])
        max_lat = float(plume_lats[i, j])
        hotspot_key = (round(max_lon, 4), round(max_lat, 4))
        if hotspot_key in seen:
            continue
        seen.add(hotspot_key)
        max_enhance.append(plume[i, j])
        max_lons.append(max_lon)
        max_lats.append(max_lat)
        max_winds.append(plume_wind[i, j])
        max_pressures.append(plume_pressure[i, j])
        polygons.append(pgon)
    return polygons, max_enhance, max_lons, max_lats, max_winds, max_pressures


def generate_results(grid_lon, grid_lat, fch4, detected_plumes, detected_plume_wind,
                     detected_plume_pressure, detected_plumes_lons, detected_plumes_lats, date_str):
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cf

    polygons, max_enhance, max_lons, max_lats, max_winds, max_pressures = unique_hotspots(
        detected_plumes,
        detected_plume_wind,
        detected_plume_pressure,
        detected_plumes_lons,
        detected_plumes_lats,
    )
    Q = quantification_mass_balance(max_enhance, max_winds, max_pressures, max_lons, max_lats)

    target = os.path.abspath(os.path.join(os.getcwd(), "assets"))
    os.makedirs(target, exist_ok=True)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5), subplot_kw={"projection": ccrs.PlateCarree()})
    fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)
    tro = ax.pcolormesh(grid_lon, grid_lat, fch4, cmap="bwr", transform=ccrs.PlateCarree())
    ax.add_feature(cf.BORDERS)
    ax.coastlines()
    ax.set_xlim(np.nanmin(grid_lon) - 0.5, np.nanmax(grid_lon) + 0.5)
    ax.set_ylim(np.nanmin(grid_lat) - 0.5, np.nanmax(grid_lat) + 0.5)
    ax.stock_img()
    ax.set_title(f"Valid TROPOMI Methane Observations on {date_str}", fontsize=8.5)
    cbar = plt.colorbar(tro, pad=0.02, orientation="horizontal")
    cbar.set_label("Column average methane mixing ratio (ppb)", fontsize=8.5)
    for pgon in polygons:
        ax.add_geometries([pgon], crs=ccrs.PlateCarree(), facecolor="None", edgecolor="black")

    maxlon = np.round(np.nanmax(grid_lon))
    minlon = np.round(np.nanmin(grid_lon))
    maxlat = np.round(np.nanmax(grid_lat))
    minlat = np.round(np.nanmin(grid_lat))
    stem = f"TROPOMI_data_{date_str}_{maxlon}_{minlon}_{maxlat}_{minlat}"
    figure_name = os.path.join("assets", f"{stem}.jpg")
    fig.savefig(figure_name, dpi=300)
    plt.close(fig)

    df = pd.DataFrame({
        "Plume ID": range(len(max_enhance)),
        "Maximum Enhancement (ppb)": max_enhance,
        "longitude": max_lons,
        "latitude": max_lats,
        "Emission rate (kg/hr)": Q,
    })
    df.to_csv(os.path.join("assets", f"plumes_{date_str}_{maxlon}_{minlon}_{maxlat}_{minlat}.csv"),
              sep=",", index=False)
    return figure_name


def find_indices(target_lon, target_lat, X, Y):
    """
    target_lat: target latitude
    target_lon: target longitude
    X: 2D longitude array (-180,180)
    Y: 2D latitude array (90,-90)
    """
    lat_diff = np.abs(Y - target_lat)
    lon_diff = np.abs(X - target_lon)
    total_diff = lat_diff + lon_diff
    ind_x, ind_y = np.unravel_index(np.argmin(total_diff), total_diff.shape)
    return int(ind_x), int(ind_y)
