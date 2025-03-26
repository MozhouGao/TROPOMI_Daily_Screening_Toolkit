import os 
import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.path import Path
import boto3
from os import listdir
from os.path import isfile, join
import pandas as pd 
import cartopy.crs as ccrs
import cartopy.feature as cf
from shapely.geometry.polygon import Polygon
from math import sqrt
from datetime import datetime, timedelta

def read_aws_keys():
    path = os.getcwd()
    aws_keys_file = os.path.abspath(f"{path}/AWS_Keys.txt")
    AWS_keys = open(aws_keys_file, "r")
    access_key = ''
    secret_key = ''
    for line in AWS_keys:
        if "access_key_id" in line: 
            access_key = line.split(":")[-1].strip()
        if "secret_access_key" in line: 
            secret_key = line.split(":")[-1].strip()
    return access_key, secret_key


def download_TROPOMI_CH4_L2_data(start_date,end_date):
    access_key,secret_key = read_aws_keys()
    session = boto3.session.Session()
    s3 = boto3.resource(
        's3',
        endpoint_url='https://eodata.dataspace.copernicus.eu',
        aws_access_key_id=access_key,
        aws_secret_access_key=secret_key,
        region_name='default'
    ) 

    bucket = s3.Bucket("eodata")

    prefixes = []
    start_date = datetime.strptime(start_date, "%Y-%m-%d")
    end_date = datetime.strptime(end_date, "%Y-%m-%d")
    for i in range((end_date - start_date).days + 1):  # Include both start and end dates
        download_date = start_date + timedelta(days=i)
        download_date = download_date.strftime("%Y/%m/%d/")
        prefix = f"Sentinel-5P/TROPOMI/L2__CH4___/{download_date}" 
        prefixes.append(prefix)

    objects = [] 
    for prefix in prefixes: 
        objects += bucket.objects.filter(Prefix=prefix)
    
    if not list(objects):
        raise FileNotFoundError(f"Could not find any files for between {start_date} and {end_date}")
    else:
        path = os.getcwd()
        target = os.path.abspath(os.path.join(path, "TROPOMI_data"))
        if not os.path.exists(target):
            os.makedirs(target)

        files = [f for f in listdir(target) if isfile(join(target, f))]

    count_download = 0 
    for obj in objects:
        L2_file_name = obj.key.split('/')[-1]
        if L2_file_name not in files: 
            os.makedirs(os.path.dirname(obj.key), exist_ok=True)
            if not os.path.isdir(obj.key):
                bucket.download_file(obj.key, fr"{target}\{L2_file_name}")
                count_download += 1
    if count_download == 0: 
        download_message = "all files are already downloaded"
    else: 
        download_message = f"{count_download} files are downloaded"
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
    idx = np.argmax(a, axis=None)
    multi_idx = np.unravel_index(idx, a.shape)
    if np.isnan(a[multi_idx]):
        nan_count = np.sum(np.isnan(a))
        # In numpy < 1.8 use idx = np.argsort(a, axis=None)[-nan_count-1]
        idx = np.argpartition(a, -nan_count-1, axis=None)[-nan_count-1]
        multi_idx = np.unravel_index(idx, a.shape)
    return multi_idx

def find_nc_filenames( path_to_dir, date, suffix=".nc" ):
    
    filenames = os.listdir(path_to_dir)
    date_str = date.strftime('%Y%m%d')
    date_list = [date_str] 


    file_list = [] 
    for filename in filenames: 
        sample_time = filename.split("___")[1].split("_")[1].split('T')[0]
        if filename.endswith(suffix) and sample_time in date_list:
            file_list.append(filename)
    return file_list

def find_nearest(array, value):     
    idx = np.unravel_index((np.abs(array - value)).argmin(),array.shape)     
    return idx

def quantification_mass_balance(enhancements,winds,pressures,lons,lats): 
    lons,lats = np.array(lons), np.array(lats)
    enhancements = np.array(enhancements)
    winds = np.array(winds) * 3.6
    pressures = np.array(pressures)
    
    L = 0.05 * 111.32 * np.cos(lats*np.pi/180)*111.31*0.05
    Ls = [sqrt(l) for l in L] 
    Ls = np.array(Ls)
    
    Q = enhancements * 5.345 * (pressures/1013) * Ls * winds * 2 
    
    return Q 

def Load_CH4(minlat, maxlat, minlon, maxlon, date, qa_pass = 0.5):

    date_str = date.strftime("%Y%m%d")
    daily_folder = "Daily_Global_TROPOMI_Concentration_Maps"
    if not os.path.exists(daily_folder):  # Check if the folder exists
        os.makedirs(daily_folder)  # Create the folder
    
    daily_file_name = f"{date_str}.nc"
    daily_file_path = os.path.join(daily_folder, daily_file_name)
    if os.path.exists(daily_file_path):
        merge_file = False
    else:
        merge_file = True 

    if merge_file:
        # define a global 0.05 x 0.05 degree grid cells
        grid_longitudes = np.arange(-180,180,0.05)
        grid_latitudes = np.arange(-90,90,0.05)
        X,Y = np.meshgrid(grid_longitudes,grid_latitudes)
        fliped_Y =  np.flipud(Y)
        global_ch4 = np.empty(X.shape)
        global_winds = np.empty(X.shape)
        global_qa = np.empty(X.shape)
        global_pressures = np.empty(X.shape)
        # lons = np.arange(minlon, maxlon + 0.05, 0.05)
        # lats = np.arange(minlat, maxlat + 0.05, 0.05)
        # grid_lon,grid_lat = np.meshgrid(lons, lats)
        # a,b= np.shape(grid_lat)
        # p = Path([(minlat,minlon), (minlat, maxlon),
        #            (maxlat, maxlon), (maxlat, minlon)]) 
        # find file path 
        data_path = os.path.abspath(f"{os.getcwd()}/TROPOMI_data")
        
        # file_list = find_nc_filenames(local_path, date)
        file_list = []
        for f in listdir(data_path):
            if isfile(join(data_path, f)) and date_str in f.split('__')[-1].split('_')[1]:
                file_list.append(rf"{os.getcwd()}/TROPOMI_data/{f}")
        
        all_ch4 = [] 
        all_qa = [] 
        all_winds = [] 
        all_pressures = [] 
        all_lats = []
        all_lons = [] 
        for file in file_list: 
            # for each file, read data 
            data = nc.Dataset(file)
            xch4 = data.groups['PRODUCT']['methane_mixing_ratio_bias_corrected'][:].data[0,:,:]
            # create mask for CH4 concentration 
            mask = 9.96921e+36
            valid_mask = xch4 != mask
            masked_xch4 = xch4[valid_mask]
            all_ch4 += list(masked_xch4)
            # create mask for lat and long
            latitudes = data.groups['PRODUCT']['latitude'][:].data[0,:,:]
            longitudes = data.groups['PRODUCT']['longitude'][:].data[0,:,:]
            masked_latitudes = latitudes[valid_mask]
            all_lats += list(masked_latitudes)
            masked_longitudes = longitudes[valid_mask]
            all_lons += list(masked_longitudes)
            # create mask for QA value
            qa_value = data.groups['PRODUCT']['qa_value'][:].data[0,:,:]
            masked_qa = qa_value[valid_mask]
            all_qa += list(masked_qa)
            # create mask for wind speed 
            u = data.groups['PRODUCT'].groups["SUPPORT_DATA"]["INPUT_DATA"]["eastward_wind"][0,:,:]
            masked_u = u[valid_mask]
            v = data.groups['PRODUCT'].groups["SUPPORT_DATA"]["INPUT_DATA"]["northward_wind"][0,:,:]
            masked_v = v[valid_mask]
            # create mask for pressure 
            surface_pressure = data.groups['PRODUCT'].groups["SUPPORT_DATA"]["INPUT_DATA"]["surface_pressure"][0,:,:]
            masked_sp = surface_pressure[valid_mask] 
            data.close()
            
            # caluclate wind speed 
            ws = (masked_u**2 + masked_v**2)**0.5 
            all_winds += list(ws)
            # surface pressure  
            surface_pressure = masked_sp / 100 
            all_pressures += list(surface_pressure)
        
        # create dataframe for storing global valid CH4 concentration  
        valid_lat = np.round(all_lats,decimals=2)
        valid_lon = np.round(all_lons,decimals=2)
        valid_lat  = valid_lat.astype(str)
        valid_lon = valid_lon.astype(str)
        xdf = pd.DataFrame(data={'xch4':all_ch4,
                                 'qa_value':all_qa,
                                 'wind_speed': all_winds,
                                 'pressure': all_pressures,
                                 'lon':valid_lon,
                                 'lat':valid_lat  })
        xdf = xdf[xdf.qa_value >= qa_pass]
        grouped_xdf = xdf.groupby(['lon', 'lat'], as_index=False).mean()
    
        # create the array for storing global valid CH4 concentration 
        for _,row in xdf.iterrows(): 
            la = float(row.lat)
            lat_idx = np.argmin(np.abs(grid_latitudes - la))
            lo = float(row.lon)
            lon_idx = np.argmin(np.abs(grid_longitudes - lo))  
            global_ch4[lat_idx,lon_idx] = row.xch4
            global_winds[lat_idx,lon_idx] = row.wind_speed
            global_qa[lat_idx,lon_idx] = row.qa_value
            global_pressures[lat_idx,lon_idx] = row.pressure
    
        flipped_ch4 = np.flipud(global_ch4) 
        flipped_qa = np.flipud(global_qa)
        flipped_winds = np.flipud(global_winds)
        flipped_pressures = np.flipud(global_pressures)

        #save the daily TROPOMI concentrations 
        # Create a new NetCDF file
        with nc.Dataset(daily_file_path, "w", format="NETCDF4") as ncfile:
            # Get dimensions
            lat_size, lon_size = flipped_ch4.shape  # Assuming same shape for all arrays
        
            # Define dimensions
            ncfile.createDimension("lat", lat_size)
            ncfile.createDimension("lon", lon_size)
        
            # Create coordinate variables
            latitudes = ncfile.createVariable("latitude", "f4", ("lat", "lon"))
            longitudes = ncfile.createVariable("longitude", "f4", ("lat", "lon"))
        
            # Create data variables
            ch4 = ncfile.createVariable("CH4", "f4", ("lat", "lon"))
            wind = ncfile.createVariable("Wind", "f4", ("lat", "lon"))
            pressure = ncfile.createVariable("Pressure", "f4", ("lat", "lon"))
            qa_valie = ncfile.createVariable("QA_value","f4", ("lat","lon"))
        
            # Assign data to variables
            latitudes[:, :] = fliped_Y
            longitudes[:, :] = X
            ch4[:, :] = flipped_ch4
            wind[:, :] = flipped_winds
            pressure[:, :] =flipped_pressures
        
            # Add metadata 
            # ncfile.description = "Extracted dataset with CH4, wind, and pressure values."
            # latitudes.units = "degrees_north"
            # longitudes.units = "degrees_east"
            # ch4.units = "ppm"
            # wind.units = "m/s"
            # pressure.units = "hPa"

    else:
        daily_data = nc.Dataset(daily_file_path)
        flipped_ch4 = daily_data.variables["CH4"][:].data
        #flipped_qa = np.flipud(global_qa)
        flipped_winds = daily_data.variables["Wind"][:].data
        flipped_pressures = daily_data.variables["Pressure"][:].data
        X = daily_data.variables["longitude"][:].data
        fliped_Y = daily_data.variables["latitude"][:].data
        daily_data.close()

    # find indices of defined region
    upper_left_x,upper_left_y = find_indices(minlon,maxlat,X,fliped_Y)
    lower_right_x,lower_right_y = find_indices(maxlon,minlat,X,fliped_Y)

    fch4 = flipped_ch4[upper_left_x:lower_right_x,upper_left_y:lower_right_y]
    #midland_qa = flipped_qa[upper_left_x:lower_right_x,upper_left_y:lower_right_y]
    fwind = flipped_winds[upper_left_x:lower_right_x,upper_left_y:lower_right_y]
    fpressure = flipped_pressures[upper_left_x:lower_right_x,upper_left_y:lower_right_y]
    grid_lons = X[upper_left_x:lower_right_x,upper_left_y:lower_right_y]
    grid_lats = fliped_Y[upper_left_x:lower_right_x,upper_left_y:lower_right_y]

    fch4[fch4 == 0] = np.nan
    return grid_lons,grid_lats,fch4,fwind,fpressure


def screening_plumes(ch4_obs,wind,pressure,grid_lons,grid_lats,threshold_delta,min_pixelcount):
    num_detected_plume = 0 
    detected_plumes = [] 
    detected_plumes_lons = []
    detected_plumes_lats = []
    patch_checks = [] 
    mean_delta_checks = []
    detected_plume_wind = []
    detected_plume_pressure = [] 
    
    x,y = np.shape(ch4_obs)
    i = 5 
    while i<x-6: 
        j = 5
        while j<y-6:
            patch = ch4_obs[i-5:i+6,j-5:j+6]
            if np.nansum(patch) > 0: 
                mean_patch = np.nanmean(patch)
                median_patch = np.nanmedian(patch)
                std_patch = np.nanstd(patch,ddof=0)
                c = (mean_patch - median_patch)/std_patch
                if std_patch > 0:
                    if c > 0.3: 
                        bgd_patch = median_patch 
                    else: 
                        bgd_patch = (2.5*median_patch) - (1.5*mean_patch) 
                    # calculate the XCH4 anomaly in the patch
                    ano_patch = patch - bgd_patch - (3 * std_patch) 
                    # calculate the XCH4 enhancement in the patch
                    delta_patch = patch - bgd_patch
                    # define a suspect plume
                    ind_patch  = np.any(ano_patch>0)
                    if ind_patch:
                        mean_delta = np.nanmean(delta_patch[delta_patch>0])
                    count_pixel = len(delta_patch[delta_patch>0])
            
                    # Record the patches (XCH4, longitude, latitude) if all requirements fulfilled.
                    if ind_patch and mean_delta>threshold_delta and count_pixel > min_pixelcount: 
                        patch_checks.append(ano_patch)
                        mean_delta_checks.append(mean_delta)
                        detected_plumes.append(delta_patch)
                        detected_plume_wind.append(wind[i-5:i+6,j-5:j+6])
                        detected_plume_pressure.append(pressure[i-5:i+6,j-5:j+6])    
                        detected_plumes_lons.append(grid_lons[i-5:i+6,j-5:j+6])
                        detected_plumes_lats.append(grid_lats[i-5:i+6,j-5:j+6])
                        num_detected_plume += 1 
            j += 1 
        i += 1 
        
    return detected_plumes, detected_plume_wind, detected_plume_pressure, detected_plumes_lons, detected_plumes_lats 

def generate_results(grid_lon,grid_lat,fch4,detected_plumes, detected_plume_wind,
                     detected_plume_pressure , detected_plumes_lons,detected_plumes_lats, date_str):
    Polygon_list = [] 
    max_enhance = [] 
    max_lons = [] 
    max_lats = []
    max_winds = [] 
    max_pressures = [] 
    for plume_coor in zip(detected_plumes_lons,detected_plumes_lats,detected_plumes, detected_plume_wind, detected_plume_pressure):
        plume_lons = plume_coor[0]
        plume_lats = plume_coor[1]
        
        plume_wind = plume_coor[3]
        plume_pressure = plume_coor[4]
        
        ulon = np.max(plume_lons)
        llon = np.min(plume_lons)
        ulat = np.max(plume_lats)
        llat = np.min(plume_lats)
        pgon = Polygon(((llon,llat),
            (llon, ulat),
            (ulon, ulat),
            (ulon, llat),
            (llon,llat)))
    
        # find hotspot info 
        plume = plume_coor[2]
        (i,j) = nanargmax(plume)
        max_e = plume[i,j]
        max_lon = plume_lons[i,j]
        max_lat = plume_lats[i,j]
        max_wind =  plume_wind[i,j]
        max_pressure = plume_pressure[i,j]
        
        if max_lon not in max_lons and max_lat not in max_lats: 
            max_enhance.append(max_e)
            max_lons.append(max_lon)
            max_lats.append(max_lat)
            max_winds.append(max_wind) 
            max_pressures.append(max_pressure)
            Polygon_list.append(pgon)

    
    # Quantification 
    Q = quantification_mass_balance(max_enhance,max_winds,max_pressures,max_lons,max_lats)
    

    path = os.getcwd()
    target = os.path.abspath(os.path.join(path, "assets"))
    if not os.path.exists(target):
        os.makedirs(target)
                         
    fig, ax = plt.subplots(1, 1, figsize=(5,5),
                       subplot_kw={'projection': ccrs.PlateCarree()})

    fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)
    tro = ax.pcolormesh(grid_lon, grid_lat,fch4,cmap= "bwr",transform=ccrs.PlateCarree())
    ax.add_feature(cf.BORDERS)
    ax.coastlines()
    ax.set_xlim(np.min(grid_lon) -0.5, np.max(grid_lon)+0.5)
    ax.set_ylim(np.min(grid_lat) -0.5, np.max(grid_lat)+0.5)
    ax.stock_img()
    ax.set_title(f"Valid TROPOMI Methane Observations on {date_str}",fontsize = 8.5 )
    cbar = plt.colorbar(tro, pad=0.02, orientation= "horizontal")
    cbar.set_label('Column average methane mixing ratio (ppb)',fontsize=8.5)
    for pgon in Polygon_list:
        ax.add_geometries([pgon], crs=ccrs.PlateCarree(),facecolor="None",
                          edgecolor='black')
    
    maxlon = np.round(np.max(grid_lon))
    minlon = np.round(np.min(grid_lon)) 
    maxlat = np.round(np.max(grid_lat)) 
    minlat = np.round(np.min(grid_lat)) 
    
    figure_name =  fr"assets/TROPOMI_data_{date_str}_{maxlon}_{minlon}_{maxlat}_{minlat}.jpg"
    plt.savefig(figure_name,dpi=300)


    # Create results table 
    IDs = range(len(max_enhance))
    df = pd.DataFrame(data = {"Plume ID":IDs,
                                "Maximum Enhancement (ppb)":max_enhance,
                                "longitude":max_lons,
                                "latitude":max_lats,
                                "Emission rate (kg/hr)":Q})
    df.to_csv(fr"assets/plumes_{date_str}_{maxlon}_{minlon}_{maxlat}_{minlat}.csv",sep=',',index= False)

    return figure_name


def create_matplotlib_figure():
    fig, ax = plt.subplots()
    ax.plot([0, 1, 2, 3], [0, 1, 4, 9])

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)
    encoded_image = base64.b64encode(buf.read()).decode("utf-8")
    return f"data:image/png;base64,{encoded_image}"

def find_indices (target_lon,target_lat,X,Y):
    '''
    target_lat: target latitude
    target_lon: target longitude 
    X: 2D longitude array (-180,180)   
    Y: 2D latitude array (90,-90)
    '''
    lat_diff = np.abs(Y - target_lat)
    lon_diff = np.abs(X - target_lon)
    total_diff = lat_diff + lon_diff
    ind_x, ind_y = np.unravel_index(np.argmin(total_diff), total_diff.shape)
    return ind_x, ind_y 
