import os 
import netCDF4 as nc
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.path import Path
from sentinelsat import SentinelAPI
import pandas as pd 
import cartopy.crs as ccrs
import cartopy.feature as cf
from datetime import timedelta
from shapely.geometry.polygon import Polygon

def download_TROPOMI_CH4_L2_data(start_date,end_date,west_lon, east_lon, south_lat, north_lat): 
    
    #create a folder to save data 
    path = os.getcwd()
    local_path = os.path.abspath(f"{path}/TROPOMI_data")
    if not os.path.exists(local_path):
        #print("create folder to save downloaded data")
        os.mkdir(local_path)

    # Query S5P Data Hub using guest login credentials
    api = SentinelAPI('s5pguest', 's5pguest', 'https://s5phub.copernicus.eu/dhus')

    # define date 
    query_date=(start_date,end_date)

    # define polygon 
    west_lon = str(west_lon)
    east_lon = str(east_lon)
    south_lat = str(south_lat)
    north_lat = str(north_lat)
    footprint = 'POLYGON((' + west_lon + ' ' + south_lat + ',' + east_lon + ' ' + south_lat + ',' + east_lon + ' ' + north_lat + ',' + west_lon + ' ' + north_lat + ',' + west_lon + ' ' + south_lat + '))'

    # query data product
    try:
        products = api.query(footprint,date=query_date,platformname="Sentinel-5 Precursor",producttype="L2__CH4___")
        if bool(len(products)>0): 
            # Before downloading, check if the file already exist in the local drive.
            downloaded_keys = [] 
            for key in products.keys():  
                file_title = products[key]['title']
                file_path = os.path.abspath(f"{local_path}/{file_title}.nc")
                #print(file_path)
                if os.path.exists(file_path):
                    downloaded_keys.append(key)
            for key in downloaded_keys:
                products.pop(key)
        else: 
            return f"Data files not available between {start_date} and {end_date}"
        if bool(len(products)>0):
            products_df = api.to_dataframe(products)
            # log data 
            file_name_list = []
            download_prefix = ""
            for ele in zip(products_df['title'],products_df['size']):
                download_prefix = download_prefix + f" Data files downloading: {ele[0]}; size: {ele[1]}"
                file_name_list.append(ele[0])
            # download all results from the search
            try:
                api.download_all(products,local_path)
                return download_prefix
            except: 
                return "Too many requests, trying to download again"
        else: 
            return "All data files are already exist in the local drive"
    except:
        return "Error connecting to SciHub server"

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

def find_nc_filenames( path_to_dir, start_date, end_date, suffix=".nc" ):
    
    filenames = os.listdir(path_to_dir)
    
    date_list = [] 
    while start_date.day <= end_date.day: 
        start_date_str = start_date.strftime('%Y%m%d')
        date_list.append(start_date_str )
        start_date += timedelta(days = 1)

    file_list = [] 
    for filename in filenames: 
        sample_time = filename.split("___")[1].split("_")[1].split('T')[0]
        if filename.endswith(suffix) and sample_time in date_list:
            file_list.append(filename)
    return file_list

def find_nearest(array, value):     
    idx = np.unravel_index((np.abs(array - value)).argmin(),array.shape)     
    return idx

def Load_CH4(minlat, maxlat, minlon, maxlon, start_date, end_date, qa_pass = 0.5): 
    lons = np.arange(minlon, maxlon + 0.05, 0.05)
    lats = np.arange(minlat, maxlat + 0.05, 0.05)
    grid_lon,grid_lat = np.meshgrid(lons, lats)
    a,b= np.shape(grid_lat)
    p = Path([(minlat,minlon), (minlat, maxlon),
               (maxlat, maxlon), (maxlat, minlon)]) 
    # find file path 
    local_path = os.path.abspath(f"{os.getcwd()}/TROPOMI_data")
    
    file_list = find_nc_filenames(local_path, start_date, end_date)
    
    all_ch4 = [] 
    for file_name in file_list: 
        # for each file, read data 
        file_path = os.path.abspath(f"{local_path}/{file_name}")
        data = nc.Dataset(file_path)
        xch4 = data.groups['PRODUCT']['methane_mixing_ratio_bias_corrected'][:].data[0,:,:]
        xch4[xch4 == 9.96921e+36] = np.nan
        latitudes = data.groups['PRODUCT']['latitude'][:].data[0,:,:]
        longitudes = data.groups['PRODUCT']['longitude'][:].data[0,:,:]
        qa_value = data.groups['PRODUCT']['qa_value'][:].data[0,:,:]
        data.close()
        # find dims of xch4 
        x,y = np.shape(xch4)
        # find data in the defined region 
        frame_lons = [] 
        frame_lats = []
        frame_ch4 = [] 
        frame_qa = [] 
        i = 0 
        while i < x: 
            j= 0
            while j < y:
                lo = longitudes[i,j]
                la = latitudes[i,j] 
                check1 = p.contains_points([(la, lo)])
                if check1[0]: 
                    frame_lons.append(lo)
                    frame_lats.append(la)
                    frame_ch4.append(xch4[i,j])
                    frame_qa.append(qa_value[i,j])
                j += 1 
            i += 1
            
        zoom_xch4 = np.zeros(shape = (a,b))  
        for val in zip(frame_lats,frame_lons,frame_ch4,frame_qa):
            if val[3] >=qa_pass:
                lo = val[1]
                la = val[0]
                yi = find_nearest(grid_lon,lo)[1]
                xi = find_nearest(grid_lat,la)[0]

                if zoom_xch4[xi,yi] > 0: 
                    zoom_xch4[xi,yi] =  (zoom_xch4[xi,yi] + val[2])/2 
                else:
                    zoom_xch4[xi,yi] = val[2]
        zoom_xch4[zoom_xch4 == 0] = np.nan
        all_ch4.append(zoom_xch4)
        
    all_ch4 = np.array(all_ch4)
    fch4 = np.nanmean(all_ch4,axis=0)
    
    return grid_lon,grid_lat,fch4


def screening_plumes(ch4_obs,grid_lons,grid_lats,threshold_delta,min_pixelcount):
    num_detected_plume = 0 
    detected_plumes = [] 
    detected_plumes_lons = []
    detected_plumes_lats = []
    patch_checks = [] 
    mean_delta_checks = []
    x,y = np.shape(ch4_obs)
    i = 5 
    while i< x-6: 
        j = 5
        while j< y-6:
            patch = ch4_obs[i-5:i+6,j-5:j+6]
            mean_patch = np.nanmean(patch)
            median_patch = np.nanmedian(patch)
            std_patch = np.nanstd(patch,ddof=1)
            c = (mean_patch - median_patch)/std_patch
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
            mean_delta = np.nanmean(delta_patch[delta_patch>0])
            count_pixel = len(delta_patch[delta_patch>0])

            # Record the patches (XCH4, longitude, latitude) if all requirements fulfilled.
            if ind_patch and mean_delta>threshold_delta and count_pixel > min_pixelcount: 
                patch_checks.append(ano_patch)
                mean_delta_checks.append(mean_delta)
                detected_plumes.append(delta_patch)
                detected_plumes_lons.append(grid_lons[i-5:i+6,j-5:j+6])
                detected_plumes_lats.append(grid_lats[i-5:i+6,j-5:j+6])
                num_detected_plume += 1 
            j += 1 
        i += 1 
        
    return detected_plumes, detected_plumes_lons, detected_plumes_lats 

def create_figures(grid_lon,grid_lat,fch4,detected_plumes, detected_plumes_lons,detected_plumes_lats, date1, date2):
    Polygon_list = [] 
    max_enhance = [] 
    max_lons = [] 
    max_lats = []
    for plume_coor in zip(detected_plumes_lons,detected_plumes_lats,detected_plumes):
        plume_lons = plume_coor[0]
        plume_lats = plume_coor[1]

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
        max_lon = plume_lons [i,j]
        max_lat = plume_lats[i,j]
        
        if max_lon not in max_lons and max_lat not in max_lats: 
            max_enhance.append(max_e)
            max_lons.append(max_lon)
            max_lats.append(max_lat)
            Polygon_list.append(pgon)


        
    fig, ax = plt.subplots(1, 1, figsize=(5,5),
                       subplot_kw={'projection': ccrs.PlateCarree()})

    fig.subplots_adjust(hspace=0, wspace=0, top=0.925, left=0.1)
    tro = ax.pcolormesh(grid_lon, grid_lat,fch4,cmap= "bwr",transform=ccrs.PlateCarree())
    ax.add_feature(cf.BORDERS)
    ax.coastlines()
    ax.set_xlim(np.min(grid_lon) -0.5, np.max(grid_lon)+0.5)
    ax.set_ylim(np.min(grid_lat) -0.5, np.max(grid_lat)+0.5)
    ax.stock_img()
    ax.set_title(f"Valid TROPOMI Methane Observations from {date1} to {date2}",fontsize = 8.5 )
    cbar = plt.colorbar(tro, pad=0.02, orientation= "horizontal")
    cbar.set_label('Column average methane mixing ratio (ppb)',fontsize=8.5)
    for pgon in Polygon_list:
        ax.add_geometries([pgon], crs=ccrs.PlateCarree(),facecolor="None",
                          edgecolor='black')
    
    maxlon = np.round(np.max(grid_lon))
    minlon = np.round(np.min(grid_lon)) 
    maxlat = np.round(np.max(grid_lat)) 
    minlat = np.round(np.min(grid_lat)) 
    
    figure_name =  fr"assets/TROPOMI_data_{date1}_{date2}_{maxlon}_{minlon}_{maxlat}_{minlat}.jpg"
    plt.savefig(figure_name,dpi=300)
    figure_path =  figure_name


    # Create results table 
    IDs = range(len(max_enhance))
    df = pd.DataFrame(data = {"Plume ID":IDs,
                                "Maximum Enhancement":max_enhance,
                                "longitude":max_lons,
                                "latitude":max_lats})
    df.to_csv(fr"assets/plumes_{date1}_{date2}_{maxlon}_{minlon}_{maxlat}_{minlat}.csv",sep=',')

    return figure_path