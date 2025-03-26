# TROPOMI Daily Screening Toolkit - Version 1.1 
### Release Notes
#### Updates:

- 1.TROPOMI Data Download: Migrated from Sentinelsat to Boto3 for data retrieval: https://documentation.dataspace.copernicus.eu/APIs/S3.html#accessing-eodata-via-aws-cli
-   To enable API access, follow the EODATA credential generation guide: https://eodata-s3keysmanager.dataspace.copernicus.eu/
-   If you do not have an account, register here: https://documentation.dataspace.copernicus.eu/Registration.html
-   Once credentials are generated, store your access and secret key in a file named AWS_Keys.txt under the workspace. Set up access credentials using the provided example.
![image](https://github.com/user-attachments/assets/fa541963-3e3e-4271-bad1-c390e040d830)

- 2.Added map-based selection.

# TROPOMI Daily Screening Toolkit - Version 1.0 
## Summary
This toolkit is developed to automatically screen the suspect methane (CH4) plumes over a user defined region based on the public accessible satellite observations (i.e., methane dry air mixing ratio from TROPOMI). Users are required to input several parameters to complete the screening process, including: the screening date(s), region boundaries, screening criteria (i.e., threshold enhancement, number of valid plume pixels. The output of each screening run is a XCH4 concentration map by highlighting the regions with high probability of detecting suspect methane plume(s). A list of the potential source locations is also available from the screening result. The noteworthy point is that this toolkit is designed neither for pinpointing the emission source at the facility or component level, nor for screening of small methane leaks (generally <25 tons/hour). For source attributions of the detected suspect plumes, the follow-up targeted fine-scale observations over the regions with high probability of detecting suspect methane plume(s) are required.

## Installation 
The following were steps taken to run TROPOMI Daily Screening Tookit
Install packages from command Shell: 
- if you are using pip environment 
  -   pip freeze > requirements.txt
  -   pip install cartopy or pip install git+https://github.com/SciTools/cartopy.git
- if you are using conda envionment 
  -   conda list -e > requirements.txt   
  -   conda install -c conda-forge cartopy
  
  
## Internal default settings
The datasets are the Level-2 user products from the TROPOMI instrument onboard the Sentinel-5P satellite. Original data files are obtained from the Sentinels Scientific Data Hub (currently known as The Copernicus Open Access Hub). The atmospheric methane concentrations are indicated by the column averaged dry air mole mixing ratio of methane (XCH4) with the unit of ppb. For each data file, the methane_mixing_ratio_bias_corrected was retrieved and remapped on a 0.05-deg by 0.05-deg map. To ensure that the remapped observations are valid, only quality assured data points (i.e., the qa_value greater than 0.5, as recommended by the data provider) were used in the plume screening and mapping. The basic workflow in this toolkit is 1) automatically downloading TROPOMI data files, 2) loading methane observations, 3) segmenting methane anomaly, 4) recording suspect methane emitting region, and 5) finalizing results.
The plume screening algorithm originally proposed by Lauvaux et al. (2022) and was simplified to screen the suspect methane plumes and highlight the regions with high probability of detecting methane ultra emitters. In the plume screening algorithm, the suspect plume screening is patch-based. The geographic span of each patch is 0.5-deg by 0.5-deg (roughly 55km by 55km at the low latitudes). For each patch, the mean, median, and standard deviation of the observed XCH4 is calculated to assess the distributions of methane observations and further find the anomaly to indicate the presence of suspect methane plumes. 

### Background XCH4 
Determination of the background XCH4 is critical for the patch-based plume screening, considering this could significantly influence the calculations of methane enhancement. In contrast to taking the value of the pixel in the vicinity of a detected plume in the upwind direction, the background XCH4 during the patch-based plume screening in this toolkit is computed according to a matrix $$((XCH4,mean- XCH4,median)/ XCH4,std)$$ indicating the skewness of the pixel distribution. For the patches with XCH4 values are strongly skewed $$((XCH4,mean- XCH4,median)/ XCH4,std > 0.3)$$, background XCH4 is the median of the patch. Otherwise, the background XCH4 is computed as. This method is commonly used for robust background estimation in noisy astronomical images analysis. The  s a tunable parameter, which is typically set at 2.5 to successfully capture some of the well-known methane emission events (Lauvaux et al., 2022).

## User defined parameters 
A well-fitted plume shape is another factor that may influence the identification of the methane emission plume. In respect of plume shape, two parameters are used to reflect the fitness of the detected anomalies to a suspect plume, which are the mean XCH4 enhancement over the pixels with positive enhancements (ΔXCH4,thr) and the number of valid pixels with positive enhancements (n). Under the conditions that these two parameters fulfilled the requirements, a suspect plume will be tagged. The users could define these two parameters based on their demand to obtain either conservative or radical results (e.g., the larger values for ΔXCH4,thr and n, the more conservative results would be). 

## Quick guide
### Select the date
Please specify an individual date or a period for daily screening Multiple days screening may cause longer data processing time. 
Note: If a period is selected, only the screening result of the last day will be displayed on the webpage. For full list of the daily screening results, check the local path: ~/TROPOMI_Daily_Screening_Toolkit-main/assets.
### Define regions
Currently the toolkit only supports screening over a rectangular region. Click <Define polygon> to confirm.
Longitude range: -180 ~ 180, latitude range: -90 ~ 90.
### Download Level-2 TROPOMI methane observations
Click <Download> to download the data files to the local path: ~/TROPOMI_Daily_Screening_Toolkit-main/TROPOMI_data.
### Start screening
Enter the Threshold delta (ΔXCH4,thr) and Minimum pixel count (n). Then click <Start screening> to kick off the daily plume screening. The screening time may vary with region size and number of days. Please do NOT hit on <Start screening> multiple times. Thanks for your patience.

## About the screening results:
If any suspect methane plumes are detected, a methane concentration map with the highlighted plume patches will be provided in the tool. If applicable, the potential locations of the detected methene emissions sources would be provided in the format of .csv file. The users could retrieve the results, including the maps and list of suspect sources, from the local path: ~/TROPOMI_Daily_Screening_Toolkit-main/assets
- The daily averaged methane mixing ratio is calculated based on multiple observations from the TROPOspheric Monitoring Instrument (TROPOMI, i.e., satellite instrument on board the Copernicus Sentinel-5 Precursor satellite. Only valid observations (i.e., observations with qa_value greater than 0.5) were used. However, data filtering by qa_value could not guarantee that all the bad data points are eliminated.
- The location of the suspect plume is determined as patch based. For each 11×11 patch, we tag the grid cell (0.05° × 0.05°) with the maximum observed XCH4 as the potential location of the detected plume. As limited by the spatial resolution of TROPOMI observations (i.e., 7 km × 5.5 km) and the completeness of the plume puzzles, the "locations" indicated by this toolkit refer to the suspect regions with higher probability of detecting methane emissions.

## Feedback
For any questions and concerns may you have about this toolkit, please feel free to get in touch with zhenyu.xing2@ucalgary.ca and mozhou.gao@ucalgary.ca.

## References
Lauvaux, T., Giron, C., Mazzolini, M., d’Aspremont, A., Duren, R., Cusworth, D., Shindell, D. and Ciais, P., 2022. Global assessment of oil and gas methane ultra-emitters. Science, 375(6580), pp.557-561.


