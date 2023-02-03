# TROPOMI daily Screening Toolkit - Version 1.0 
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
  
  


