#-------------------------------------------------------------------------------
# Name          Multi Site Extract Precip Tool
# Description:  Tool to get extract PXR2/PXR4 precip data
#               Takes a json file of station name and location
#               Outputs dss files with precip timeseries
# Author:       Chandler Engel
#               US Army Corps of Engineers
#               Cold Regions Research and Engineering Laboratory (CRREL)
#               Chandler.S.Engel@usace.army.mil
# Created:      31 January 2020
# Updated:      05 February 2020
#               
#-------------------------------------------------------------------------------
import PET
import json
import numpy as np

data_dir = r'C:\Users\RDCRLCSE\Documents\Python Scripts\PXRextract\PXR'
    
ds_pxr2, ds_pxr4 = PET.loadPXRdata(data_dir)

duration = 24 #storm duration, hours for selecting overall average storm intensity
T = 100 #return period, years

Tp=12   #time to peak, hours
total_duration = 24 #storm duration, sets length of storm in output series. Usually the same as duration

with open('Import_Locations.txt', 'r') as filehandle: #load site locations
    sites = json.load(filehandle)

for site in sites:  #loop through each site
    site_coord = site[1]    #get site coords, lat long
 
    installation = site[0]
    output_file = 'Global_Precip_'+installation+"_metric.dss"
    
    ####PXR2 alternating block#######
    intensities = PET.extract_PXR2_duration_intensities(ds_pxr2,site_coord,T)
    precip = PET.create_24_hr_alternating_block_hyetograph(intensities)
    
    cum_rainfall_dist = np.cumsum(precip)
    Part_F = '100-yr PXR2'
        
    PET.export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,60,output_file)
    
    ##############PXR4 chicago method########
    V=PET.temporal_dist_PXR4(ds_pxr4,site_coord,total_duration,T,Tp)
            
    cum_rainfall_dist = V
    Part_F = '100-yr PXR4'
        
    PET.export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,6,output_file)
    
    #############PXR2 depth PXR4 distribution###########
    intensity_PXR2 = PET.extract_intensity_from_PXR2(ds_pxr2,site_coord,duration,T)
    V_unit=PET.create_unit_rainfall_dist(V)
    PXR2_PXR4_dist = PET.gen_cum_precip(V_unit,intensity_PXR2.values[0]*24)
    Part_F = '100-yr PXR2 Depth PXR4 Dist'
    
    PET.export_rainfall_toDSS(installation,Part_F,PXR2_PXR4_dist,6,output_file)
    

    