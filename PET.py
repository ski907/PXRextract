#-------------------------------------------------------------------------------
# Name          Precip Extraction Tool
# Description:  Tool to obtain 100-year 24-hr precip intensity and temporal
#               distribution globally from PXR data
#               PXR data from https://github.com/lrntct/pxr
# Author:       Chandler Engel
#               US Army Corps of Engineers
#               Cold Regions Research and Engineering Laboratory (CRREL)
#               Chandler.S.Engel@usace.army.mil
# Created:      31 January 2020
# Updated:      05 February 2020
#               
#-------------------------------------------------------------------------------

import xarray as xr
import os
import numpy as np


    
def convert_lon(longitude):
    """convert negative longitude into 360 longitude
    """
    return xr.where(longitude < 0, longitude + 360, longitude)
    # return np.where(longitude < 0, longitude + 360, longitude)
    
def lnt(T):
    return -np.log(1 - 1/T)

def y_gev_nonzero(T, shape):
    """
    Lu, L.-H., & Stedinger, J. R. (1992).
    Variance of two- and three-parameter GEV/PWM quantile estimators:
    formulae, confidence intervals, and a comparison.
    Journal of Hydrology, 138(1–2), 247–267.
    http://doi.org/10.1016/0022-1694(92)90167-T
    """
    return (1 - lnt(T)**shape) / shape    

def loadPXRdata(data_dir):
    """Loads PXR netcdf files, expects them both to be in the same dir
    
    """
    
    PXR2 = f'pxr2-2.1.0.nc'
    PXR4 = f'pxr4-2.1.0.nc'
    
    ds_pxr2 = xr.open_dataset(os.path.join(data_dir, PXR2))
    ds_pxr4 = xr.open_dataset(os.path.join(data_dir, PXR4))
    
    return ds_pxr2, ds_pxr4

def extract_intensity_from_PXR2(ds,site_coord,duration,T):
        
    ds_cont_extract = ds.sel(latitude=site_coord[0], longitude=convert_lon(site_coord[1]), method='nearest')
    duration_index = ds_cont_extract['duration']==duration
    loc = ds_cont_extract.location[duration_index]
    scale = ds_cont_extract.scale[duration_index]
    
    y = y_gev_nonzero(T,-0.114)
    
    intensity = loc + scale * y
    
    return intensity

def extract_PXR2_duration_intensities(ds,site_coord,T):
    """Extract all duration intensities for a given location
    and recurrance interval"""
    
    ds_cont_extract = ds.sel(latitude=site_coord[0], longitude=convert_lon(site_coord[1]), method='nearest')
    
    y = y_gev_nonzero(T,-0.114)
    
    intensities = ds_cont_extract.location.values + ds_cont_extract.scale.values * y
    
    return intensities

def create_24_hr_alternating_block_hyetograph(intensities):
    
    durations = [1,2,3,4,6,8,10,12,18,24,48,72,96,120,144,192,240,288,360]
    
    precip = np.zeros(24)
    
    dt = np.diff(durations)
    
    cum_depth = intensities*durations
    inc_depth_irr = np.diff(cum_depth)/dt
    
    inc_depth = np.array(cum_depth[0])
    
    for i,d in enumerate(dt):
        inc_depth = np.append(inc_depth,np.ones(d)*inc_depth_irr[i])
  
    alt_block_index = [12,13,11,14,10,15,9,16,8,17,7,18,6,19,5,20,4,21,3,22,2,23,1,0]    
        
    for i,ind in enumerate(alt_block_index):
        precip[ind] = inc_depth[i]
    #test below
    precip = np.append(precip,precip[-1])       
    
    return precip
    

def extract_intensity_from_PXR4(ds,site_coord,duration,T):
    
    a, alpha, b, beta, y = getPXR4_params(ds,site_coord,T)
    
    intensity = a*duration**alpha + b*duration**beta*y
    
    return intensity

def getPXR4_params(ds,site_coord,T):
    
    ds_cont_extract = ds.sel(latitude=site_coord[0], longitude=convert_lon(site_coord[1]), method='nearest')
    
    a = ds_cont_extract.a
    alpha = ds_cont_extract.alpha
    b = ds_cont_extract.b
    beta = ds_cont_extract.beta
    y = y_gev_nonzero(T,-0.114)
    
    return a, alpha, b, beta, y
    

def temporal_dist_PXR4(ds,site_coord,total_duration,T,Tp):
    """Compute hyetograph from PXR4 intensity function, user specified
    location, total storm duration (hours), recurrance interval (years), and
    time to peak (hours)    instantaneous intensity model.
    """
    
    r = Tp/total_duration
    a, alpha, b, beta, y = getPXR4_params(ds,site_coord,T)
    time_steps = np.linspace(0,total_duration,total_duration*10+1) #6 minute time steps
    
    V=np.zeros(len(time_steps))
    
    for i, t in enumerate(time_steps):
        if t < Tp:
            V[i] = -1*(Tp-t)*(a*((Tp-t)/r)**alpha+b*(((Tp-t)/r)**beta)*y)
        elif t > Tp:
            V[i] = (t-Tp)*(a*((t-Tp)/(1-r))**alpha+b*(((t-Tp)/(1-r))**beta)*y)
        elif t == Tp:
            V[i] = 0
        if t == 0:
            V0= -V[i].copy()
        V[i] = V[i]+V0
       
    return V

def create_unit_rainfall_dist(V):
    
    V_unit=V/V[-1]
    
    return V_unit

def gen_cum_precip(unit_rainfall_dist,total_precip):
    
    cum_precip = []
    
    for x in unit_rainfall_dist:
        cum_precip.append(x*total_precip)
        
    return cum_precip


def export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,dt,output_file):
    
    #from datetime import datetime
    from pydsstools.heclib.dss import HecDss
    from pydsstools.core import TimeSeriesContainer,UNDEFINED

    precip = np.diff(cum_rainfall_dist) #convert cumulative rainfall to incremental
    precip = np.insert(precip,0,0)
    
    if dt < 60:
        time_string = str(dt)+"MIN"
    elif dt == 60:
        time_string = "1HOUR"

    dss_file = output_file
    pathname = r"/"+installation+r"/"+installation+r"/PRECIP_INC//"+time_string+r"/"+Part_F+r"/"
    tsc = TimeSeriesContainer()
    tsc.pathname = pathname
    tsc.startDateTime = "01JAN2000 00:00:00"
    tsc.numberValues = len(precip)
    tsc.units = "mm"
    tsc.type = "PER-CUM"
    tsc.interval = 1
    tsc.values = precip




    fid = HecDss.Open(dss_file,version=6)
    fid.deletePathname(tsc.pathname)
    fid.put_ts(tsc)
    ts = fid.read_ts(pathname)
    fid.close()


    
def main():
    
    data_dir = r'C:\Users\RDCRLCSE\Documents\Python Scripts\PXRextract\PXR'
    
    ds_pxr2, ds_pxr4 = loadPXRdata(data_dir)

    duration = 24 #storm duration, hours for selecting overall average storm intensity
    T = 100 #return period, years

    Tp=12   #time to peak, hours
    total_duration = duration #storm duration, sets length of storm in output series. Usually the same as duration
    
    site_coord = [44.564624, -72.761110]    #get site coords, lat long
 
    installation = 'Freefall'
    output_file = 'Global_Precip_'+installation+"_metric.dss"
    
    ####PXR2 alternating block#######
    intensities = extract_PXR2_duration_intensities(ds_pxr2,site_coord,T)
    precip = create_24_hr_alternating_block_hyetograph(intensities)
    
    cum_rainfall_dist = np.cumsum(precip)
    Part_F = '100-yr PXR2'
        
    export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,60,output_file) #60 minute blocks
    
    print(intensities)
    
    ##############PXR4 chicago method########
    V=temporal_dist_PXR4(ds_pxr4,site_coord,total_duration,T,Tp)
            
    cum_rainfall_dist = V
    Part_F = '100-yr PXR4'
        
    export_rainfall_toDSS(installation,Part_F,cum_rainfall_dist,6,output_file) #six minute blocks
    
    #############PXR2 depth PXR4 distribution###########
    intensity_PXR2 = extract_intensity_from_PXR2(ds_pxr2,site_coord,duration,T)
    V_unit=create_unit_rainfall_dist(V)
    PXR2_PXR4_dist = gen_cum_precip(V_unit,intensity_PXR2.values[0]*24)
    Part_F = '100-yr PXR2 Depth PXR4 Dist'
    
    export_rainfall_toDSS(installation,Part_F,PXR2_PXR4_dist,6,output_file)
    
if __name__ == "__main__":
    main()