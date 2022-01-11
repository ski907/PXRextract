#-------------------------------------------------------------------------------
# Name          Get Centroids Tool
# Description:  Tool to get centroid coordinates from input feature class AOI
#               Outputs in a json file in the format that multi_site_extract.py
#               expects for computing precip and exporting dss timeseries
# Author:       Chandler Engel
#               US Army Corps of Engineers
#               Cold Regions Research and Engineering Laboratory (CRREL)
#               Chandler.S.Engel@usace.army.mil
# Created:      31 January 2020
# Updated:      05 February 2020
#               
#-------------------------------------------------------------------------------

#use something like execfile(r'c:\my\script.py') to run from ArcMap Python window

import arcpy
import json
arcpy.CheckOutExtension("Spatial")

def setworkspace(directory):
    arcpy.env.workspace = directory

def calculate_centroid(shape,outdir):
    """Returns a point feature from a provided polygon
       
    """
    
    out_file_string = outdir+r"/centroid_of_"+shape    
    arcpy.FeatureToPoint_management(in_features=shape, out_feature_class=out_file_string, point_location="CENTROID")
    arcpy.AddXY_management(in_features=out_file_string)
    
    cursor = arcpy.da.SearchCursor(out_file_string, ['POINT_X', 'POINT_Y'])
    for row in cursor:
        lon = row[0]
        lat = row[1]
        
    return [lat, lon]

def main():        
    setworkspace(r"C:\Users\u4rrecse\Documents\RAS OCONUS\OCONUS_DEM_Extents_For_Requested_Sites_20200818\split")
    outdir = r"C:\Users\u4rrecse\Documents\RAS OCONUS\OCONUS_DEM_Extents_For_Requested_Sites_20200818\split\TempWorking"
    
    #"C:\Users\u4rrecse\Documents\Global Precip\Precip Extraction Tool\get_centroids.py"
    
    featureclasses = arcpy.ListFeatureClasses()
    
    locations = []
    for fc in featureclasses:
        #loc_name = fc[0:-23] #strips off last characters from batch sent JAN2020, this is a little too hardcoded
        loc_name = fc[0:-4]
        locations.append([loc_name, calculate_centroid(fc,outdir)])
        
    with open(r'C:\Users\u4rrecse\Documents\RAS OCONUS\OCONUS_DEM_Extents_For_Requested_Sites_20200818\split\TempWorking\OCONUS_Locations.txt', 'w') as filehandle:
        json.dump(locations,filehandle)
    
if __name__ == "__main__":
    main()