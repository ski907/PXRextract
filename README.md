﻿# PXR Precipitation Extraction Tool
## Summary
The PXR Precip Extraction Tool (PET) is used to extract precipitation data from the Parameterized eXtreme Rain (PXR) dataset (Courty et al. 2019) and can produce a DSS database of precipitation that can be used to define the precipitation boundary condition of a HEC-RAS or HEC-HMS hydrology model.

## Background
The PET was written to automate the generation of 24-hour 100-year storm precipitation boundary conditions for HEC-RAS models of **OCONUS** watersheds as part of a USACE climate project in early 2020. The tool was created because [NOAA Atlas 14](https://hdsc.nws.noaa.gov/hdsc/pfds/pfds_map_cont.html) type point precipitation data and recommended synthetic storms are typically not available or easily obtained in OCONUS locations. Use of a common source of precipitation data allowed for a standardized method of precipitation estimation for OCONUS sites. 

The tool extracts rainfall intensity data from global duration-frequency curves in the the Parameterized eXtreme Rain dataset (Courty et al. 2019) which are derived from ERA5 reanalysis data. 

**There are two datasets that the PET can work with, PXR2 and PXR4.** 

**PXR2** is a 31 km resolution dataset of GEV parameters (location ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\mu}) and scale ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\sigma}), the shape parameter ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\kappa}) is assumed constant globally) for a list of discrete durations. Intensity is a funtion of duration ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\d}) and return period ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\T}) , ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\i(d,T)=\mu_d%2B\sigma_dy})  where ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\y}) is a quantile function. The location and scale parameters are selected from a discrete list of values corresponding to durations ranging from 1 to 360 hours. 

**PXR4** is a 31 km resolution four parameter fit of the same data, but the location and scale parameters are assumed to be continuous functions of duration ![formula](https://render.githubusercontent.com/render/math?math=\bbox[white,2pt]{\large\i(d,T)={ad^\alpha%2Bbd^\beta}y})

*Note that in Courty et al. (2019) eqn. 7 contains a typographical error and the lone minus sign should be an addition sign. This was confirmed via correspondence with the lead author.

Caution should be used for using the PXR data for short duration storms (i.e. < 24-hr) because ERA5 may not fully capture the intensity of short duration convective storms- see Courty et al. (2019) for more detail. A limited comparison of the PXR data and NOAA Atlas 14 was conducted by CRREL in 2020, which showed that the PXR2 areal precitipation values tended to be less than the NOAA Atlas 14 point data.

The PET was developed using version 2.1.0 of the NetCDF files for PXR2 and PXR4 which were obtained from the following location: [Parametrized eXtreme Rain (PXR) | Zenodo](https://zenodo.org/record/3351812)

Courty, L. G., Wilby, R. L., Hillier, J. K., & Slater, L. J. (2019). Intensity-duration-frequency curves at the global scale. _Environmental Research Letters_, _14_(8), 084045. http://dx.doi.org/10.1088/1748-9326/ab370a

## Use
Disclaimer: This tool was developed to rapidly create precip boundary conditions for a large number of OCONUS sites as part of a fast moving study. While the use of the tool is relatively straightforward, it was not created with user experience as a key motivator. This documentation was created to describe the use of the PET as-is, basically as it stood when it was initially created in January/February 2020. Contributions to improve the user experience or anything else are welcome.

### Point estimates
At the most basic level, the PET can be used to extract precipitation depth as a function of duration and return period from the PXR datasets (PXR2, PXR4) at a location specified with decimal latitude and longitude.

The tool was developed to generate several types of hyetographs:

### 24-hr alternating block hyetograph 
Assembles an alternating block hyetograph (see Chow 1964 for details)
Uses durations available in the PXR2 dataset. 
The only variables are location and return period
It is currently hardcoded to be centered at the 12-hr mark
**This is the method used in previous climate RAS models**

Chow, V.T. (1964). Handbook of applied hydrology: a compendium of water-resources technology.

### Chicago Storm Method hyeotgraph
Uses the 4 fitted parameters from the PXR4 dataset to create a smooth hyetograph using the Chicago Storm method ([Synthetic Storm Pattern for Drainage Design | Journal of the Hydraulics Division | Vol 83, No 4 (ascelibrary.org)](https://ascelibrary.org/doi/abs/10.1061/JYCEAJ.0000104) [Theory - Derivation Of The Chicago Storm (alanasmith.com)](http://www.alanasmith.com/theory-Derivation-Chicago-Storm.htm))

### Hybrid PXR2 Depth and PXR4 Chicago Storm Distribution
This approach uses the depth magnitude from PXR2 for a specified storm duration, but applies it using a temporal distribution based on the PXR4 fit parameters. PXR2 depths tend to be higher then PXR4 depths for the 100-year 24-hr storm.

## Dependencies
PET is written in python 3
The tool uses the following libraries: **xarray, pydsstools, numpy, netcdf4, scipy**

the pydsstools library is required to export the data to DSS for use in HEC software.See [gyanz/pydsstools: Python library for simple HEC-DSS functions (github.com)](https://github.com/gyanz/pydsstools) for installation information. Note installing from source is not recommended (required MS Visual  C++ and Fortan compliers), the more reliable method is to install from one of the provided wheel files. See [Problems with the installation of the library · Issue #23 · gyanz/pydsstools (github.com)](https://github.com/gyanz/pydsstools/issues/23)

For computing hyetographs for multiple sites, there are two additional scripts that ingest watershed polygons and produce DSS databases. These scripts (get_centroids.py, multi_site_extract.py) require access to arcpy (and a spatial analyist licence) and the json library. Arcpy is used only to extract the centroids listed in the JSON file that feeds multi_site_extract.py - a user without access to arcpy could generate an input file of watershed centroid locations in the same format (example in file directory) and still run multi_site_extract.py 

## Workflows
### Single Site
The PET is point based and only needs a decimal latitude and longitude to extract data. The only script needed to do this is PET.py. PET can be imported into another script or used at the command line, but it may be sufficient to simply edit the main() function in PET.py as needed and run the script.

The user must specify the location of the NetCDF files containing the PXR data. The data can be obtained from: [Parametrized eXtreme Rain (PXR) | Zenodo](https://zenodo.org/record/3351812)

The user specifies:
 
**Duration in hours** e.g. 24

**Return period in years** e.g 100

**Time to peak in hours** e.g 12 (note- the alternating block hyetograph method is hardcoded to peak at 12 hours- Time to peak only affects the Chicago Storm Method)

**Site coordinates in decimal degrees** [lat, long]

**Installation Name** String that will be used for the output file

There are example workflows in the main() function of PET.py for the three hydrograph methods discussed above. Note that in the main() function the output file is the same for all three hyetograph options, so all three outputs will be written to the same DSS, but with different DSS "Part F"s which allow for differentiation within the database file.

### Multiple Sites
The PET is point based and only needs a decimal latitude and longitude to extract data. However, a user may want to process a number of watersheds at once. Some auxiliary scripts were written to facilitate batch processing of multiple watersheds.

#### Obtain centroids of AOI and write to JSON text file
The get_centroids.py script is a simple tool that uses arcpy to process a folder of waterhshed shapefiles and outputs a text file in JSON format. This tool requires arcpy and a Spatial Analyist license. The main() function gives an example of how to run the tool. 

Because this tool requires arcpy, it needs to run on the Python tha comes with ArcMap. It may be easiest to edit the main() function in an IDE to point to the right local paths for input/output, but run the script from within an ArcMap Python prompt using a command like:

execfile(r'c:\my_directory\get_centroids.py')

The JSON text file format for the centroid output is simply:

\[\["Installation1 Name", \[Centroid1 Latitude, Centroid1 Longitude]], 
\["Installation2 Name", \[Centroid2 Latitude, Centroid2 Longitude]],...] 

lat/lons are in decimal degrees. The user doesn't need to use get_centroids.py or arcpy to create this file, it could be done manually or with another script. It would be nice to make this a little less cumbersome/ not rely on arcpy, but it does work fine if you have ArcMap. 

#### Run PET.py tools in a loop with multi_site_extract.py
The multi_site_extact.py script simply parses installation names and lat/lon locations from a JSON formatted text file. The get_centroids.py script will generate this file from a folder of watershed shapefiles, but the text file can be generated in other ways. See the **Import_Locations.txt** file in the directory for an example of the file format. The multi_site_extract.py script simply runs the PET functions on each location on the list, creating a unique DSS file for each installation. The main() function is pretty self explanatory.

## Single Site Example
Generate a 24-hour 100-year alternating block hyetograph DSS file for use in HEC RAS at a location in Vermont, USA

Download the applicable PXR datasets (NetCDF files) from the location provided by Courty et al. (2019):

![PXR download page](https://github.com/ski907/PXRextract/blob/main/images/Pasted%20image%2020220110120144.png?raw=true)

![data_dir](https://github.com/ski907/PXRextract/blob/main/images/Pasted%20image%2020220110124344.png?raw=true)
Open PET.py in you IDE of choice. Update **data_dir** to point to the location of the pxr2 and pxr4 files.

![duration, return period](https://github.com/ski907/PXRextract/blob/main/images/Pasted%20image%2020220110125200.png?raw=true)
Check the **duration** and **T** (return period) values and adjust as necessary. **Tp** (time to peak) can be modified, but note for the alternating block hyetograph the peak is hardcoded to be centered at 12 hours, so modification to **Tp** won't do anything. 

![site coordinates](https://github.com/ski907/PXRextract/blob/main/images/Pasted%20image%2020220110125237.png?raw=true)
Update **site_coord** to the location of interest \[lat, lon] in decimal degrees
Update **installation** to a string identifying the site location. It will be used in the output file name

Run the script (comment out other methods below, i.e. Chicago, if desired)
The tool will import the NetCDF, extract the GEV fit parameters (location and scale) for the location specified, calculate the rainfall intensities for a range of durations at the specified return period, assemble a hyetograph using alternating blocks and export that to a DSS file. The units are metric. The output can be visualized in DSSVue or by importing it into HEC-RAS as a rainfall boundary condition.

![files in explorer](https://github.com/ski907/PXRextract/blob/main/images/Pasted%20image%2020220110134050.png?raw=true)

![DSSVue](https://github.com/ski907/PXRextract/blob/main/images/Pasted%20image%2020220110133946.png?raw=true)


If you have questions please feel free to contact me: Chandler.S.Engel@usace.army.mil
