# SMAR_EUS_RemoteSensing
Spatial Root Zone Soil Moisture output from the Soil Moisture Analytical Relationship model for the Shale Hills Critical Zone Observatory, which is featured in a manuscript submitted to the journal Remote Sensing.  Example code for calibrating SMAR with NASA's AMSR-E near-surface satellite data and in situ root zone soil moisture data is also included.

SMAR output for the Eastern United States Region totals 68.9 GB, but I would be happy to provide this output upon request (email: dcb5006@gmail.com). 

(/SMAR_EUS_RemoteSensing/Figure8.JPG)

### Acronyms ###
EUS: Eastern United States
RZSM: Root Zone Soil Moisture (in units of relative soil moisture content: dimensionless)
SMAR: Soil Moisture Analytical Relationship
SCAN: Soil Climate Analysis Network: https://www.wcc.nrcs.usda.gov/scan/

######## Files ##########
1) SMAR Parameter Maps:
AMSRE_a_map.tif
AMSRE_b_map.tif
AMSRE_FC_map.tif
AMSRE_WL_map.tif
SMOS_a_map.tif
SMOS_b_map.tif
SMOS_FC_map.tif
SMOS_WL_map.tif
SMAP_a_map.tif
SMAP_b_map.tif
SMAP_FC_map.tif
SMAP_WL_map.tif
These are EUS regional maps of SMAR parameters that were created with neural network models.  The neural networks were trained on calibrated SMAR parameters associated with SMAR model-fitting using data from SCAN/SNOTEL across the EUS region.  Covariates include soil property maps (http://www.soilinfo.psu.edu/index.cgi?soil_data&conus&data_cov) and USGS Elevation (https://nationalmap.gov/small_scale/mld/elev100.html).  The maps match the CONUS soil property map resolution of 1 km.

2) SMAR_optim_example_ShaleHills.R:
This code inputs in situ RZSM data and AMSR-E near-surface satellite data (in the file: "Moisture_data.csv") as part of an optimization routine to calibrate the SMAR model. A Bayesian-MCMC approach is used for the optimization. MCMC chains are partially informed by estimating SMAR parameters with soil texture (provided by the file: "SH_RealtimeSite_Data.csv"). The code was written in R.


References: Baldwin et al. 2019. Estimating Root Zone Soil Moisture Across the Eastern United States with Passive Microwave Satellite Data and a Simple Hydrologic Model. Remote Sensing: (https://www.mdpi.com/2072-4292/11/17/2013)
