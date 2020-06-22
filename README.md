# eDITH
MATLAB code and data supporting Carraro L, Mächler E, Wüthrich R, Altermatt F. Environmental DNA allows upscaling spatial patterns of biodiversity in freshwater ecosystems, *Nature Communications* (accepted), 2020.

## Main scripts
- `ANALYZE_DATA.m`: analyzes data and produces figures and tables shown in the manuscript and extended data.
- `MAIN_MCMC.m`: calibrates the eDITH model for all 50 genera with data from all 61 sites. Results from this script are stored in `results/all`.
- `MAIN_MCMC_subsampling.m`: calibrates the eDITH model for all 50 genera with data from subsets of the overall sampling sites. Results from this script are stored in `results/subsampling_...`.

## Folders
- `additional_scripts`: contains additional scripts called by the main scripts.
	- `BuildCovariateMatrix.m`: evaluates morphological, geological and land cover covariates. Its results are stored in `data/Covariates.m`.
	- `cmap_DTM.m`: defines terrain colormap.
	- `CreateNetwork.m`: extracts the river network and calculates morphological properties. Its results are stored in `results/ThurData.m`.
	- `hydrology.m`: calculates hydrological properties of the catchment. Its results are stores in `results/ThurHydrology.m`. 
	- `RegioCov.m`: evaluates geographical covariates.
- `data`: contains data required by the scripts.
	- `Clipped_Geology.asc`: raster map of geological classes covering the area of interest (processed via ArcGIS, based on data from the Swiss Federal Institute of Topography, accessible at https://www.geocat.ch/geonetwork/srv/eng/md.viewer#/full_view/ca917a71-dcc9-44b6-8804-823c694be516/tab/complete).
	- `colZissou.mat`: colormap for Figure 5.
	- `Covariates.mat`: covariate matrix (obtained via `additional_scripts/BuildCovariateMatrix.m`).
	- `DTM.asc`: digital terrain model of the area of interest (processed via ArcGIS, based on data from the Swiss Federal Institute of Topography: https://www.geocat.ch/geonetwork/srv/eng/md.viewer#/full_view/7f6ba579-1917-4868-966a-ba91ac157612/tab/complete).
	- `DTMad8.asc`: raster map of drainage areas (obtained via Taudem subroutine run on ArcGIS).
	- `DTMfel.asc`: pit-filled digital terrain model (obtained via Taudem subroutine run on ArcGIS).
	- `DTMp.asc`: flow direction matrix (obtained via Taudem subroutine run on ArcGIS).
	- `GenusData.m`: eDNA and kicknet data. Data obtained as reported in Mächler et al., 2019 https://onlinelibrary.wiley.com/doi/pdf/10.1002/edn3.33.
	- `landcover.asc`: raster map of land cover classes covering the area of interest (processed via ArcGIS, based on data from the Swiss Federal Institute of Topography, accessible at https://www.geocat.ch/geonetwork/srv/eng/md.viewer#/full_view/cfbd4793-4225-4743-942b-d9b97acfbfcc/tab/complete).
	- `SitesFunWorks_coordinates.xlsx`: coordinates of sites where eDNA and kicknet sampling were performed.
	- `StageDischarge.xlsx`: stage-discharge relationships (provided by the Swiss Federal Office for the Environment, available at https://www.hydrodaten.admin.ch/en/2181.html; https://www.hydrodaten.admin.ch/en/2303.html; https://www.hydrodaten.admin.ch/en/2305.html; https://www.hydrodaten.admin.ch/en/2374.html; https://www.hydrodaten.admin.ch/en/2414.html).
	- `ThurData.m`: data on the catchment morphology and river network extraction (produced by `additional_scripts/CreateNetwork.m`).
	- `ThurHydrology.m`: data on hydrological variables across the catchment (produced by `additional_scripts/hydrology.m`).
	- `ThurQ_jun2016.xlsx`: discharge time series (provided by the Swiss Federal Office for the Environment, available at https://www.hydrodaten.admin.ch/en/2181.html; https://www.hydrodaten.admin.ch/en/2303.html; https://www.hydrodaten.admin.ch/en/2305.html; https://www.hydrodaten.admin.ch/en/2374.html; https://www.hydrodaten.admin.ch/en/2414.html).
- `functions`: MATLAB functions required by the scripts.
	- `DrawRiverMap.m`: draws thematic river maps.
	- `eval_accuracy.m`: evaluates the accuracy of eDITH predictions with respect to kicknet data.
	- `eval_likelihood.m`: evaluates likelihood, given a set of parameters and read number data for a given genus.
	- `eval_model.m`: runs the eDITH model.
	- `GOF_NBtest.m`: performs the goodness-of-fit test.
	- `multiple_boxplot.m`: draws multiple boxplots on the same figure. Copied from: Ander Biguri (2020). multiple_boxplot.m https://www.mathworks.com/matlabcentral/fileexchange/47233-multiple_boxplot-m), MATLAB Central File Exchange. Retrieved June 22, 2020.
	- `neigh.m`: transforms D8 flow direction values in row/column movements.
	- `v2struct.m`: assembles and disassembles MATLAB structures (used to import data into functions). Copied from: Adi Navve (2020). Pack & Unpack variables to & from structures with enhanced functionality (https://www.mathworks.com/matlabcentral/fileexchange/31532-pack-unpack-variables-to-from-structures-with-enhanced-functionality), MATLAB Central File Exchange. Retrieved June 22, 2020.
- `results`: stores results produced by scripts `MAIN_MCMC.m` and `MAIN_MCMC_subsampling.m`.
	- `all`: results produced by `MAIN_MCMC.m`. Files' names correspond to the genera. 
	-  `subsampling_X_Y`: results produced by `MAIN_MCMC_subsampling.m`, where `X` is the number of validation sites and `Y` the run number. Files' names correspond to the genera.
	- `subsampling_results.m`: synthesis of results from `subsampling_X_Y`.
	