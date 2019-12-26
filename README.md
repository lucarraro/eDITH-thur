# eDITH
MATLAB code and data supporting Carraro L, Mächler E, Wüthrich R, Altermatt F. Upscaling spatial patterns of biodiversity in freshwater ecosystems, *Nature Communications*, under review as of December 2019.

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
	- `Clipped_Geology.asc`: raster map of geological classes covering the area of interest (processed via ArcGIS, based on data from the Swiss Federal Institute of Topography).
	- `Covariates.mat`: covariate matrix (obtained via `additional_scripts/BuildCovariateMatrix.m`).
	- `DTM.asc`: digital terrain model of the area of interest (processed via ArcGIS, based on data from the Swiss Federal Institute of Topography).
	- `DTMad8.asc`: raster map of drainage areas (obtained via Taudem subroutine run on ArcGIS).
	- `DTMfel.asc`: pit-filled digital terrain model (obtained via Taudem subroutine run on ArcGIS).
	- `DTMp.asc`: flow direction matrix (obtained via Taudem subroutine run on ArcGIS).
	- `GenusData.m`: eDNA and kicknet data. Data obtained as reported in Mächler et al., 2019 <doi:10.1002/edn3.33>.
	- `landcover.asc`: raster map of land cover classes covering the area of interest (processed via ArcGIS, based on data from the Swiss Federal Institute of Topography).
	- `SitesFunWorks_coordinates.xlsx`: coordinates of sites where eDNA and kicknet sampling were performed.
	- `StageDischarge.xlsx`: stage-discharge relationships (provided by the Swiss Federal Office for the Environment).
	- `ThurData.m`: data on the catchment morphology and river network extraction (produced by `additional_scripts/CreateNetwork.m`).
	- `ThurHydrology.m`: data on hydrological variables across the catchment (produced by `additional_scripts/hydrology.m`).
	- `ThurQ.xlsx`: discharge time series (provided by the Swiss Federal Office for the Environment).
- `functions`: MATLAB functions required by the scripts.
	- `DrawRiverMap.m`: draws thematic river maps.
	- `eval_accuracy.m`: evaluates the accuracy of eDITH predictions with respect to kicknet data.
	- `eval_likelihood.m`: evaluates likelihood, given a set of parameters and read number data for a given genus.
	- `eval_model.m`: runs the eDITH model.
	- `GOF_NBtest.m`: performs the goodness-of-fit test.
	- `multiple_boxplot.m`: draws multiple boxplots on the same figure.
	- `neigh.m`: transforms D8 flow direction values in row/column movements.
	- `v2struct.m`: assembles and disassembles MATLAB structures (used to import data into functions).
- `results`: stores results produced by scripts `MAIN_MCMC.m` and `MAIN_MCMC_subsampling.m`.
	- `all`: results produced by `MAIN_MCMC.m`. Files' names correspond to the genera. 
	-  `subsampling_X_Y`: results produced by `MAIN_MCMC_subsampling.m`, where `X` is the number of validation sites and `Y` the run number. Files' names correspond to the genera.
	- `subsampling_results.m`: synthesis of results from `subsampling_X_Y`.
	