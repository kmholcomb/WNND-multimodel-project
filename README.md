## Prediction of West Nile virus disease by climate region in the United States (2015-2020)
This repository contains R scripts for fitting and scoring models to predict annual West Nile virus neuroinvasive disease (WNND) cases in climate regions of the United States (2015-2021). We also used machine learning models to identify important regional climate/weather factors for WNND prediction, using normalized anomalies in seasonal or monthly temperature and precipitation. We fitted ten models, ranging in complexity, to predict annual WNND cases across a hexagonal grid. 

### Repository structure
All included code can be run with R with the required packages indicated at the top of each script. Analyses were performed in RStudio (v. 2022.07.1+554) using R (v. 4.2.0) on a 64-bit Windows 10 laptop (16 GB RAM, 11th Gen Intel(R) Core(TM) i7-1185G7 @ 3.00GHz processor) or using R (v. 4.2.0) through a compute cluster environment (33 compute nodes with 44 cores each running at 2.8Ghz & 384 GB RAM each). 

All R scripts and an R project file are in the **R_files** folder. The R scripts will place processed files (.RDS files) in a separate **data** folder (you will need to manually create this folder). Machine learning models are in Rmarkdown (.RMD) files for visualization of fit/prediction metrics and variable importance.

### Data sources
The following data sources were used in the analysis and are publically available for download (you need to download them). Scripts for processesing these data to the hexagonal grid (and saving as .RDS files for use in model fitting) are included here.
- Human WNND case data are visible via CDC's ArboNET [disease map](https://wwwn.cdc.gov/arbonet/maps/ADB_Diseases_Map/index.html), and researchers can request it from ArboNET by emailing <dvbid2@cdc.gov>. The corresponding processing script is _case_demographics_hex.R_.
- Demographic data from the 2010 census are openly available from the [U.S. Census Bureau](https://www.census.gov/programs-surveys/decennial-census/data/datasets.2010.html) for download. The corresponding processing scripts are _census_data.R_ and _case_demographics_hex.R_.
- Monthly temperature (minimum and mean) and total precipitation data are openly available from [PRISM](https://www.prism.oregonstate.edu/recent/) for download. The corresponding processing script is _prism_temp_precip.R_.
- Land use data are openly available from the [Multi-Resolution Land Characteristics Consortium (MRLC)](https://www.mrlc.gov/data/nlcd-2011-land-cover-conus) for download. The corresponding processing script is _landcover_hex.R_.

### Hexagonal grid
We created a grid over the contiguous United States (equal area projection) using 200 km diameter hexagons (~34,640 km^2 area).

Code for making the grid and assigning climate regions to hexagons are in the _hex_grids.R_ script.

### Model fitting and scoring
Model fitting, prediction, and calculation of logarithmic scores for each of the models are found in the following scripts. Scripts often use multi-core computation, but comments indicate how to adjust for using a single core.
- Always Absent model: _AA_hex.R_
- Negative binomial models (hex, region, and national scales): _NB_models.R_
- Autoregressive model (AR(1)): _clim_ar1_hex.R_
- Autoregressive climate model (AR(1) Climate): _clim_ar1_hex.R_
- Machine learning models (random forest and neural network models): _NN_hex_month.RMD_, _NN_hex_season.RMD_, _RF_hex_month.RMD_, and _RF_hex_season.RMD_

Comparison of logarithmic scores between models is outlined in the _comparison_models_hex.R_ script.

Bayesian regression analyses used to statisitically identify differences in model performance are illustrated in the _regression_hex.R_ script.
