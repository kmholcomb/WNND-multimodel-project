## Prediction of West Nile virus disease by climate region

This repository contains R scripts for fitting and scoring models to predict annual West Nile virus neuroinvasive disease (WNND) cases in climate regions of the United States (2015-2021). We fitted ten models, ranging in complexity, to predict WNND cases across a hexagonal grid

### Repository structure
All included code can be run with R. Analyses were performed using R (v. 4.2.0) through RStudio (v. 2022.07.1).

All R scripts and an R project file are in the **R_files** folder. There is also an empty **data** folder where processed files (.RDS) created by scripts are placed.

### Hexagonal grid
We created a grid over the contiguous United States (equal area projection) using 200 km diameter hexagons (~34,640 km2 area).

Code for making the grid and assigning climate regions to hexagons are in the hex_grids.R script.

### Model fitting and scoring
Model fitting, prediction, and calculation of logarithmic scores for each of the models are found in the following scripts.
- Always Absent model: AA_hex.R
- Negative binomial models (hex, region, and national scales): NB_models.R
- Autoregressive model (AR(1)): clim_ar1_hex.R
- Autoregressive climate model (AR(1) Climate): clim_ar1_hex.R
- Machine learning models (random forest and neural network models using anomalies in either seasonal or monthly climate variables): NN_hex_month.RMD, NN_hex_season.RMD, RF_hex_month.RMD, and RF_hex_season.RDM

Comparison of logarithmic scores between models is outlined in the comparison_models_hex.R script.

Bayesian regression analyses used to statisitically identify are illustrated in the sup_hex.R script.

### Data sources
The following data sources were used in the analysis. Scripts for processesing these data to the hexagonal grid are also included.
- Human WNND case data (2005-2021) are visible via CDC's ArboNET [disease map](https://wwwn.cdc.gov/arbonet/maps/ADB_Diseases_Map/index.html), and researchers can request it from ArboNET by emailing <dvbid2@cdc.gov>. The corresponding processing script is case_demographics_hex.R.
- Demographic data from the 2010 census are openly available from the [U.S. Census Bureau](https://www.census.gov/programs-surveys/decennial-census/data/datasets.2010.html) for download. The corresponding processing scripts are census_data.R and case_demographics_hex.R.
- Monthly temperature (minimum and mean) and total precipitation data are openly available from [PRISM](https://www.prism.oregonstate.edu/recent/) for download. The corresponding processing script is prism_temp_precip.R.
- Land use data are openly available from the [Multi-Resolution Land Characteristics Consortium (MRLC)](https://www.mrlc.gov/data/nlcd-2011-land-cover-conus) for download. The corresponding processing script is landcover_hex.R.
