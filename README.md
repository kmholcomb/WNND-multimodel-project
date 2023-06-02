## Prediction of West Nile virus disease by climate region

This repository contains R scripts for fitting and scoring models to predict West Nile virus neuroinvasive disease (WNND) cases in climate regions of the United States (2015-2021). We fitted ten models, ranging in complexity, to predict WNND cases across a hexagonal grid

### Model fitting and scoring
Models
AA
NB
AR1
ML

Scoring

### Data sources
The following data sources were used in the analysis. Scripts for processesing these data to the hexagonal grid are also included.
- Human WNND case data are visible via CDC's ArboNET [disease map](https://wwwn.cdc.gov/arbonet/maps/ADB_Diseases_Map/index.html), and researchers can request it from ArboNET by emailing <dvbid2@cdc.gov>. 
- Demographic data from the 2010 census are openly available from the [U.S. Census Bureau](https://www.census.gov/programs-surveys/decennial-census/data/datasets.2010.html) for download. 
- Monthly temperature (minimum and mean) and total precipitation data are openly available from [PRISM](https://www.prism.oregonstate.edu/recent/) for download. 
- Land use data are openly available from the [Multi-Resolution Land Characteristics Consortium (MRLC)](https://www.mrlc.gov/data/nlcd-2011-land-cover-conus) for download.
