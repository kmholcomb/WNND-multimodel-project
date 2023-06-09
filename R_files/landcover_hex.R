###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Extract proportion of each county in different categories
# Using 2011 National Land Cover Database: https://www.mrlc.gov/data/nlcd-2011-land-cover-conus
# Legend: https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
###

library(tidyr)
library(dplyr)
library(sf)
library(tigris)

lc <- raster::raster("~/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img") #image to raster
cts <- tigris::counties(state = unique(fips_codes$state[!fips_codes$state %in% c("AK","HI","PR","AS","GU","MP","UM","VI")]),
                        resolution = "20m", class="sf") %>% st_transform(., st_crs(lc)) #same proj as lc raster
h <- st_make_grid(cts, cellsize = 2e5, square=FALSE) #hexagonal grid, cell size in m^2
h_cts <- h[cts] #crop to US outline


library(parallel) #if not parallelizing, change 'mclapply' below to 'lapply' and remove 'mc.cores' argument
lu.prop <- mclapply(length(h_cts), function(x) { #for each hexagon
  print(x) #monitoring progress through loop
  lu <- raster::extract(lc, st_as_sf(h_cts[x])) #extract land use codes (see legend) for that hex
  data.frame(fips = x, urb = sum(lu[[1]] %in% c(21,22,23,24)) / length(lu[[1]]) * 100, #proportion of grid cells in hex that are developed (open space to high intensity)
             crops = sum(lu[[1]] == 82) / length(lu[[1]]) * 100, #prop of grid cells in hex that are cultivated crops
             wetlands = sum(lu[[1]] == 95) / length(lu[[1]]) * 100) #prop of grid cells with are wetlands
  }, mc.cores = 48)

lu.df <- do.call(rbind.data.frame, lu.prop)

saveRDS(lu.df, "../data/lu.df_hex.rds") #R object with processes land use data per hex
