###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Compile weather data per hexagon
# Extract mean & min temperature, total precipitation by month and season; calculate long-term normals and use to normalize weather data (num SD away from mean)
# PRISM data download portal: https://www.prism.oregonstate.edu/recent/
# NOTE: at time of processing, Oct-Dec of 2021 data was still provisional
###

library(raster)
library(parallel)
library(sf)
library(dplyr)
library(tigris)

cts <- counties(state = unique(fips_codes$state[!fips_codes$state %in% c("AK","HI","PR","AS","GU","MP","UM","VI")]),
                resolution = "20m", class="sf") %>% st_transform(., "EPSG:5070")
h <- st_make_grid(cts, cellsize = 2e5, square=FALSE) #cell size in m^2
conus <- h[cts] #crop to US outline


#### Calculate Normals (1999-2021) - mean and SD
ym <- expand.grid(c("01","02","03","04","05","06","07","08","09","10","11","12"), 1999:2021) #set up object for looping through
ym <- rbind(c("12", 1998), ym) #add Dec of 1998 for seasonal calculations

#mean temperature
yr.stack_t <- lapply(1:nrow(ym), function(x) {
  if(ym$Var2[x] == 2021 & (ym$Var1[x] %in% c("10","11","12"))) {
    raster(paste('../data/PRISM_monthly/', ym[x,"Var2"],'/PRISM_tmean_provisional_4kmM3_', ym[x,"Var2"], ym[x,"Var1"],'_bil.bil',sep="")) 
  } else {
    raster(paste('../data/PRISM_monthly/', ym[x,"Var2"],'/PRISM_tmean_stable_4kmM3_', ym[x,"Var2"], ym[x,"Var1"],'_bil.bil',sep="")) 
  }
})
yr.stack_t <- stack(yr.stack_t) #stack list of rasters

#min temperature
yr.stack_tn <- lapply(1:nrow(ym), function(x) {
  if(ym$Var2[x] == 2021 & (ym$Var1[x] %in% c("10","11","12"))) {
    raster(paste('../data/PRISM_monthly/', ym[x,"Var2"],'/PRISM_tmin_provisional_4kmM3_', ym[x,"Var2"], ym[x,"Var1"],'_bil.bil',sep=""))
  } else {
    raster(paste('../data/PRISM_monthly/', ym[x,"Var2"],'/PRISM_tmin_stable_4kmM3_', ym[x,"Var2"], ym[x,"Var1"],'_bil.bil',sep="")) 
  }
})
yr.stack_tn <- stack(yr.stack_tn) #stack list of rasters

#total precipitation
yr.stack_p <- lapply(1:nrow(ym), function(x) {
  if(ym$Var2[x] == 2021 & (ym$Var1[x] %in% c("10","11","12"))) {
    raster(paste('../data/PRISM_monthly/', ym[x,"Var2"],'/PRISM_ppt_provisional_4kmM3_', ym[x,"Var2"], ym[x,"Var1"],'_bil.bil',sep=""))
  } else {
    raster(paste('../data/PRISM_monthly/', ym[x,"Var2"],'/PRISM_ppt_stable_4kmM3_', ym[x,"Var2"], ym[x,"Var1"],'_bil.bil',sep="")) 
  }
})
yr.stack_p <- stack(yr.stack_p) #stack list of rasters

tp_norm <- mclapply(1:length(conus), function(x) { #for each hexagon in CONUS, calc normal temp, precip
  
  print(round(x/length(conus), 3)) #monitoring progress through hexes (% of way done)
  
  ##normal calcs - mean of mean temp, mean of min temp, and mean of total precip
  #mean temperature calcs
  t <- raster::extract(yr.stack_t, st_as_sf(conus[x]), df=TRUE) #all temps for the year (Dec-Nov)
  print(1) #monitor progress
  
  #seasonal mean temperature (long term mean and SD)
  t.w <- c(tmean.w = mean(unlist(t[,head(grep("12_|01_|02_", names(t)),-1)]), na.rm =T), #mean of Dec-Feb (winters), not including Dec 2021
           tsd.w = sd(unlist(t[,head(grep("12_|01_|02_", names(t)),-1)]), na.rm =T)) #SD
  t.sp <- c(tmean.sp = mean(unlist(t[,grep("03_|04_|05_", names(t))]), na.rm =T), #mean of Mar-May (springs)
            tsd.sp = sd(unlist(t[,grep("03_|04_|05_", names(t))]), na.rm =T))
  t.su <- c(tmean.su = mean(unlist(t[,grep("06_|07_|08_", names(t))]), na.rm =T), #mean of June-Aug (summers)
            tsd.su = sd(unlist(t[,grep("06_|07_|08_", names(t))]), na.rm =T))
  t.f <- c(tmean.f = mean(unlist(t[,grep("09_|10_|11_", names(t))]), na.rm =T), #mean of Sep-Nov (falls)
           tsd.f = sd(unlist(t[,grep("09_|10_|11_", names(t))]), na.rm =T))
  
  #monthly mean temperature (long term mean and SD)
  t.m <- sapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(x) {
    if(x != "12") {
      mean(unlist(t[,grep(paste(x,"_",sep=""), names(t))]), na.rm =T) #mean on each month for the 22 year period
    } else {
      mean(unlist(t[,head(grep("12_", names(t)),-1)]), na.rm =T) #not Dec 1998)
    }
  })
  names(t.m) <- paste('tmean_', names(t.m), sep="") #more descriptive names
  
  tsd.m <- sapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(x) {
    if(x != "12") {
      sd(unlist(t[,grep(paste(x,"_",sep=""), names(t))]), na.rm =T)
    } else {
      sd(unlist(t[,head(grep("12_", names(t)),-1)]), na.rm =T)
    }
  })
  names(tsd.m) <- paste('tmeansd_', names(tsd.m), sep="") #more descriptive names
  
  #min temperature calcs
  tn <- raster::extract(yr.stack_tn, st_as_sf(conus[x]), df=TRUE) #all monthly min temps for all years (Jan-Dec)
  print(2) #monitor progress
  
  #seasonal min temperature (long term mean and SD)
  tn.w <- c(tmin.w = mean(unlist(tn[,head(grep("12_|01_|02_", names(tn)),-1)]), na.rm =T), #mean of Dec-Feb (winters), not including Dec 2020
            tnsd.w = sd(unlist(tn[,head(grep("12_|01_|02_", names(tn)),-1)]), na.rm =T)) #SD
  tn.sp <- c(tmin.sp = mean(unlist(tn[,grep("03_|04_|05_", names(tn))]), na.rm =T), #mean of Mar-May (springs)
             tnsd.sp = sd(unlist(tn[,grep("03_|04_|05_", names(tn))]), na.rm =T))
  tn.su <- c(tmin.su = mean(unlist(tn[,grep("06_|07_|08_", names(tn))]), na.rm =T), #mean of June-Aug (summers)
             tnsd.su = sd(unlist(tn[,grep("06_|07_|08_", names(tn))]), na.rm =T))
  tn.f <- c(tmin.f = mean(unlist(tn[,grep("09_|10_|11_", names(tn))]), na.rm =T), #mean of Sep-Nov (falls)
            tsd.f = sd(unlist(tn[,grep("09_|10_|11_", names(tn))]), na.rm =T))
  
  #monthly min temperature (long term mean and SD)
  tn.m <- sapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(x) {
    if(x != "12") {
      mean(unlist(tn[,grep(paste(x,"_",sep=""), names(tn))]), na.rm =T) #mean on each month for the 22 year period
    } else {
      mean(unlist(tn[,head(grep("12_", names(tn)),-1)]), na.rm =T) #not Dec 1998)
    }
  })
  names(tn.m) <- paste('tmin_', names(tn.m), sep="") #more descriptive names
  
  tn.sd <- sapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(x) {
    if(x != "12") {
      sd(unlist(tn[,grep(paste(x,"_",sep=""), names(tn))]), na.rm =T)
    } else {
      sd(unlist(tn[,head(grep("12_", names(tn)),-1)]), na.rm =T)
    }
  })
  names(tn.sd) <- paste('tminsd_', names(tn.sd), sep="") #more descriptive names
  
  #total precipitation calcs
  p <- raster::extract(yr.stack_p, st_as_sf(conus[x]), df=TRUE) #all total precip for all years (Jan-Dec)
  print(3) #monitoring progress
  
  #seasonal total precipitation (long term mean and SD)
  ppt.w <- sapply(0:21, function(x, pa = p[,grep("12_|01_|02_", names(p))]) { #seasonal total (winter - Dec-Feb)
    sum(colMeans(pa[,c((1+(3*x)):(3+(3*x)))], na.rm = T))
  })
  ppt.sp <- sapply(0:21, function(x, pa = p[,grep("03_|04_|05_", names(p))]) { #seasonal total (spring - Mar-May])
    sum(colMeans(pa[,c((1+(3*x)):(3+(3*x)))], na.rm = T))
  })
  ppt.su <- sapply(0:21, function(x, pa = p[,grep("06_|07_|08_", names(p))]) { #seasonal total (summer - Jun-Aug)
    sum(colMeans(pa[,c((1+(3*x)):(3+(3*x)))], na.rm = T))
  })
  ppt.f <- sapply(0:21, function(x, pa = p[,grep("09_|10_|11_", names(p))]) { #seasonal total (fall - Sep-Nov)
    sum(colMeans(pa[,c((1+(3*x)):(3+(3*x)))], na.rm = T))
  })
  
  #monthly total precipitation (long term mean and SD)
  ppt.m <- sapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(x) {
    if(x != "12") {
      mean(unlist(colMeans(p[,grep(paste(x,"_",sep=""), names(p))], na.rm =T))) #mean for tot precip/month the 22 year period
    } else {
      mean(unlist(colMeans(p[,head(grep("12_", names(p)),-1)], na.rm =T))) #not Dec 1998
    }
  })
  names(ppt.m) <- paste('ppt_', names(ppt.m), sep="") #more descriptive names
  
  pptsd.m <- sapply(c("01","02","03","04","05","06","07","08","09","10","11","12"), function(x) {
    if(x != "12") {
      sd(unlist(colMeans(p[,grep(paste(x,"_",sep=""), names(p))], na.rm =T))) #mean for tot precip/month the 22 year period
    } else {
      sd(unlist(colMeans(p[,head(grep("12_", names(p)),-1)], na.rm =T))) #not Dec 1998
    }
  })
  names(pptsd.m) <- paste('pptsd_', names(pptsd.m), sep="") #more descriptive names
  
  #all together into single data frame
  data.frame(fips = x, t(t.w), t(t.sp), t(t.su), t(t.f), t(t.m), t(tsd.m), 
             t(tn.w), t(tn.sp), t(tn.su), t(tn.f), t(tn.m), t(tn.sd), 
             ppt.w = mean(ppt.w), pptsd.w = sd(ppt.w), ppt.sp = mean(ppt.sp), pptsd.sp = sd(ppt.sp),
             ppt.su = mean(ppt.su), pptsd.su = sd(ppt.su), ppt.f = mean(ppt.f), pptsd.f = sd(ppt.f),
             t(ppt.m), t(pptsd.m))
  
}, mc.cores = detectCores())
saveRDS(do.call(rbind.data.frame, tp_norm), "../data/temp.precip_norm_hex.rds")


### Anomalies in seasonal and monthly amounts
# seasonal and monthly mean temp, min temp, and total precip
# NOTE: 2021 has some provisional data at time of processing
tp_ms <- lapply(1999:2021, function(x) { #for each year
  
  yr.stack_t <- if(x != 2021) { #mean temperature for that year
    stack(raster(paste('../data/PRISM_monthly/', x-1,'/PRISM_tmean_stable_4kmM3_', x-1, '12_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '01_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '02_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '03_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '04_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '05_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '06_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '07_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '08_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '09_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '10_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '11_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '12_bil.bil',sep="")))
    
  } else {
    stack(raster(paste('../data/PRISM_monthly/', x-1,'/PRISM_tmean_stable_4kmM3_', x-1, '12_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '01_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '02_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '03_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '04_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '05_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '06_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '07_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '08_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_stable_4kmM3_', x, '09_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_provisional_4kmM3_', x, '10_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_provisional_4kmM3_', x, '11_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmean_provisional_4kmM3_', x, '12_bil.bil',sep="")))
  }
  
  yr.stack_p <- if(x != 2021) { #total precipitation for that year
    stack(raster(paste('../data/PRISM_monthly/', x-1,'/PRISM_ppt_stable_4kmM3_', x-1, '12_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '01_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '02_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '03_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '04_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '05_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '06_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '07_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '08_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '09_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '10_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '11_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '12_bil.bil',sep="")))
  } else {
    stack(raster(paste('../data/PRISM_monthly/', x-1,'/PRISM_ppt_stable_4kmM3_', x-1, '12_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '01_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '02_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '03_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '04_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '05_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '06_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '07_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '08_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_stable_4kmM3_', x, '09_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_provisional_4kmM3_', x, '10_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_provisional_4kmM3_', x, '11_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_ppt_provisional_4kmM3_', x, '12_bil.bil',sep="")))
  }
  
  yr.stack_t2 <- if(x != 2021) { #min temperature for that year
    stack(raster(paste('../data/PRISM_monthly/', x-1,'/PRISM_tmin_stable_4kmM3_', x-1, '12_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '01_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '02_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '03_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '04_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '05_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '06_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '07_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '08_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '09_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '10_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '11_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '12_bil.bil',sep="")))
  } else {
    stack(raster(paste('../data/PRISM_monthly/', x-1,'/PRISM_tmin_stable_4kmM3_', x-1, '12_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '01_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '02_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '03_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '04_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '05_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '06_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '07_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '08_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_stable_4kmM3_', x, '09_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_provisional_4kmM3_', x, '10_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_provisional_4kmM3_', x, '11_bil.bil',sep="")),
          raster(paste('../data/PRISM_monthly/', x,'/PRISM_tmin_provisional_4kmM3_', x, '12_bil.bil',sep="")))
  }
  
  print(x) #monitor progress
  
  conus <- st_transform(conus, st_crs(yr.stack_t)) #same CRS as raster layer for extracting data
  
  g <- lapply(1:length(conus), function(z) { #calculate seasonal and monthly weather data for that year for each hexagon
    print(round(z/length(conus), 3)) #prop through hexes
    
    t <- raster::extract(yr.stack_t, st_as_sf(conus[z]), df=TRUE) #all mean temps for the year (Dec-Nov)
    tmean.w <- mean(unlist(t[,2:4]), na.rm =T) #mean of Dec-Feb (winter)
    tmean.sp <- mean(unlist(t[,5:7]), na.rm =T) #mean of Mar-May (spring)
    tmean.su <- mean(unlist(t[,8:10]), na.rm =T) #mean of June-Aug (summer)
    tmean.f <- mean(unlist(t[,11:13]), na.rm =T) #mean of Sep-Nov (fall)
    
    tmean.m <- colMeans(t[,3:14], na.rm = TRUE) #mean per month of current year (col 1 is ID col, col 2 is prev Dec)
    nms = stringr::str_match(names(tmean.m), "_([0-9]+)_")[,2] #pull out yearmonth from name
    names(tmean.m) <- paste('tmean_', substr(nms,5,6), sep="") #conciser names
    
    t2 <- raster::extract(yr.stack_t2, st_as_sf(conus[z]), df=TRUE) #all min temps for the year (Dec-Nov)
    tmin.w <- mean(unlist(t2[,2:4]), na.rm =T) #mean min of Dec-Feb (winter)
    tmin.sp <- mean(unlist(t2[,5:7]), na.rm =T) #mean min of Mar-May (spring)
    tmin.su <- mean(unlist(t2[,8:10]), na.rm =T) #mean min of June-Aug (summer)
    tmin.f <- mean(unlist(t2[,11:13]), na.rm =T) #mean min of Sep-Nov (fall)
    
    tmin.m <- colMeans(t2[,3:14], na.rm = TRUE) #mean per month of current year (col 1 is ID col, col 2 is prev Dec)
    nms = stringr::str_match(names(tmin.m), "_([0-9]+)_")[,2] #pull out yearmonth from name
    names(tmin.m) <- paste('tmin_', substr(nms,5,6), sep="") #conciser names
    
    p <- raster::extract(yr.stack_p, st_as_sf(conus[z]), df=TRUE) #all precip for the year (Dec-Nov)
    p_tot.w <- sum(colMeans(p[,2:4], na.rm = T)) #total precip accumulation (winter, DJF)
    p_tot.sp <- sum(colMeans(p[,5:7], na.rm = T)) #total precip accumulation (spring, MAM)
    p_tot.su <- sum(colMeans(p[,8:10], na.rm = T)) #total precip accumulation (summer, JJA)
    p_tot.f <- sum(colMeans(p[,11:13], na.rm = T)) #total precip accumulation (fall, SON)
    
    ppt.m <- colMeans(p[,3:14], na.rm = T) #monthly total precip (mean for all grid cells per county), (col 1 is ID col, col 2 is prev Dec)
    nms = stringr::str_match(names(ppt.m), "_([0-9]+)_")[,2] #pull out yearmonth from name
    names(ppt.m) <- paste('ppt_', substr(nms,5,6), sep="") #conciser names
    
    data.frame(fips = z, year = x, 
               t(tmean.m), tmean.w = tmean.w, tmean.sp = tmean.sp, tmean.su = tmean.su, tmean.f = tmean.f,
               t(tmin.m), tmin.w = tmin.w, tmin.sp= tmin.sp, tmin.su = tmin.su, tmin.f = tmin.f,
               t(ppt.m), ppt.w = p_tot.w, ppt.sp = p_tot.sp, ppt.su = p_tot.su, ppt.f = p_tot.f)
  }) #, mc.cores = detectCores())
  
  do.call(rbind.data.frame, g) #all hexes for that year together
}) 
saveRDS(do.call(rbind.data.frame, tp_ms),"../data/temp.precip_ms_hex.rds") #all years into single df, saved as R object


## Lag variables by a year -> useful when fitting AR(1) Climate candidate models
#tp_ms <- readRDS("../data/temp.precip_ms_hex.rds")
lags <- tp_ms
lags$year <- lags$year + 1 #bump year by one

lnames <- names(lags)[grep("tmean|tmin|ppt", names(lags))]
names(lags)[grep("tmean|tmin|ppt", names(lags))] <- paste("l.1y", lnames, sep="_") #append l.1_y to vars (lag of 1 year)

tpa2 <- merge(tp_ms, lags, by=c("fips","year")) #merge back together
saveRDS(tpa2, "../data/temp.precip_ms_hex2.rds")


## Calculate normalized anomalies (subtract mean / SD) for month and season
dat.av <- readRDS("../data/temp.precip_norm_hex.rds") #long term mean

dat <- readRDS("../data/temp.precip_ms_hex2.rds") #monthly and seasonal data (current and previous year)
dat.m <- dat[,-grep("w|sp|su|.f", names(dat))] #monthly tmean, ppt per county (that year and previous)
dat.s <- dat[,grep("fips|year|w|sp|su|f", names(dat))] #seasonal mean tmean, ppt per county (that year and previous)

#month
dat.m <- merge(dat.m, dat.av[,c(1, grep("_", names(dat.av)))], by="fips") #merge month avs/SD
std.tmean <- (dat.m[,3:14] - dat.m[,75:86])/dat.m[,87:98] #(obs - mean)/SD for current year tmean
std.tmin <- (dat.m[,15:26] - dat.m[,99:110])/dat.m[,111:122] #(obs - mean)/SD for current year tmin
std.ppt <- (dat.m[,27:38] - dat.m[,123:134])/dat.m[,135:146] #(obs - mean)/SD for current year ppt
std.tmean_l <- (dat.m[,39:50] - dat.m[,75:86])/dat.m[,87:98] #(obs - mean)/SD for previous year tmin
std.tmin_l <- (dat.m[,51:62] - dat.m[,99:110])/dat.m[,111:122] #(obs - mean)/SD for previous year tmean
std.ppt_l <- (dat.m[,63:74] - dat.m[,123:134])/dat.m[,135:146] #(obs - mean)/SD for previous year ppt
std.dat.m <- cbind(dat.m[,c("fips","year")], std.tmean, std.tmin, std.ppt, std.tmean_l, std.tmin_l, std.ppt_l) #bind back together with fips and year
names(std.dat.m)[-c(1:2)] <- paste(stringr::str_split(names(std.dat.m)[-c(1:2)], ".x", simplify = T)[,1], "_std", sep="") #append for clarity that values standardized
saveRDS(std.dat.m, "../data/prism_temp.precip_m_std_hex.rds")

#season
dat.s <- merge(dat.s, dat.av[,c(1,grep(".w|.sp|.su|.f", names(dat.av)))], by="fips") #merge seasonal avs/SD
std.tmean <- (dat.s[,3:6] - dat.s[,c(27,29,31,33)])/dat.s[,c(28,30,32,34)] #(obs - mean)/SD for current year tmean
std.tmin <- (dat.s[,7:10] - dat.s[,c(35,37,39,41)])/dat.s[,c(36,38,40,42)] #(obs - mean)/SD for current year tmin
std.ppt <- (dat.s[,11:14] - dat.s[,c(43,45,47,49)])/dat.s[,c(44,46,48,50)] #(obs - mean)/SD for current year ppt
std.tmean_l <- (dat.s[,15:18] - dat.s[,c(27,29,31,33)])/dat.s[,c(28,30,32,34)] #(obs - mean)/SD for previous year tmean
std.tmin_l <- (dat.s[,19:22] - dat.s[,c(35,37,39,41)])/dat.s[,c(36,38,40,42)] #(obs - mean)/SD for previous year tmin
std.ppt_l <- (dat.s[,23:26] - dat.s[,c(43,45,47,49)])/dat.s[,c(44,46,48,50)] #(obs - mean)/SD for previous year ppt
std.dat.s <- cbind(dat.s[,c("fips","year")], std.tmean, std.tmin, std.ppt, std.tmean_l, std.tmin_l, std.ppt_l) #bind back together with fips and year
names(std.dat.s)[-c(1:2)] <- paste(stringr::str_split(names(std.dat.s)[-c(1:2)], ".x", simplify = T)[,1], "_std", sep="") #append for clarity that values standardized
saveRDS(std.dat.s, "../data/prism_temp.precip_s_std_hex.rds")
