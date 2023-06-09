###
# Author: Karen Holcomb
# Last Updated: June 8, 2023
##
# Fit AR(1) and AR(1) Climate candidate models
# Final AR(1) Climate model per region selected based on mean logarithmic score across prediction years (2015-2021)
###

library(tidyr)
library(dplyr)
library(forecast)
library(parallel)

## load in data
dat.clim <- readRDS("prism_temp.precip_s_std_hex.rds") #seasonal climate data on hex scale - normalized anomalies for current and previous year
df_all <- readRDS("../data/WNND_hex_2005-2021.rds") #WNND cases per year per hex
h_clim <- readRDS("h_clim_proj.rds") #climate regions

data <- merge(df_all, h_clim %>% as.data.frame() %>% dplyr::select(-geometry), by=c("fips")) #add in climate region

data <- merge(data, dat.clim, by = c("fips","year")) %>% #case data info and climate info
  mutate(ln_cases = log(tot_count + 1)) %>% #log count of cases per hex
  as.data.frame() %>% drop_na(ln_cases) #rm no data hexes


## function for quantile method of calculating logarithmic score per model
q23 <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99) #23 quantiles to calculate from models
quant_to_nbinom <- function(quants, values) {
  if (any(values < 0)) stop('Values for negatitve binomial must be >= 0')
  if (any(duplicated(quants))) stop('Repeated quantiles')
  # select reasonable start values
  start <- c(size = 1, mu = mean(values))
  
  nbinom_optim_fn <- function(par, quants, values) {
    tryCatch(
      { 
        quants_theo_hi <- pnbinom(values, size=par[[1]], mu=par[[2]])
        quants_theo_lo <- pnbinom(values - 1, size=par[[1]], mu=par[[2]])
        sum((pmin(quants_theo_hi, quants) - pmax(quants_theo_lo, quants))^2)
      },
      warning = function(w) {
        NaN
      }
    )
  }
  opt_par <- optim(par = start, fn = nbinom_optim_fn, quants=quants, values=values)$par
  return(list(mu=opt_par[['mu']], size=opt_par[['size']])) #returns optimized 'mu' and 'size' parameters
}


## loop through hexes and years to fit AR(1) and AR(1) Climate candidate models, calculate logarithmic scores, and output score and coefficients
# Note: if not using parallel package, change mclapply() to lapply() and remove mc.cores argument
AR_clim_hex <- mclapply(unique(data$fips), function(x, n=1e5) {
  print(x) #monitor progress
  loc.dat <- data[data$fips == x,] #subset to single hex
  
  yr_loop <- lapply(2015:2021, function(z) { #for each prediction year
    
    ## Fit AR(1) and predict next year
    if(all(loc.dat$ln_cases[loc.dat$year < z] == 0)) { #if never any cases, algorithm doesn't work so use NAs
      ar1 <- NA #ar1 coefficient NA since not fit
      b0 <- NA #intercept/mean NA since not fit
      xreg = NA #coef for clim variable NA since not fit
      ests_lar1 <- rep(0, n) #all 0s for never cases hexes
      
    } else { #if at least 1 case has previously been reported, can fit an Ar(1) model using Arima() from forecast package
      fit_l <- Arima(loc.dat$ln_cases[loc.dat$year < z], order = c(1,0,0), method="ML") 
      
      if(class(fit_l) == "Arima") {
        ar1 <- coef(fit_l)[1] #ar1 coefficient
        b0 <- coef(fit_l)[2] #intercept/mean
        xreg <- NA #climate var
        
        pred = predict(fit_l, n.ahead = 1) #predict next year
        ests_lar1 = rnorm(n, mean = as.numeric(pred$pred), 
                          sd = as.numeric(pred$se)) + sample(as.numeric(fit_l$residuals), n, replace = TRUE) #bootstrap samples based on normal approx + training residuals
        ests_lar1 = exp(ests_lar1) - 1 #log counts back to raw counts for calculating log scores
        ests_lar1[ests_lar1 < 0] = 0 #change any negative predictions to be 0
      } else { #propagate NAs from above
        ar1 <- NA #ar1 coef
        b0 <- NA #intercept/mean
        xreg <- NA
        ests_lar1 <- NA 
      }
    }
    
    # Calculate logscore: quantiles -> NB -> logscore
    fit_pars <- quant_to_nbinom(quants = q23, values = quantile(ests_lar1, q23)) #parameters to fit NB based on quantiles
    o_count <- loc.dat$tot_count[loc.dat$year == z] #observed number of cases
    
    p_count <- pnbinom(o_count, mu = fit_pars[["mu"]], size = fit_pars[["size"]], lower.tail = TRUE) - #prob of less than equal to obs cases
      pnbinom(o_count-1, mu = fit_pars[["mu"]], size = fit_pars[["size"]], lower.tail = TRUE) #prob of less than eqaul to obs cases - 1
    logscore_nb <-  log(p_count) #probability to logscore
    
    # output results for AR(1) model - hex, year, climate region, model name, logscore, coefficients
    e.log <- data.frame(location = x, clim.regn = unique(loc.dat$clim.reg), yr.pred = z, mod = "ar1", 
                        score_nb = logscore_nb, xreg = xreg, ar1 = ar1, b0 = b0, row.names = NULL) 
    
    
    ## Fit candidate AR(1) Climate models for that hex and year
    # each candidate model uses a single climate covariate (normalized anomalies in seasonal mean temp, min temp, total precip) for current and previous year (1-year lag)
    cvars <- data[data$fips == x & data$year == z,]
    e2.log <- lapply(1:length(grep("tmean|ppt|tmin", names(cvars))), function(y) { #for each column with a climate covariate
      if(all(loc.dat$ln_cases[loc.dat$year < z] == 0)) { #if never any cases, algorithm doesn't work so use NAs
        ar1 <- NA #ar1 coefficient NA since not fit
        b0 <- NA #intercept/mean NA since not fit
        xreg = NA #coef for climate variable NA since not fit
        ests <- rep(0, n)
      } else {
        fit <- Arima(loc.dat$ln_cases[loc.dat$year < z], order = c(1,0,0), xreg = loc.dat[loc.dat$year < z, c(4+y)], 
                     method="ML") #climate covariates start in the 5th column of the dataset
        if(class(fit) == "Arima") {
          ar1 <- coef(fit)[1] # ar1 coefficient
          b0 <- coef(fit)[2] # intercept/mean
          xreg = coef(fit)[3] # coefficient for climate varariable (Note: on log(case+1) scale)
          
          pred = predict(fit, n.ahead = 1, newxreg = cvars[,c(4+y)]) #predict next year using corresponding climate covariate value
          ests = rnorm(n, mean = as.numeric(pred$pred), 
                       sd = as.numeric(pred$se)) + sample(as.numeric(fit_l$residuals), n, replace = TRUE) #bootstrap samples based on normal approx + training residuals
          ests = exp(ests) - 1 #logscale back to raw counts
          ests[ests < 0] = 0 #change any negative predictions to be 0
        } else { #propagate NAs from above
          ar1 <- NA #ar1 coef
          b0 <- NA #intercept/mean
          preg <- NA
          ests <- NA
        }
      }
      
      # Calculate logscore for each AR(1) Climate models: quantiles -> NB -> logscore
      fit_pars <- quant_to_nbinom(quants = q23, values = quantile(ests, q23)) #parameters to fit NB based on quantiles
      o_count <- loc.dat$tot_count[loc.dat$year == z]
      
      p_count <- pnbinom(o_count, mu = fit_pars[["mu"]], size = fit_pars[["size"]], lower.tail = TRUE) - #prob of less than equal to obs cases
        pnbinom(o_count-1, mu = fit_pars[["mu"]], size = fit_pars[["size"]], lower.tail = TRUE) #prob of less than equal to obs cases - 1
      logscore_nb <-  log(p_count) #probability to logarithmic score
      
      # output results for candidate AR(1) Climate model - hex, year, climate region, 
      #   model name (w/ climate covariate), logscore, coefficients
      data.frame(location = x, clim.regn = unique(loc.dat$clim.reg), yr.pred = z, mod = paste("ar1", names(cvars)[7+y], sep="-"), 
                 score_nb = logscore_nb, xreg = xreg, ar1 = ar1, b0 = b0, row.names = NULL) 
    })
    
    arc_log <- do.call(rbind.data.frame, e2.log) #data frame of AR(1) Climate candidate models
    
    rbind.data.frame(e.log, arc_log) #put AR(1) and AR(1) Climate candidate models for that hex and year together
  })
  
  do.call(rbind.data.frame, yr_loop) #put all years together
  
}, mc.cores = 5)

ar_mods <- do.call(rbind.data.frame, AR_clim_hex)

saveRDS(ar_mods, "../data/AR_clim_allcandidates.rds") #list into a data frame for it all

ar1 <- ar_mods %>% subset(mod == "ar1") %>% mutate(fips = as.numeric(as.character(fips))) #just AR(1) models
saveRDS(ar1, "../data/ar1_ls.rds") #logarithmic scores for AR(1) models

arC <- ar_mods %>% subset(mod != "ar1") %>% #AR(1) Climate candidate models
  mutate(mscore_nb = pmax(-10, mscore_nb)) %>% #truncate to -10 for calculating best log score
  group_by(clim.regn) %>% summarise(mscore = mean(mscore_nb)) %>% #mean logarithmic score across all years for that climate region
  arrange(desc(mscore_nb), .by_group = TRUE) %>% slice_head(n = 1) %>% #keep single best scoring model
  right_join(., ar_mods, by = c("clim.regn", "mod")) #get rows from all candidate models for just the best Climate models

saveRDS(arC %>% mutate(fips = as.numeric(as.character(fips))), "../data/AR_clim_final.rds") #AR(1) Climate models
