###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Fit negative binomial (NB) models as benchmark models - fit on hex (NB-hex), region (NB-region), and national (NB-nation) scales
# NB-hex and NB-nation fit with stan models (improved way to deal with places with 0 cases)
# NB-region fit using NB GLM
###

library(dplyr)
require(rstan) #Note, may need the development version from GitHub (not one on CRAN)
##remove.packages(c("rstan", "StanHeaders"))
##install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

#stan model code
stan_nb <-
  "
  data {
    int<lower=0> N;
    int<lower=0> y[N];
  }

  parameters {
    real<lower=0> mu;     // mean
    real<lower=0> invroot_phi;    // dispersion
  }

  model {
    mu ~ normal(0,1);
    invroot_phi ~ normal(0, 1);
    y ~ neg_binomial_2(mu, 1/sqrt(invroot_phi));
  }
  
  generated quantities {
    int<lower=0> y_pred;
    y_pred = neg_binomial_2_rng(mu, 1/sqrt(invroot_phi));
  }
"

df_all <- readRDS("../data/WNND_hex_2005-2021.rds") #counts per year by hex, with fips indicated already

##NB-hex: fit NB model per hex
#Note: this takes a while to run
nb_hex <- lapply(unique(df_all$fips), function(x) {
  loc.dat <- df_all[df_all$fips == x,]
  print(x)
  
  stan_fit <- stan(model_code = stan_nb, data = list(N = nrow(loc.dat[loc.dat$year < y,]), 
                                                     y = loc.dat$tot_count[loc.dat$year < y]),
                   chains = 4, iter = 1000) #initial fit so don't need to re-compile every time
  
  yr_loop <- lapply(2015:2021, function(y) { #for each prediction year used in analysis
      stan_fit <- stan(model_code = stan_nb, data = list(N = nrow(loc.dat[loc.dat$year < y,]),
                                                         y = loc.dat$tot_count[loc.dat$year < y]),
                       chains = 4, iter = 3000, warmup = 1500, fit = stan_fit, #use previous fit to prevent re-compiling
                       control = list(adapt_delta = 0.85)) #control to prevent divergent transitions
    
    ## Calculate logarithmic score for each model (ls = log(P(obs)))
      #prob assigned to obs number of cases (P(obs) = P(X<=obs) - P(X<=(obs-1))), using median value for parameters from Stan fit
    ls <- pnbinom(loc.dat$tot_count[loc.dat$year == y], size = 1/sqrt(quantile(rstan::extract(stan_fit, 'invroot_phi')[[1]], probs=0.5)), 
                  mu = quantile(rstan::extract(stan_fit, 'mu')[[1]], probs=0.5), lower.tail = T) - 
      pnbinom(loc.dat$tot_count[loc.dat$year == y] - 1, size = 1/sqrt(quantile(rstan::extract(stan_fit, 'invroot_phi')[[1]], probs=0.5)), 
              mu = quantile(rstan::extract(stan_fit, 'mu')[[1]], probs=0.5), lower.tail = T) 
    
    data.frame(fips = x, year = y, score_nb = log(ls))
  })
  do.call(rbind.data.frame, yr_loop)
})

nb_h <- do.call(rbind.data.frame, nb_hex) #log scores for NB-hex model per hex and prediction year


##NB-nation: fit NB model based on all hexagons in US
nb_nation <- do.call(rbind.data.frame, lapply(2015:2021, function(y) { #for each prediction year used in RF
  
  stan_fit <- stan(model_code = stan_nb, data = list(N = nrow(N = nrow(df_all[df_all$year < y,]), 
                                                              y = df_all$tot_count[df_all$year < y])),
                     chains = 4, iter = 3000, warmup = 1500, fit = stan_fit, #use previous fit to prevent re-compiling cores = 4, 
                     control = list(adapt_delta = 0.85))
  
  ## Calculate logarithmic score for each model (ls = log(P(obs)))
   #prob assigned to obs number of cases (P(obs) = P(X<=obs) - P(X<=(obs-1))), using median value for parameters from Stan fit
  p <- pnbinom(df_all$tot_count[df_all$year == y], size = 1/sqrt(quantile(rstan::extract(stan_fit, 'invroot_phi')[[1]], probs=0.5)), 
                mu = quantile(rstan::extract(stan_fit, 'mu')[[1]], probs=0.5), lower.tail = T) - 
    pnbinom(df_all$tot_count[df_all$year == y] - 1, size = 1/sqrt(quantile(rstan::extract(stan_fit, 'invroot_phi')[[1]], probs=0.5)), 
            mu = quantile(rstan::extract(stan_fit, 'mu')[[1]], probs=0.5), lower.tail = T) 
  
  data.frame(fips = unique(df_all$fips), year = y, score_nb_n = log(p))
}))


##NB-region: fit NB model using all hexagons in each region
df2 <- merge(df_all, h_clim %>% as.data.frame() %>% select(-geometry), by=c("fips")) #add in regions

nb_region <- do.call(rbind.data.frame, lapply(unique(df2$clim.regn), function(x) {
  loc.dat <- df2[df2$clim.regn == x & df2$year >= 2005,] %>% mutate(fips = factor(fips))
  
  yr.loop <- lapply(2015:2021, function(y) { #for each prediction year
    mod <- MASS::glm.nb(tot_count ~ fips, data = loc.dat[loc.dat$year < y,])
    pred <- predict(mod, newdata = data.frame(fips = unique(loc.dat$fips))) #fit mean per fips
    
    ests_nb <- sapply(1:length(pred), function(z) { #calc prob of obs case per fips, unique mean, same dispersal (theta)
      #prob of obs cases from fit distribution P(obs) - P(obs-1)
      pnbinom(loc.dat$tot_count[loc.dat$year == y & loc.dat$fips == unique(loc.dat$fips)[z]], size=mod$theta, mu=exp(pred[z]), lower.tail = T) -
        pnbinom(loc.dat$tot_count[loc.dat$year == y & loc.dat$fips == unique(loc.dat$fips)[z]] - 1, size=mod$theta, mu=exp(pred[z]), lower.tail = T)
    })
    
    data.frame(fips = unique(loc.dat$fips), clim.regn = x, year = y, score_nb_r = log(ests_nb))
  })
  do.call(rbind.data.frame, yr.loop)
}))

#save logscores from all NB models
nb_ls <- merge(nb_h, h_clim %>% as.data.frame() %>% select(-geometry), by=c("fips")) %>% #add in climate region
  left_join(., nb_region %>% mutate(fips = as.numeric(as.character(fips)))) %>%
  left_join(., nb_nation %>% mutate(fips = as.numeric(as.character(fips))))

saveRDS(nb_ls, "../data/NB_ls.rds") #save as R object so do not need to re-run
