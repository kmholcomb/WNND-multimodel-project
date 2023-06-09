###
# Author: Karen Holcomb
# Last Updated: June 2, 2023
##
# Comparison models for prediction of WNND cases: NB-hex, NB-region, NB-nation, AR1(1), AR(1) Climate, Always Absent
###

library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

## Logarithmic scores for all baseline/comparison models (calculated from other scripts)
nb_ls <- readRDS("../data/NB_ls.rds") #logscore for three NB models, from NB_models.R script
ar1 <- readRDS("../data/ar1_ls.rds") %>%
  select(-mod) %>% rename(ar1_nb = mscore_nb) #logscore for AR(1), from clim_ar1_hex.R
ar1clim <- readRDS("../data/ar_clim_final.rds") %>% 
  select(-mod) %>% rename(arC_nb = mscore_nb) #logscore for AR(1) Climate model, from clim_ar1_hex.R
aa <- readRDS("../data/AA_ls.rds") #logscore from Always Absent model, from AA_hex.R script
comps_ls <- left_join(nb_ls, ar1) %>% left_join(., ar1clim) %>% left_join(., aa) #all baseline comparison models together


## Logarithmic score from RF models
ls_preds_s <- readRDS("../Data/ls_preds_seasonal.rds") #seasonal RF
ls_preds_m <- readRDS("../Data/ls_preds_monthly.rds") #monthly RF

mls_s <- lapply(seq_along(ls_preds_s), function(x) {
  ls_preds_s[[x]] %>% group_by(location, horizon, clim.regn, year) %>% 
    summarise(mscore_RFs = mean(logscore_nb)) %>% 
    left_join(., ln_cases_country[[x]][,c("ln_cases","location","year","pdens","urb","crops","wetlands")] %>% 
                mutate(location = as.numeric(as.character(location))))
})
mls_m <- lapply(seq_along(ls_preds_m), function(x) {
  ls_preds_m[[x]] %>% group_by(location, horizon, clim.regn, year) %>% 
    summarise(mscore_RFm = mean(logscore_nb)) #mean logscore for each location's model (multiple windows) for each year predicted
})

mls_df_yr <- merge(do.call(rbind.data.frame, mls_s), do.call(rbind.data.frame, mls_m), 
                   by=c("location", "horizon", "clim.regn", "year"))


## Logarithmic score from NN models
ls_preds_nn.s <- readRDS("../data/NN_list_seasonal.rds")[["ls_preds"]] %>% rename(ls_nn.s_nb = logscore_nb) #seasonal NN
ls_preds_nn.m <- readRDS("../data/NN_list_month.rds")[["ls_preds"]] %>% rename(ls_nn.m_nb = logscore_nb) #monthly NN

ls_preds_nn <- left_join(ls_preds_nn.s, ls_preds_nn.m, by = c("year", "location", "clim.regn"))

ls_all <- mls_df_yr %>% left_join(., comps_ls, by=c("location"="fips", "year", "clim.regn")) %>%
  left_join(., ls_preds_nn, by=c("location", "year", "clim.regn")) 

saveRDS(ls_all, "../data/logscore_allmods.rds") #save logscores for all models to use in regression stan models

## Compare mean logscore per year for all models
ls_all2 <- ls_all %>% dplyr::select(location:year, score_nb, score_nb_r, score_nb_n, ar1_nb, arC_nb,
                                    score_aa, mscore_RFs, mscore_RFm, ls_nn.s_nb, ls_nn.m_nb) %>%
  pivot_longer(score_nb:ls_nn.m_nb, names_to="mod", values_to="logscore") %>%
  group_by(clim.regn, year, mod) %>% summarise(av_ls = mean(pmax(-10, logscore)))

l2 <- ls_all %>% dplyr::select(location:year, score_nb, score_nb_r, score_nb_n, ar1_nb, arC_nb,
                               score_aa, mscore_RFs, mscore_RFm, ls_nn.s_nb, ls_nn.m_nb) %>%
  pivot_longer(score_nb:ls_nn.m_nb, names_to="mod", values_to="logscore") %>%
  group_by(clim.regn, mod) %>% summarise(av_ls = mean(pmax(-10, logscore), na.rm=T))


## Make median ensemble - median of location probabilities per year
ens <- ls_all %>% dplyr::select(location, clim.regn:year, score_nb, score_nb_r, score_nb_n, ar1_nb, arC_nb,
                                mscore_RFs, mscore_RFm, ls_nn.s_nb, ls_nn.m_nb) %>% #all models, but the AA model
  pivot_longer(score_nb:ls_nn.m_nb, names_to="mod", values_to="loc_ls") %>%
  mutate(loc_prob = exp(loc_ls)) %>% #back to probability scale
  group_by(location, clim.regn, year) %>% summarise(ens_prob = median(loc_prob), #median probability for ensemble
                                                    ens_ls = log(ens_prob)) %>% #logscore for ensemble
  ungroup()

ens_yr <- ens %>% group_by(clim.regn, year) %>% summarise(av_ls = mean(pmax(-10, ens_ls))) %>% mutate(mod = "ens")

ens_all <- ens %>% group_by(clim.regn) %>% summarise(av_ls = mean(pmax(-10, ens_ls))) %>% mutate(mod = "ens")

ls_all2 <- rbind(ls_all2, ens_yr) #add in ensemble to dataframe for model scores for each year

l2 <- rbind(l2, ens_all) #add in ensemble to dataframe for overall model scores


## Plotting
ls_all2$clim.regn_plot <- as.numeric(factor(ls_all2$clim.regn, levels = c("Northwest", "NRockiesandPlains", "UpperMidwest", "OhioValley", "Northeast",
                                                                           "West", "Southwest","South", "Southeast")))
l2$clim.regn_plot <- as.numeric(factor(l2$clim.regn, levels = c("Northwest", "NRockiesandPlains", "UpperMidwest", "OhioValley", "Northeast",
                                                                          "West", "Southwest","South", "Southeast")))

region_full = c("1"="Northwest", "2"="N Rockies & Plains", "3"="Upper Midwest", "4"="Ohio Valley", "5"="Northeast", "6"="West", 
                "7"="Southwest", "8"="South", "9"="Southeast")

o <- l2 %>% group_by(clim.regn_plot) %>% arrange(desc(av_ls), .by_group = TRUE) %>% 
  mutate(ord = 1:n())
o_best = o[o$ord == 1,] %>% mutate(clim_mod = paste0(clim.regn, mod)) #looking at top models and relative scoring

ls_best <- ls_all2 %>% mutate(cr_mod = paste0(clim.regn, mod)) %>%
  subset(cr_mod %in% o_best$clim_mod) %>% select(-cr_mod) %>%
  mutate(mod = factor(mod, levels = c("score_nb_r", "arC_nb", "ar1_nb", "score_nb", "mscore_nbs", "ls_nn.s_nb", 
                                      "mscore_nbm" , "ls_nn.m_nb", "score_nb_u", "ens", "score_aa")))


ls_all2$mod <- factor(ls_all2$mod, levels = c("score_nb_r", "arC_nb", "ar1_nb", "score_nb", "mscore_nbs", "ls_nn.s_nb", 
                                             "mscore_nbm" , "ls_nn.m_nb", "score_nb_u", "ens", "score_aa"))

mypal <- c("#e6194B","#800000","#f58231","#808000","#3cb44b","#42d4f4","#469990","#000075","#911eb4","black","darkgrey")

linescore = ggplot(ls_all2[ls_all2$mod != "score_aa",]) + 
  geom_line(aes(x=year, y=av_ls, color = mod), linewidth = 0.3, alpha = 0.5) + #, linetype = "dashed"
  geom_line(data = ls_best, aes(x=year, y=av_ls, color = mod), linewidth = 0.7) +
  facet_wrap(~clim.regn_plot, nrow=2, labeller = labeller(clim.regn_plot=region_full)) +
  labs(y="mean logarithmic score per year", x="year") +
  scale_color_manual(values = mypal) +
  theme_cowplot(10) + panel_border(color = "black") + theme(legend.position = "none")

ggsave("meanlogscore_line.tiff", 
       plot_grid(linescore) +
         annotate("text", x = 0.82, y = 0.46, label = "Model", hjust = 0 , vjust = 0, size = 3.5) +
         annotate("text", x = rep(0.85, 10), y = seq(0.43, 0.15, length = 10),
                  label = c("NB-region", "AR(1)-Climate", "AR(1)", "NB-hex", "RF-season", "NN-season", 
                             "RF-month", "NN-month", "NB-nation", "Ensemble"), size = 3,
                  hjust = 0, vjust = 0) +
         annotate("segment", x = rep(0.82, 10), xend = rep(0.84, 10), y = seq(0.435, 0.155, length = 10), 
                  yend = seq(0.435, 0.155, length = 10), color = mypal[-length(mypal)]),
       width = 6.5, height = 4.5, compression = "lzw")

linescore2 = ggplot(ls_all2) + 
  geom_line(aes(x=year, y=av_ls, color = mod), linewidth = 0.3, alpha = 0.7) + 
  facet_wrap(~clim.regn_plot, nrow=2, labeller = labeller(clim.regn_plot=region_full)) +
  labs(y="mean logarithmic score per year", x="year") +
  scale_color_manual(values = mypal) +
  theme_cowplot(10) + panel_border(color = "black") + theme(legend.position = "none")

ggsave("meanlogscore_line_AA.tiff", 
       plot_grid(linescore2) +
         annotate("text", x = 0.82, y = 0.46, label = "Model", hjust = 0 , vjust = 0, size = 3.5) +
         annotate("text", x = rep(0.85, 11), y = seq(0.43, 0.1, length = 11),
                  label = c("NB-region", "AR(1)-Climate", "AR(1)", "NB-hex", "RF-season", "NN-season", 
                            "RF-month", "NN-month", "NB-nation", "Ensemble", "Always Absent"), size = 3, # 
                  hjust = 0, vjust = 0) +
         annotate("segment", x = rep(0.82, 11), xend = rep(0.84, 11), y = seq(0.435, 0.115, length = 11), 
                  yend = seq(0.435, 0.115, length = 11), color = mypal),
       width = 6.5, height = 4.5, compression = "lzw")
