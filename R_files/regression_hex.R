###
# Author: Karen Holcomb
# Last Updated: June 8, 2023
##
# Fit GLM models with Stan to compare model performance
#  using logarithmic score by hexagon and year per model, WNNV history
##

library(rstanarm)
library(tidyr)
library(dplyr)
library(emmeans)

ls_all <- readRDS("../data/logscore_allmods.rds") #hex-year data for all models

ens <- ls_all %>% mutate(loc_prob = exp(mean_ls)) %>% #back to probability scale
  subset(mod != "score_aa") %>% #drop Always Absent model from calculation of median ensemble
  group_by(location, clim.regn, year) %>% summarise(ens_prob = median(loc_prob), #median probability for ensemble
                                                    ens_ls = log(ens_prob)) %>% #logscore for ensemble
  ungroup() %>% rename(mean_ls = ens_ls) %>% mutate(mod = "ens") %>% select(-ens_prob)

ls_all <- rbind(ls_all, ens) %>%
  mutate(sup = -1 * pmax(-10, mean_ls)) %>% #no -Inf for complete misses
  mutate(sup = pmax(1e-10, sup), #no 0s for Gamma distribution (min sup > 0 -> 1.522998e-09)
         location = factor(location), year = factor(year)) #categorical variables


# Few metrics of historic cases for regression modeling
df_all <- readRDS("../data/WNND_hex_2005-2021.rds") #WNNV cases per year per hex

da_med <- do.call(rbind.data.frame, lapply(unique(df_all$fips), function(x) {
  do.call(rbind.data.frame, lapply(2005:2021, function(y) {
    data.frame(fips = x, year = y, 
               med_hist = median(da$tot_count[da$fips == x & da$year < y]), #historic median (2005-previous year)
               hist_sum = sum(da$tot_count[da$fips == x & da$year < y])) #number of cases (2005-previous year)
  }))
}))

da <- left_join(df_all, da_med) %>% 
  mutate(fips = as.character(fips), year = as.character(year))

ls_all <- left_join(ls_all, da, by = c("location" = "fips", "year")) #all into a single datafrmae

saveRDS(ls_all, "../data/sup_cases_regdat.rds") #save dataframe of all data for use in regression models below


## Fit regression models (using Stan on surprisal scale) for estimating differences in prediction
ls_all <- readRDS("../data/sup_cases_regdat.rds") %>%
  mutate(lcases = log(tot_count + 1),
         hist_cases = hist_sum > 0,
         obs_cases = tot_count > 0)

# Model performance overall
gamma_log_mod <- stan_glm(sup ~ mod, data = ls_all, family = Gamma(link = "log"), chains = 4, iter = 1000, cores = 4)

# pairwise comparisons of model performance
comps <- as.data.frame(contrast(emmeans(gamma_log_mod, ~mod), method = "pairwise"))
NBh_comps <- comps[grep("score_nb$|score_nb -", comps$contrast),] #contrasts relative to NB_hex ("baseline" model)
NBh_comps[grep("score_nb -", NBh_comps$contrast), 2:4] <- -1 * NBh_comps[grep("score_nb -", NBh_comps$contrast), 2:4] 
    #reverse direction of comps so all estimates of models relative to NB_hex
NBh_comps <- NBh_comps %>% arrange(estimate) %>%
  mutate(mod = sapply(stringr::str_split(contrast, " - "), function(x) x[x != "score_nb"]), #which model is compared
         nord = 1:n()) #order by diff in estimates
  
emmplot <- ggplot(NBh_comps[NBh_comps$mod != "score_aa",]) + #drop AA model since so much larger in magnitude of diff than others
  geom_hline(yintercept = 0, col = "darkgrey") +
  geom_point(aes(x=nord, y=estimate), size = 1) + 
  geom_errorbar(aes(x=nord, ymin=lower.HPD, ymax=upper.HPD), width = 0.5) +
  scale_x_continuous("model", breaks = 1:9, 
                     labels = c("NB-region", "Ensemble", "AR(1)-Climate", "AR(1)", "RF-season", "NN-season", "RF-month", 
                                "NN-month", "NB-nation")) +
  labs(y = "median difference in estimated\nmarginal mean vs. NB-hex (95% HPD)") +
  theme_cowplot(8) + panel_border(color="black", size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

write.csv(file = "../HPD_diff.csv", comps %>% mutate(sig = sign(lower.HPD) == sign(upper.HPD)) %>% 
            subset(sig) %>% mutate_at(2:4, round, 3), row.names = F) #csv of comparisons with HPD != 0
#rounded to 3 sigfigs, cleaned names of models and columns

# Scatter plot of observed cases vs. logscore (NB-hex as example since baseline model) -> similar performance across all models
sup_case <- ggplot(ls_all[ls_all$mod == "score_nb",]) +
  geom_hline(yintercept = 0, col = "darkgrey") +
  geom_point(aes(x=log10(exp(lcases)), y=log10((-1*mean_ls) + 1)), shape = 16, alpha = 0.2) + #convert to log10 cases
  scale_x_continuous("observed WNND cases", breaks = log10(c(0, 1, 10, 100, 1000) +1), labels = c(0, 1, 10, 100, 1000),
                     expand = expansion(0,0.06)) +
  scale_y_continuous("surprisal (per hexagon and year)", breaks = log10(c(0, 1, 5, 25) +1), labels = c(0, 1, 5, 25)) +
  theme_minimal_grid(8)

toplots2 <- plot_grid(emmplot, sup_case, nrow = 1, rel_widths = c(0.9, 1), labels = "AUTO")


# Model performance based on model, historical cases, and observed cases
# Note: ran this stan model on a cluster, took ~2hrs
gamma_log2 <- stan_glm(sup ~ lcases + mod*hist_cases*obs_cases, data = ls_all,
                       family = Gamma(link = "log"), cores = 2, chains = 4, iter = 1000)

nd = expand.grid(mod = unique(gamma_log2$data$mod), lcases = log(c(0, 1, 5)+1), hist_cases = c(TRUE, FALSE))
nd2 = nd %>%  mutate(obs_cases = lcases > 0)

preds2 <- posterior_epred(gamma_log2, newdata = nd2, draws = 500) #500 samples from posterior for predictions
q_pp2 <- do.call(rbind.data.frame, lapply(1:ncol(preds2), 
                                         function(x) c(quantile(preds2[,x], c(0.025, 0.5, 0.975)),#median and 95% PI for all nd predictions
                                                       mean(preds2[,x]), sd(preds2[,x])/sqrt(length(preds2[,x]))))) #mean and se for CIs
names(q_pp2) <- c("l95", "med", "u95", "mean", "se")

preds_df2 <- cbind(nd2, q_pp2)
# preds_df2 %>% arrange(med) #order of models, re-arrange levels below based on this

preds_df2$mod <- factor(preds_df2$mod, levels = c("score_nb_r", "ens", "ls_ar1Clim", "ls_ar1", "score_nb", "mscoreRF_nbs", 
                                                  "ls_nn.s_nb", "mscoreRF_nbm", "ls_nn.m_nb", "score_nb_u", "score_aa")) #order to emmeans

hs_namers = c("FALSE" = "no historical cases reported", "TRUE" = "historical case(s) reported")

q2 <- ggplot(preds_df2[preds_df2$mod != "score_aa",]) + #predicted surprisal w/ SE per model across year and location predicted -> from stan models
  geom_point(aes(x=mod, y=med, color=factor(lcases)), size = 1) +
  geom_errorbar(aes(x=mod, ymin = l95, ymax=u95, col=factor(lcases)), width = 0.5) + 
  facet_wrap(~hist_cases, labeller = labeller(hist_cases = hs_namers), scale = "free_y") +
  scale_x_discrete(labels = c("NB-region", "Ensemble", "AR(1)-Climate", "AR(1)", "NB-hex", "RF-season", 
                              "NN-season", "RF-month", "NN-month", "NB-nation")) + 
  labs(x = "model", y = "median estimated surprisal (95% PI)", color = "Obs cases", title = "Q2") +
  scale_color_discrete("observed\ncases", labels = c(0,1,5)) +
  theme_cowplot(8) + panel_border(color="black", size = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), legend.position = "right")

ggsave("mod_comps.tiff" ,
       plot_grid(toplots2, q2 + labs(title = NULL), nrow = 2, labels = c("","C")),
       width = 6.5, height = 6.5, units = "in", compression = "lzw")


# Model performance by region
gamma_log_reg <- stan_glm(sup ~ mod * clim.regn, data = ls_all, family = Gamma(link = "log"), chains = 4, 
                          iter = 1000, cores = 4)
cr_pairs <- pairs(emmeans(gamma_log_reg, ~mod*clim.regn), by = "clim.regn")

df_regn <- do.call(rbind.data.frame,
                   list(as.data.frame(subset(cr_pairs, clim.regn == "West")), as.data.frame(subset(cr_pairs, clim.regn == "Northwest")), 
                        as.data.frame(subset(cr_pairs, clim.regn == "NRockiesandPlains")), 
                        as.data.frame(subset(cr_pairs, clim.regn == "South")), as.data.frame(subset(cr_pairs, clim.regn == "UpperMidwest")),
                        as.data.frame(subset(cr_pairs, clim.regn == "Southwest")), as.data.frame(subset(cr_pairs, clim.regn == "OhioValley")), 
                        as.data.frame(subset(cr_pairs, clim.regn == "Southeast")), as.data.frame(subset(cr_pairs, clim.regn == "Northeast"))))

NBh_regn <- df_regn[grep("score_nb$|score_nb -", df_regn$contrast),] #contrasts relative to NB_hex ("baseline" model)
NBh_regn <- NBh_regn %>% group_by(clim.regn) %>% 
  mutate(mod = sapply(stringr::str_split(contrast, " - "), function(x) x[x != "score_nb"]), #which model is compared
         comp_dir = sapply(stringr::str_split(contrast, " - "), function(x) which(x == "score_nb")), #order of contrast
         est_ord = ifelse(comp_dir == 1, -1 * estimate, estimate), #reverse direction so all "mod - score_nb" contrast
         lHPD_ord = ifelse(comp_dir == 1, -1 * lower.HPD, lower.HPD),
         uHPD_ord = ifelse(comp_dir == 1, -1 * upper.HPD, upper.HPD)) %>% arrange(est_ord) %>%
  mutate(nord = 1:n()) #order by diff in estimates

NBh_regn$clim.regn_plot <- as.numeric(factor(NBh_regn$clim.regn, levels = c("Northwest", "NRockiesandPlains", "UpperMidwest", "OhioValley", "Northeast",
                                                                            "West", "Southwest","South", "Southeast")))
region_full = c("1"="Northwest", "2"="N Rockies & Plains", "3"="Upper Midwest", "4"="Ohio Valley", "5"="Northeast", "6"="West", 
                "7"="Southwest", "8"="South", "9"="Southeast")

mypal2 <- c("#e6194B","#800000","#f58231","#3cb44b","#42d4f4","#469990","#000075","#911eb4","darkgrey") 

NBh_regn$mod <- factor(NBh_regn$mod, levels = c("score_nb_r", "ls_ar1Clim", "ls_ar1", "score_nb", "mscoreRF_nbs", 
                                                "ls_nn.s_nb", "mscoreRF_nbm" , "ls_nn.m_nb", "score_nb_u", "ens", "score_aa"),
                       labels = c("NB-region", "AR(1)-Climate", "AR(1)", "NB-hex", "RF-season", 
                                  "NN-season", "RF-month" , "NN-month", "NB-nation", "Ensemble", "Always Absent"))

reg_emm <- ggplot(NBh_regn[NBh_regn$mod != "Always Absent",]) + #drop AA model since so much larger in magnitude of diff than others
  geom_hline(yintercept = 0, col = "#808000") + 
  geom_point(aes(x=nord, y=est_ord, col=mod), size = 1) + 
  geom_errorbar(aes(x=nord, ymin=lHPD_ord, ymax=uHPD_ord, col=mod), width = 0.5, linewidth = 0.5) +
  facet_wrap(~clim.regn_plot, nrow=2, labeller = labeller(clim.regn_plot=region_full)) + 
  labs(y = "median difference in estimated \nmarginal mean vs. NB-hex (95% HPD)", x = "") +
  scale_color_manual("Model", values = mypal2) +
  theme_cowplot(8) + panel_border(color="black", size = 0.5) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = c(0.81,0.22))

ggsave("region_emmeans.tiff", reg_emm,
       width = 6.5, height = 6.4/2, units = "in", compression = "lzw")

write.csv(file = "../HPD_diff_regn.csv", df_regn %>% mutate(sig = sign(lower.HPD) == sign(upper.HPD)) %>% 
            subset(sig) %>% mutate_at(3:5, round, 3), row.names = F)
#rounded to 3 sigfigs, cleaned names of models and columns
