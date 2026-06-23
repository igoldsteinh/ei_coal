# summarise simulation results
# from different models EI Coal, BDMM and PhydynR
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(stringr)
library(tidybayes)
library(posterior)
library(fs)
library(gridExtra)
library(ggplot2)
library(scales)
library(cowplot)
source(here::here("src", "utility_functions.R"))
sim_dict <- read_csv(here::here("data", "sim_data", "sim_dict.csv"), 
                     show_col_types = FALSE)

# mcmc diags --------------------------------------------------------------
mcmc_sim_bdmm <- read_csv(here::here("scripts",
                                     "BEAST2",
                                 "mcmc_diagnostics_bdmm_sim9v2.csv"))

fails_bdmm <- mcmc_sim_bdmm %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)

mcmc_sim_phydynR <- read_csv(here::here("scripts",
                                     "phydynR",
                                     "mcmc_diagnostics_phydynR_sim9.csv"))

fails_phydynR <- mcmc_sim_phydynR %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail),
            max_rhat = max(rhat)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100 | max_rhat > 1.1)

# read in the rt quantiles ------------------------------------------------
width = 0.95
rt9 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim9_allseeds_rt_quantiles.csv")) %>%
  filter(.width == width)

bdmm_rt9 <- read_csv(here::here("scripts", 
                                "BEAST2", 
                                "bdmm_sim9_allseeds_rt_quantilesv2.csv")) %>%
  filter(.width == width) %>%
  filter(address != "beast_controliso50_simnum1_oldversion.log")

phydyn_rt9 <- read_csv(here::here("scripts", 
                                "phydynR", 
                                "phdynR_sim9_allseeds_rt_quantiles.csv")) %>%
  filter(.width == width)

# calculate metrics -------------------------------------------------------
metrics9 = NULL
sim_numvals9 = unique(rt9$sim_num_val)
for (i in 1:length(sim_numvals9)) {
  sim_num = sim_numvals9[i]
  sub_frame = rt9 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf", 
           sim = 9)
  metrics9 = bind_rows(metrics9, metrics)
}

bdmm_metrics9 = NULL
bdmm_sim_numvals9 = unique(bdmm_rt9$sim_num_val)
for (i in 1:length(bdmm_sim_numvals9)) {
  sim_num = bdmm_sim_numvals9[i]
  sub_frame = bdmm_rt9 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "bdmm", 
           sim = 9)
  bdmm_metrics9 = bind_rows(bdmm_metrics9, metrics)
}

phydynR_metrics9 = NULL
phydynR_sim_numvals9 = unique(phydyn_rt9$sim_num_val)
for (i in 1:length(phydynR_sim_numvals9)) {
  sim_num = phydynR_sim_numvals9[i]
  sub_frame = phydyn_rt9 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "phydynR", 
           sim = 9)
  phydynR_metrics9 = bind_rows(phydynR_metrics9, metrics)
}


all_metrics <- bind_rows(metrics9, bdmm_metrics9, phydynR_metrics9) %>%
  group_by(model) %>%
  summarise(
    median_env = median(mean_env),
    env_lower = quantile(mean_env, 0.025),
    env_higher = quantile(mean_env, 0.975),
    median_dev = median(mean_dev),
    dev_lower = quantile(mean_dev, 0.025),
    dev_higher = quantile(mean_dev, 0.975),
    medianMCIW = median(MCIW),
    MCIW_lower = quantile(MCIW, 0.025),
    MCIW_upper = quantile(MCIW, 0.975),
    num_sims = n()
  )

library(xtable)
xtable(all_metrics)