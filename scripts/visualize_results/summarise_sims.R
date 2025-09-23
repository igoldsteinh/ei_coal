# summarise simulation results
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
mcmc_sim1 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim1.csv"))
mcmc_sim2 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim2.csv"))

mcmc_sim3 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim3.csv"))
mcmc_sim4 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim4.csv"))

mcmc_sim5 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim5.csv"))
mcmc_sim6 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim6.csv")) %>%
  filter(!(sim_num_val %in% c(22,28,32)))
mcmc_sim6_local <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim6_local.csv"))
mcmc_sim6 <- bind_rows(mcmc_sim6, mcmc_sim6_local)
mcmc_sim7 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim7.csv"))
mcmc_sim8 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim8.csv")) %>%
  filter(sim_num_val != 14)
mcmc_sim8_local <- read_csv(here::here("results", 
                                       "my_mcmc_summaries",
                                       "mcmc_diagnostics_eicdf_coal_sim8_local.csv"))
mcmc_sim8 <- bind_rows(mcmc_sim8, mcmc_sim8_local)

mcmc_sim9 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim9.csv"))
mcmc_sim10 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim10.csv"))

mcmc_sim11 <- read_csv(here::here("results",
                                 "my_mcmc_summaries",
                                 "mcmc_diagnostics_eicdf_coal_sim11.csv"))

mcmc_sim12 <- read_csv(here::here("results",
                                  "my_mcmc_summaries",
                                  "mcmc_diagnostics_eicdf_coal_sim12.csv"))


fails1 <- mcmc_sim1 %>% 
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100) 
fails2 <- mcmc_sim2 %>% 
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100) 

fails3 <- mcmc_sim3 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)
fails4 <- mcmc_sim4 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)

fails5 <- mcmc_sim5 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)
fails6 <- mcmc_sim6 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)

fails7 <- mcmc_sim7 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)
fails8 <- mcmc_sim8 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)

fails9 <- mcmc_sim9 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)
fails10 <- mcmc_sim10 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)

fails11 <- mcmc_sim11 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)

fails12 <- mcmc_sim12 %>%
  group_by(sim_num_val) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)
# read in the rt quantiles ------------------------------------------------
width = 0.95
rt1 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim1_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)
width = 0.95
rt2 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim2_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)

rt3 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim3_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)
rt4 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim4_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)

rt5 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim5_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)
rt6 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim6_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) %>%
  filter(!(sim_num_val %in% c(22,28,32)))
rt6_local <- read_csv(here::here("results",
                                 "my_generated_quantities",
                                 "ei_cdf_sim6_allseeds_rt_quantiles_local.csv"))
rt6 <- bind_rows(rt6, rt6_local) %>% filter(.width == width)

rt7 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim7_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)
rt8 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim8_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width) 
debug <- rt8 %>% filter(sim_num_val == 44)
rt9 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim9_allseeds_rt_quantiles.csv")) %>%
  filter(.width == width)

debugging <- rt9 %>% filter(.width == 0.95) %>% filter(sim_num_val == 2)
rt10 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim10_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)

rt11 <- read_csv(here::here("results", 
                           "my_generated_quantities", 
                           "ei_cdf_sim11_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)

rt12 <- read_csv(here::here("results", 
                            "my_generated_quantities", 
                            "ei_cdf_sim12_allseeds_rt_quantiles.csv")) %>% 
  filter(.width == width)
# calculate metrics -------------------------------------------------------
metrics1 = NULL
sim_numvals1 = unique(rt1$sim_num_val)
for (i in 1:length(sim_numvals1)) {
  sim_num = sim_numvals1[i]
  sub_frame = rt1 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 1)
  metrics1 = bind_rows(metrics1, metrics)
}
metrics2 = NULL
sim_numvals2 = unique(rt2$sim_num_val)
for (i in 1:length(sim_numvals2)) {
  sim_num = sim_numvals2[i]
  sub_frame = rt2 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf", 
           sim = 2)
  metrics2 = bind_rows(metrics2, metrics)
}
metrics3 = NULL
sim_numvals3 = unique(rt3$sim_num_val)
for (i in 1:length(sim_numvals3)) {
  sim_num = sim_numvals3[i]
  sub_frame = rt3 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 3)
  metrics3 = bind_rows(metrics3, metrics)
}
metrics4 = NULL
sim_numvals4 = unique(rt4$sim_num_val)
for (i in 1:length(sim_numvals4)) {
  sim_num = sim_numvals4[i]
  sub_frame = rt4 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 4)
  metrics4 = bind_rows(metrics4, metrics)
}
metrics5 = NULL
sim_numvals5 = unique(rt5$sim_num_val)
for (i in 1:length(sim_numvals5)) {
  sim_num = sim_numvals5[i]
  sub_frame = rt5 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 5)
  metrics5 = bind_rows(metrics5, metrics)
}
metrics6 = NULL
sim_numvals6 = unique(rt6$sim_num_val)
for (i in 1:length(sim_numvals6)) {
  sim_num = sim_numvals6[i]
  sub_frame = rt6 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 6)
  metrics6 = bind_rows(metrics6, metrics)
}
metrics7 = NULL
sim_numvals7 = unique(rt7$sim_num_val)
for (i in 1:length(sim_numvals7)) {
  sim_num = sim_numvals7[i]
  sub_frame = rt7 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf", 
           sim = 7)
  metrics7 = bind_rows(metrics7, metrics)
}
metrics8 = NULL
sim_numvals8 = unique(rt8$sim_num_val)
for (i in 1:length(sim_numvals8)) {
  sim_num = sim_numvals8[i]
  sub_frame = rt8 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf", 
           sim = 8)
  metrics8 = bind_rows(metrics8, metrics)
}
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
metrics10 = NULL
sim_numvals10 = unique(rt10$sim_num_val)
for (i in 1:length(sim_numvals10)) {
  sim_num = sim_numvals10[i]
  sub_frame = rt10 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf", 
           sim = 10)
  metrics10 = bind_rows(metrics10, metrics)
}
metrics11 = NULL
sim_numvals11 = unique(rt11$sim_num_val)
for (i in 1:length(sim_numvals11)) {
  sim_num = sim_numvals11[i]
  sub_frame = rt11 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 11)
  metrics11 = bind_rows(metrics11, metrics)
}
metrics12 = NULL
sim_numvals12 = unique(rt12$sim_num_val)
for (i in 1:length(sim_numvals12)) {
  sim_num = sim_numvals12[i]
  sub_frame = rt12 %>% filter(sim_num_val == sim_num) 
  metrics = rt_metrics(sub_frame, value = value, upper = .upper, lower = .lower) %>% 
    mutate(seed = sim_num,
           model = "EIcdf",
           sim = 12)
  metrics12 = bind_rows(metrics12, metrics)
}
all_metrics <- bind_rows(metrics1, metrics2, metrics3,
                         metrics4,metrics5, metrics6, metrics7, metrics8,
                         metrics9, metrics10, metrics11, metrics12) %>%
  group_by(sim) %>%
  summarise(
    meanmean_env = mean(mean_env),
    sd_env = sd(mean_env),
    meanmean_dev = mean(mean_dev),
    sd_dev = sd(mean_dev),
    meanMCIW = mean(MCIW),
    sdMCIW = sd(MCIW),
    mean_diffMASV = mean(abs(MASV - true_MASV)),
    sd_diffMASV = sd(abs(MASV - true_MASV)),
    num_sims = n()
  )
all_metrics2 <- bind_rows(metrics1, metrics2, metrics3,
                         metrics4,
                         metrics5, metrics6, metrics7, metrics8,
                         metrics9, metrics10, metrics11, metrics12) %>%
  group_by(sim) %>%
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
xtable(all_metrics2)