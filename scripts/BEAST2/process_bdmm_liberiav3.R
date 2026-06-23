# file for converting beast2 .log file
# to csv file for ebola data fit
library(tidyr)
library(dplyr)
library(tibble) 
library(beastio)
library(readr)
# ebola bdmm log file from BEAST2 (see bdmm_liberia_ebolav3.xml)
trace <- beastio::readLog(here::here("scripts", "BEAST2", "bdmm_liberia_ebolav3.log"), burnin = 0.1, as.mcmc = FALSE) # as dataframe for ggplot
df <- as.data.frame(trace)

df <- df %>%
  dplyr::select(Sample, posterior, migrationRateSPEpi.E_to_I, becomeUninfectiousRateSPEpi.I, ends_with("_I") & starts_with("samplingProportionSPEpi") |
           ends_with("I_to_E") & starts_with("ReAmongDemesSPEpi"))

new_names <- c(".iteration", "posterior", "gamma", "nu",
               paste0("samp_prop[", 0:9, "]"),
               paste0("rt_t_values[", 0:11, "]"))
names(df) <- new_names
df <- df %>%
  mutate(.chain = 1)
write_csv(df, here::here("scripts", "BEAST2", "fit_bdmm_liberia_ebolav3.csv"))
