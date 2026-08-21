# file for processing eicdf_coal model outputs
# en masse
# for the purposes of checking comptime runs to see if they converged enough
# so we can calculate average ess/ sec
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
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  sim_val = 9
} else {
  sim_val <- as.integer(args[1])
}
true_samp_time = 153
lump_val = 7
# read in posterior samples -----------------------------------------------
posterior_suffix <- paste0("ei_cdf_sim",
                             sim_val, "_simnum")
file_list <- list.files(path = path("results", "my_generated_quantities"),  pattern = posterior_suffix)
full_posterior <- map(file_list, ~read_csv(here::here("results", "my_generated_quantities", .x)))


# summarise mcmc diagnostics ----------------------------------------------
full_posterior_converted <- map(full_posterior, ~.x %>% 
                                  mutate(.iteration = row_number(),
         .chain = 1) %>% 
     dplyr::select(-log_likelihood, -actual_iteration) %>%
  as_draws())
full_mcmc_summary <- map(full_posterior_converted, ~summarise_draws(.x, "ess_basic","ess_bulk", "ess_tail") %>% 
  dplyr::select(variable, ess_basic, ess_bulk, ess_tail))
full_mcmc_summary <- map2(file_list, full_mcmc_summary, ~.y %>%
                            mutate(address = .x) %>% 
                            mutate(sim_num_val = as.numeric(stringr::str_extract(.x, stringr::regex("(\\d+)(?!.*\\d)")))))
full_mcmc_summary <- bind_rows(full_mcmc_summary)


fails <- full_mcmc_summary %>%
  group_by(sim_num_val) %>%
  filter(sim_num_val <= 20) %>%
  summarise(min_ess = min(ess_basic),
            min_bulk = min(ess_bulk),
            min_tail = min(ess_tail)) %>%
  filter(min_ess < 100 | min_bulk < 100 | min_tail < 100)
