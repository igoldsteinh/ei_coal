# file for processing eicdf_coal model outputs
# in the case of the alternate skyline prior for R_{u}
# with only five intervals, change times counted from backwards time 0
# to mimic BDMM skyline prior
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
lump_val    = 30.6
n_sim_lumps = 5
# read in posterior samples -----------------------------------------------
posterior_suffix <- paste0("ei_cdf_skyline_sim",
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
# save mcmc summary
write_csv(full_mcmc_summary, here::here("results", "my_mcmc_summaries", 
                                         paste0("mcmc_diagnostics_eicdf_coal_skyline_sim", sim_val, "_local.csv")))
# summarise rt posterior samples ---------------------------------------------
# create grid of true rt values
tree_file_name = sim_dict %>% 
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata <- read_csv(here::here("data", 
                                      "sim_data", 
                                      tree_file_name), show_col_types = FALSE) 
# this is the actual sample time from the simulated data, because coalescence occurs at some time before true time 0 (in forward time)
# ie in reality samples are taken at time 34 in forward time, but the last coalescence occurs at time 6 in forward time
# we can't infer before time 6, so we adjust the time scale so that now time 6 is time 0
# hopefully this means that all subsequent lists remain in the order in which the files are ordered in the list
sim_num_vals <- map(file_list, ~as.numeric(stringr::str_extract(.x, stringr::regex("(\\d+)(?!.*\\d)"))))
samp_time <- map(sim_num_vals, ~fitted_simdata %>% 
                   filter(sim == .x) %>% 
                   summarise(max_time = max(time)) %>% pull())
state_file_name <- sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(state_frame_file)
state_frame <- read_csv(here::here("data",
                                   "sim_data", 
                                   state_file_name), show_col_types = FALSE)
state_frame_list <- map(samp_time, ~state_frame %>%
  filter(time <= .x) %>%
  mutate(reverse_time = .x - time,
         new_time = time,
         lump = pmax(n_sim_lumps - 1L - floor(reverse_time / lump_val), 0L)))

true_rt_list <- map(state_frame_list, ~.x %>% 
  rename("true_rt" = "rt") %>%
  dplyr::select(time, new_time, reverse_time, lump, true_rt) %>%
    mutate(grid_time = ifelse(new_time < max(new_time), ceiling(new_time * 2)/2, new_time),
           diff_time = abs(grid_time - new_time)) %>%
    group_by(grid_time) %>%
    filter(diff_time == min(diff_time)) %>%
    ungroup() %>%
    filter(grid_time <= max(new_time)) %>%
    dplyr::select(lump, time, new_time, reverse_time, grid_time, true_rt)) 
# create timevarying rt quantiles
my_posterior_timevarying_quantiles <- map(full_posterior, ~.x %>% 
  mutate(.iteration = row_number(),
         .chain = 1) %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  dplyr::select( name, value) %>%
  make_timevarying_posterior_quantiles())

my_posterior_rt <- map2(my_posterior_timevarying_quantiles, true_rt_list,  ~.x %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>% 
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(.y, by = c("lump")) %>%
  dplyr::select(lump, time, new_time, reverse_time, value, true_rt, .lower, .upper, .width) %>%
  distinct()) %>%
  map2(file_list,  ~.x %>%
         mutate(address = .y) %>% 
         mutate(sim_num_val = as.numeric(stringr::str_extract(.y, stringr::regex("(\\d+)(?!.*\\d)"))),
                sim = sim_val)) %>%
  bind_rows() 

write_csv(my_posterior_rt, here::here("results", 
                                      "my_generated_quantities", 
                                      paste0("ei_cdf_sim_skyline", sim_val,  "_allseeds_rt_quantiles_local.csv")))
