# file for converting beast log files fit to simulated data
# into summaries and plots showing model performance
# log files were generated using fit_bdmm_control.xml as the template file
# outputs are mcmc diagnostic file, rt quantile file, and plots of rt credible intervals
library(tidyr)
library(tibble) 
library(beastio)
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
lump_val       = 30.6  # interval width in days (153 / 5 intervals, matching 4 change points)
n_sim_lumps    = 5     # number of Rt intervals

# read in posterior samples -----------------------------------------------
posterior_suffix <- "*.log"
file_list <- list.files(path = path("BEAST2", "control_xml2"), pattern = posterior_suffix)

full_posterior <- map(file_list, ~beastio::readLog(here::here("BEAST2", "control_xml2", .x), burnin = 0.1, as.mcmc = FALSE))
new_names <- c(".iteration", "posterior", "nu", "gamma",
               "rt_t_values[0]", "rt_t_values[1]", "rt_t_values[2]",
               "rt_t_values[3]", "rt_t_values[4]")

full_posterior_df <- map(full_posterior, ~as.data.frame(.x)) %>%
  map(~.x %>%
        dplyr::select(1, 2, 8, 9, 14, 17, 20, 23, 25)) %>%
  map(~setNames(.x, new_names)) %>%
  map(~.x %>%
        mutate(.chain = 1))
# summarise mcmc diagnostics ----------------------------------------------
full_posterior_converted <- map(full_posterior_df, ~.x %>%
                                  dplyr::select(-posterior) %>%
                                  as_draws())
full_mcmc_summary <- map(full_posterior_converted, ~summarise_draws(.x, "ess_basic","ess_bulk", "ess_tail") %>%
                           dplyr::select(variable, ess_basic, ess_bulk, ess_tail))
full_mcmc_summary <- map2(file_list, full_mcmc_summary, ~.y %>%
                            mutate(address = .x) %>%
                            mutate(sim_num_val = as.numeric(stringr::str_extract(.x, stringr::regex("(\\d+)(?!.*\\d)")))))
full_mcmc_summary <- bind_rows(full_mcmc_summary)
# save mcmc summary
write_csv(full_mcmc_summary, here::here("BEAST2",
                                        paste0("mcmc_diagnostics_bdmm_sim", sim_val, "v2.csv")))
# summarise rt posterior samples ---------------------------------------------
# create grid of true rt values
tree_file_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata <- read_csv(here::here("data",
                                      "sim_data",
                                      tree_file_name), show_col_types = FALSE)
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
my_posterior_timevarying_quantiles <- map(full_posterior_df, ~.x %>%
                                            pivot_longer(-c(.iteration, .chain)) %>%
                                            dplyr::select(name, value) %>%
                                            make_timevarying_posterior_quantiles())

my_posterior_rt <- map2(my_posterior_timevarying_quantiles, true_rt_list,  ~.x %>%
                          filter(name == "rt_t_values") %>%
                          rename(lump = time) %>%
                          dplyr::select(lump, value, .lower, .upper, .width, .point, .interval) %>%
                          right_join(.y, by = c("lump")) %>%
                          dplyr::select(lump, time, new_time, reverse_time, value, true_rt, .lower, .upper, .width) %>%
                          distinct()) %>%
  map2(file_list,  ~.x %>%
         mutate(address = .y) %>%
         mutate(sim_num_val = as.numeric(stringr::str_extract(.y, stringr::regex("(\\d+)(?!.*\\d)"))),
                sim = sim_val)) %>%
  bind_rows()

write_csv(my_posterior_rt, here::here("BEAST2",
                                      paste0("bdmm_sim", sim_val, "_allseeds_rt_quantilesv2.csv")))
# visualize results
# all credit to Damon Bayer for plot functions
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme(legend.position = "bottom"))

make_rt_plot <- function(seed_val) {
  my_posterior_rt %>%
    filter(sim_num_val == seed_val) %>%
    ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
    geom_lineribbon() +
    geom_point(aes(time, true_rt), color = "coral1") +
    scale_y_continuous("Rt", label = comma) +
    scale_x_continuous(name = "Time") +
    ggtitle(str_c("EI Coal Posterior Rt Scenario ", sim_val, " Simnum ", seed_val)) +
    my_theme
}
ggsave2(filename = here::here("BEAST2", paste0("bdmm_rtplots_sim", sim_val, "v2.pdf")),
        plot = my_posterior_rt %>%
          distinct(sim_num_val) %>%
          arrange(sim_num_val) %>%
          pull(sim_num_val) %>%
          map(make_rt_plot) %>%
          marrangeGrob(ncol = 1, nrow = 1),
        width = 12,
        height = 8)
