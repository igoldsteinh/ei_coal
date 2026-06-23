# visualize results from lbr
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(kableExtra))
library(ape)
library(ggtree)
library(phylodyn)
library(treeio)

source(here::here("src", "utility_functions.R"))

my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal(),
  theme())


# load data and model res -------------------------------------------------
model_res <- read_csv(here::here("scripts",
                                 "BEAST2",
                                 "fit_bdmm_liberia_ebolav3.csv"))
lump_val     = 30
n_bdmm_lumps = 12  # 11 change points → 12 intervals
date_crosswalk <- read_csv(here::here("data", "real_data", "ebola_lbrmcc_2014_date_crosswalkv2.csv")) %>%
  mutate(lump = pmax(n_bdmm_lumps - 1L - floor(reverse_times / lump_val), 0L))




# model diags -------------------------------------------------------------


DA_samples <-model_res %>%
  as_draws()
DA_mcmc_summary <- summarise_draws(DA_samples, "ess_basic","ess_bulk", "ess_tail") %>%
  dplyr::select(variable, ess_basic, ess_bulk, ess_tail) %>% filter(variable != "actual_iteration")
DA_mcmc_summary

log_lik_plot <- DA_samples %>%
  ggplot(aes(x = .iteration, y = posterior)) +
  geom_line() +
  theme_bw() +
  ggtitle(paste0("Log Likelihood LBR "))
log_lik_plot


# rt ----------------------------------------------------------------------
bdmm_my_posterior_timevarying_quantiles <- model_res %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  dplyr::select( name, value) %>%
  make_timevarying_posterior_quantiles()

bdmm_my_posterior_rt <- bdmm_my_posterior_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>%
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(date_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, forward_times, reverse_dates, reverse_times, value, .lower, .upper, .width) %>%
  distinct() %>%
  arrange(reverse_dates)


bdmm_plot_posterior_rt <- bdmm_my_posterior_rt %>%
  ggplot(aes(x = reverse_dates, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("BDMM") +
  my_theme +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank()) +
  ylab("Ru") +
  xlab("Time")

bdmm_plot_posterior_rt


# my model ----------------------------------------------------------------

# load data and model res -------------------------------------------------
model_res <- read_csv(here::here("results",
                                 "my_generated_quantities",
                                 "ess_joint_hetchron_inseq_samples_simLBR_version4_seed123467.csv"))
lump_val = 7
date_crosswalk <- read_csv(here::here("data", "real_data", "ebola_lbrmcc_2014_date_crosswalkv2.csv")) %>%
  mutate(lump = floor(forward_times/lump_val))




# model diags -------------------------------------------------------------


DA_samples <-model_res %>%
  mutate(.iteration = row_number(),
         .chain = 1) %>%
  filter(.iteration > 500) %>%
  as_draws()
DA_mcmc_summary <- summarise_draws(DA_samples, "ess_basic","ess_bulk", "ess_tail") %>%
  dplyr::select(variable, ess_basic, ess_bulk, ess_tail) %>% filter(variable != "actual_iteration")
DA_mcmc_summary

log_lik_plot <- DA_samples %>%
  ggplot(aes(x = .iteration, y = log_likelihood)) +
  geom_line() +
  theme_bw() +
  ggtitle(paste0("Log Likelihood LBR "))
log_lik_plot


# rt ----------------------------------------------------------------------
my_posterior_timevarying_quantiles <- model_res %>%
  mutate(.iteration = row_number(),
         .chain = 1) %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  dplyr::select( name, value) %>%
  make_timevarying_posterior_quantiles()

my_posterior_rt <- my_posterior_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>%
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(date_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, forward_times, reverse_dates, reverse_times, value, .lower, .upper, .width) %>%
  distinct() %>%
  arrange(reverse_dates)


my_plot_posterior_rt <- my_posterior_rt %>%
  ggplot(aes(x = reverse_dates, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EI Coal") +
  my_theme +
  theme(text = element_text(size = 20),
        legend.position = c(0.78,0.88),
        legend.background = element_blank()) +
  ylim(c(0, 2)) +
  ylab("Ru") +
  xlab("Time")

my_plot_posterior_rt


# phydyn ------------------------------------------------------------------
# note this file was not saved to the repo
# as it is quite large
# load data and model res -------------------------------------------------
phydyn_res <- read_csv(here::here("scripts",
                                 "phydynR",
                                 "phydynR_ebola_res.csv"))
lump_val = 7
date_crosswalk <- read_csv(here::here("data", "real_data", "ebola_lbrmcc_2014_date_crosswalkv2.csv")) %>%
  mutate(lump = floor(forward_times/lump_val))




# model diags -------------------------------------------------------------

phydyn_res <- phydyn_res %>%
  group_by(.chain) %>%
  filter(.chain != 4) %>%
  mutate(.iteration = row_number()) %>%
  filter(.iteration > 25000)
DA_samples <- phydyn_res %>%
  as_draws()
DA_mcmc_summary <- summarise_draws(DA_samples, "ess_basic","ess_bulk", "ess_tail", "rhat") %>%
  dplyr::select(variable, ess_basic, ess_bulk, ess_tail, rhat) %>% filter(variable != "actual_iteration")
DA_mcmc_summary

log_lik_plot <- DA_samples %>%
  ggplot(aes(x = .iteration, y = nu, color = as.factor(.chain))) +
  geom_line() +
  theme_bw() +
  ggtitle(paste0("Log Likelihood LBR "))
log_lik_plot


# rt ----------------------------------------------------------------------
my_posterior_timevarying_quantiles <- phydyn_res %>%
  mutate(.iteration = row_number(),
         .chain = 1) %>%
  pivot_longer(-c(.iteration, .chain)) %>%
  dplyr::select( name, value) %>%
  make_timevarying_posterior_quantiles()

my_posterior_rt <- my_posterior_timevarying_quantiles %>%
  filter(name == "rt_t_values") %>%
  rename(lump = time) %>%
  dplyr::select(lump, value, .lower, .upper, .width,.point, .interval) %>%
  right_join(date_crosswalk, by = c("lump")) %>%
  dplyr::select(lump, forward_times, reverse_dates, reverse_times, value, .lower, .upper, .width) %>%
  distinct() %>%
  arrange(reverse_dates)


phydyn_plot_posterior_rt <- my_posterior_rt %>%
  ggplot(aes(x = reverse_dates, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Phydyn") +
  my_theme +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank()) +
  ylim(c(0, 2)) +
  ylab("Ru") +
  xlab("Time")

phydyn_plot_posterior_rt



# visualize mcc tree ------------------------------------------------------

mcc_tree <- read.nexus(here::here("scripts", "ebola_scripts", "liberia_mcctreev2"))


lbr_list <- summarize_phylo(mcc_tree)
fields <- strsplit(mcc_tree$tip.label, "\\|")
sample_dates <- sapply(fields, function(x) x[6])
max_date = max(sample_dates)
tree_plot<- ggtree(mcc_tree, mrsd = max_date, as.Date = TRUE) +
  ggtitle("Liberia (2014 - 2015) MCC Tree") +
  theme_tree2() +
  xlab("Time") +
  theme(text = element_text(size = 20))

# combined ----------------------------------------------------------------
combined_plot <- tree_plot + my_plot_posterior_rt + bdmm_plot_posterior_rt + phydyn_plot_posterior_rt
ggsave(here::here("figures", "liberia_plot_allmodelsv3.pdf"),
       combined_plot,
       width = 14, height = 6, units = "in")
