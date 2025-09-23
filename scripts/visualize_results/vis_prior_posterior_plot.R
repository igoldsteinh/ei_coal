# visualize prior posterior
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
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))
source(here::here("src", "utility_functions.R"))
my_theme <- list(
  scale_fill_brewer(name = "Credible Interval Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal_grid(),
  theme())
sim_val = 1
sim_num_val = 1
state_sim_val = 1
true_samp_time = 153
lump_val = 7
sim_dict = read_csv(here::here("data", "sim_data", "sim_dict.csv"), show_col_types = FALSE)
samp_time <- true_samp_time
state_file_name <- sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(state_frame_file)
state_frame <- read_csv(here::here("data",
                                   "sim_data", 
                                   state_file_name), show_col_types = FALSE) %>%
  filter(sim == state_sim_val) %>%
  filter(time <= true_samp_time) %>%
  mutate(reverse_time = abs(time - true_samp_time),
         new_time = abs(reverse_time - samp_time),
         lump = floor(new_time/lump_val)) %>%
  filter(reverse_time <= samp_time)
# data massaging
true_rt <- state_frame %>% 
  rename("true_rt" = "rt") %>%
  dplyr::select(time, new_time, reverse_time, lump, true_rt) 
true_init_rt = true_rt %>% filter(new_time == min(new_time)) %>% pull(true_rt)
true_init_rt = true_rt %>% filter(new_time == min(new_time)) %>% pull(true_rt)
# plot prior posterior ----------------------------------------------------
model_res <- read_csv(here::here("results", 
                                 "my_generated_quantities", 
                                 paste0("ei_cdf_sim", sim_val, "_simnum", sim_num_val, ".csv")))
gamma = 1/4
nu = 1/7

fixed_names <- c("gamma",
                 "nu",
                 "rt_init",
                 "rw_sigma",
                 "e0",
                 "i0"
)

true_vals <- c(gamma,
               nu,
               true_init_rt,
               NA,
               0,
               1
)
true_frame <- data.frame(name = fixed_names, true_value = true_vals)
fixed_posterior <- model_res %>%
  dplyr::select(gamma, nu, `rt_t_values[0]`, rw_sigma, e0, i0) %>%
  rename(rt_init = `rt_t_values[0]`,
         ) %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  left_join(true_frame, by = "name") %>%
  mutate(type = "posterior") %>%
  mutate(name = ifelse(name == "rt_init", "Initial Ru", 
                       ifelse(name == "rw_sigma", "RW SD",
                              ifelse(name == "e0", "E(0)",
                                     ifelse(name == "i0", "I(0)",
                                            ifelse(name == "gamma", "Gamma",
                                                   ifelse(name == "nu", "Nu", name)))))))

# read in fixed prior samps
fixed_prior_samp <- read_csv(here::here("results",
                                        "my_generated_quantities",
                                        paste0("prior_samps_sim",sim_val,".csv"))) %>%
  dplyr::select(gamma, nu, rt_init, rw_sigma, e0, i0) %>%
  pivot_longer(cols = everything(), names_to = "name", values_to = "value") %>%
  left_join(true_frame, by = "name") %>%
  mutate(type = "prior") %>%
  mutate(name = ifelse(name == "rt_init", "Initial Ru", 
                       ifelse(name == "rw_sigma", "RW SD",
                              ifelse(name == "e0", "E(0)",
                                     ifelse(name == "i0", "I(0)",
                                            ifelse(name == "gamma", "Gamma",
                                                   ifelse(name == "nu", "Nu", name)))))))
priors_and_posteriors <- rbind(fixed_posterior %>% dplyr::select(name, value, type, true_value), fixed_prior_samp) %>%
  mutate(type = ifelse(type == "prior", "Prior", "Posterior"))
param_plot <- priors_and_posteriors %>%
  rename(Type = type) %>%
  ggplot(aes(value, Type, fill = Type)) +
  stat_halfeye(normalize = "xy")  +
  geom_vline(aes(xintercept = true_value), linetype = "dotted", size = 1) + 
  facet_wrap(. ~ name, scales = "free_x") +
  theme_bw() +
  theme(text = element_text(size = 20),
        legend.position = c(0.25,0.15),
        legend.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()) +
  ggtitle("Fixed Parameter Prior/Posterior Distributions (Fixed Scenario)") + 
  xlab("Value") + 
  ylab("Type")
ggsave(here::here( "figures", "prior_posterior_plot.pdf"),
       param_plot,
       width = 15, height = 7)
