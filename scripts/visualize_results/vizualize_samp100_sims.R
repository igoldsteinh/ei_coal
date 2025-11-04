# visualize n=100 simulations (3 rt curves, two sampling schemes, two sample sizes)
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
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ape))
source(here::here("src", "utility_functions.R"))
my_theme <- list(
  scale_fill_brewer(name = "CI Width",
                    labels = ~percent(as.numeric(.))),
  guides(fill = guide_legend(reverse = TRUE)),
  theme_minimal(),
  theme())
sim_num = 1
lump_val = 7
sim_dict = read_csv(here::here("data",
                               "sim_data",
                               "sim_dict.csv"))
# trees -------------------------------------------------------------------
sim_val = 2
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim2 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim2$time) - min(fitted_simdata_sim2$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim2 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim2 <- tree_sim2 %>%
  ggtree() + 
  ggtitle("Fixed Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim4
sim_val = 4
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim4 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim4$time) - min(fitted_simdata_sim4$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim4 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim4 <- tree_sim4 %>%
  ggtree() + 
  ggtitle("Fixed Heterochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 6
sim_val = 6
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim6 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim6$time) - min(fitted_simdata_sim6$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim6 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim6 <- tree_sim6 %>%
  ggtree() + 
  ggtitle("Increase Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))

# sim 8
sim_val = 8
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim8 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim8$time) - min(fitted_simdata_sim8$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim8 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim8 <- tree_sim8 %>%
  ggtree() + 
  ggtitle("Increase Heterochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 10
sim_val = 10
sim_num  = 3
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim10 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim10$time) - min(fitted_simdata_sim10$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim10 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim10 <- tree_sim10 %>%
  ggtree() + 
  ggtitle("Control Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 12
sim_val = 12
sim_num = 3
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim12 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim12$time) - min(fitted_simdata_sim12$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim12 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim12 <- tree_sim12 %>%
  ggtree() + 
  ggtitle("Control Heterochronous") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# rt plots ----------------------------------------------------------------
# sim 2
sim_num_val_val = 1
sim_val = 2
my_posterior_rt_sim2 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim2_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)

sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)

my_plot_posterior_rt_sim2 <- my_posterior_rt_sim2 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Fixed Isochronous") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("") +
  theme(text = element_text(size = 20),
        legend.position = "none")
# sim 4
sim_num_val_val = 1
sim_val = 4
my_posterior_rt_sim4 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim4_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)

sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)

my_plot_posterior_rt_sim4 <- my_posterior_rt_sim4 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Fixed Heterochronous") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("") +
  theme(text = element_text(size = 20),
        legend.position = c(0.75,0.9),
        legend.background = element_blank())
# sim 6
sim_num_val_val = 1
sim_val = 6
my_posterior_rt_sim6 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim6_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim6 <- my_posterior_rt_sim6 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Increase Isochronous") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("") +
  theme(text = element_text(size = 20),
        legend.position = "none")
# sim 8
sim_num_val_val = 1
sim_val = 8
my_posterior_rt_sim8 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim8_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim8 <- my_posterior_rt_sim8 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Increase Heterochronous") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("") +
  theme(text = element_text(size = 20),
        legend.position = "none")
# sim 10
sim_num_val_val = 3
sim_val = 10
my_posterior_rt_sim10 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim10_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim10 <- my_posterior_rt_sim10 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Control Isochronous") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")
# sim 12
sim_num_val_val = 3
sim_val = 12
my_posterior_rt_sim12 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim12_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)

sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim12 <- my_posterior_rt_sim12 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("Control Heterochronous") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")
# 100 samp plot ------------------------------------------------------------
fixed_plot <- (tree_plot_sim2 + 
                 my_plot_posterior_rt_sim2 + tree_plot_sim4 + 
                 my_plot_posterior_rt_sim4) + 
  plot_layout(ncol = 4, nrow = 1)
increase_plot <- (tree_plot_sim6 + my_plot_posterior_rt_sim6 + tree_plot_sim8 + 
                 my_plot_posterior_rt_sim8) + 
  plot_layout(ncol = 4, nrow = 1)
control_plot <- (tree_plot_sim10 + 
                 my_plot_posterior_rt_sim10 + tree_plot_sim12 + 
                 my_plot_posterior_rt_sim12) + 
  plot_layout(ncol = 4, nrow = 1)
samp100_plot <- fixed_plot / increase_plot / control_plot + 
  plot_annotation(title = "Simulated Trees and Posterior Ru Estimates (N=100)",
                               theme = theme(plot.title = element_text(size = 20)))
ggsave(here::here("figures", "samp100_plot.pdf"), 
       samp100_plot, 
       width = 20, height = 15,units = "in")