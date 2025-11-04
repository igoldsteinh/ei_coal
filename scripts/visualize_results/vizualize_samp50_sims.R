# visualize n=50 simulations 
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
sim_val = 1
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim1 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim1$time) - min(fitted_simdata_sim1$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim1 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))
samp_time = 153
tree_plot_sim1 <- tree_sim1 %>%
  ggtree() + 
  ggtitle("Fixed Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, samp_time)) +
  theme(text = element_text(size = 20))
# sim3
sim_val = 3
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim3 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim3$time) - min(fitted_simdata_sim3$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim3 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim3 <- tree_sim3 %>%
  ggtree() + 
  ggtitle("Fixed Heterochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 5
sim_val = 5
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim5 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim5$time) - min(fitted_simdata_sim5$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim5 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim5 <- tree_sim5 %>%
  ggtree() + 
  ggtitle("Increase Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 7
sim_val = 7
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim7 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim7$time) - min(fitted_simdata_sim7$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim7 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim7 <- tree_sim7 %>%
  ggtree() + 
  ggtitle("Increase Heterochronous") + 
  theme_tree2() + 
  scale_x_continuous("", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 9
sim_val = 9
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim9 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim9$time) - min(fitted_simdata_sim9$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim9 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim9 <- tree_sim9 %>%
  ggtree() + 
  ggtitle("Control Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# sim 11
sim_val = 11
tree_data_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata_sim11 <- read_csv(here::here("data", 
                                           "sim_data",
                                           tree_data_name), 
                                show_col_types = FALSE) %>%
  filter(sim == sim_num)
samp_time <- max(fitted_simdata_sim11$time) - min(fitted_simdata_sim11$time)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
tree_sim11 <- read.tree(here::here("data", 
                                  "sim_data", 
                                  paste0(sim_name, "_simnum", sim_num, ".tree")))

tree_plot_sim11 <- tree_sim11 %>%
  ggtree() + 
  ggtitle("Control Heterochronous") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))
# rt plots ----------------------------------------------------------------
# sim 1
sim_num_val_val = 1
sim_val = 1
my_posterior_rt_sim1 <- read_csv(here::here("results", 
                                                     "my_generated_quantities", 
                                                     "ei_cdf_sim1_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim1 <- my_posterior_rt_sim1 %>% 
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
# sim 3
sim_num_val_val = 1
sim_val = 3
my_posterior_rt_sim3 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim3_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)

sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim3 <- my_posterior_rt_sim3 %>% 
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
        legend.position = c(0.8,0.9),
        legend.background = element_blank())
# sim 5
sim_num_val_val = 1
sim_val = 5
my_posterior_rt_sim5 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim5_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim5 <- my_posterior_rt_sim5 %>% 
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
# sim 7
sim_num_val_val = 1
sim_val = 7
my_posterior_rt_sim7 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim7_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)

sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim7 <- my_posterior_rt_sim7 %>% 
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
# sim 9
sim_num_val_val = 1
sim_val = 9
my_posterior_rt_sim9 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim9_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim9 <- my_posterior_rt_sim9 %>% 
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
 # sim 11
sim_num_val_val = 1
sim_val = 11
my_posterior_rt_sim11 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim11_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
sim_name = sim_dict %>%
  filter(sim_id == sim_val) %>%
  pull(sim)
my_plot_posterior_rt_sim11 <- my_posterior_rt_sim11 %>% 
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
# 50 samp plot ------------------------------------------------------------
fixed_plot <- (tree_plot_sim1 + 
                 my_plot_posterior_rt_sim1 + tree_plot_sim3 + 
                 my_plot_posterior_rt_sim3) + 
  plot_layout(ncol = 4, nrow = 1)
increase_plot <- (tree_plot_sim5 + 
                 my_plot_posterior_rt_sim5 + tree_plot_sim7 + 
                 my_plot_posterior_rt_sim7) + 
  plot_layout(ncol = 4, nrow = 1)
control_plot <- (tree_plot_sim9 + 
                 my_plot_posterior_rt_sim9 + tree_plot_sim11 + 
                 my_plot_posterior_rt_sim11) + 
  plot_layout(ncol = 4, nrow = 1)
samp50_plot <- fixed_plot / increase_plot / control_plot + 
  plot_annotation(title = "Simulated Trees and Posterior Ru Estimates (N=50)",
                               theme = theme(plot.title = element_text(size = 20)))
ggsave(here::here( "figures", "samp50_plot.pdf"),
       samp50_plot,
       width = 20, height = 15, dpi = 300)
