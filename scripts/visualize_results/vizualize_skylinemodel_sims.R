# code used to visualize skyline prior comparison
# between EI Coalescent and bdmm
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
  ggtitle("Control Iso") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, ceiling(samp_time))) +
  theme(text = element_text(size = 20))

# eicoal ----------------------------------------------------------------
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
  ggtitle("EI Coal") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")


# eicoal_skyline ------------------------------------------------------------------

eicoal_skyline <- read_csv(here::here("results", 
                                      "my_generated_quantities", 
                                  "ei_cdf_sim_skyline9_allseeds_rt_quantiles_local.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
eicoal_skyline_plot_posterior_rt_sim9 <- eicoal_skyline %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("EI Coal Skyline") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")


# bdmm  -------------------------------------------------------------------
bdmm_rt9 <- read_csv(here::here("scripts", 
                                "BEAST2", 
                                "bdmm_sim9_allseeds_rt_quantilesv2.csv")) %>%
  filter(sim_num_val == sim_num_val_val) %>%
  filter(address != "beast_controliso50_simnum1_oldversion.log")
bdmm_plot_posterior_rt_sim9 <- bdmm_rt9 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("BDMM") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")


# 50 samp plot ------------------------------------------------------------
control_plot <- (tree_plot_sim9 + 
                 my_plot_posterior_rt_sim9 +
                   bdmm_plot_posterior_rt_sim9 + eicoal_skyline_plot_posterior_rt_sim9) + 
  plot_annotation(title = "Skyline vs GMRF RT Models",
                               theme = theme(plot.title = element_text(size = 20)))
ggsave(here::here( "figures", "control50_skylinemodels_plot.pdf"),
       control_plot,
       width = 15, height = 15, dpi = 300)
