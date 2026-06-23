# Visualize bdmm with 5 intervals vs 22 intervals
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

# bdmm  -------------------------------------------------------------------
bdmm_rt9_5 <- read_csv(here::here("scripts", 
                                "BEAST2", 
                                "bdmm_sim9_allseeds_rt_quantilesv2.csv")) %>%
  filter(sim_num_val == sim_num_val_val) %>%
  filter(address != "beast_controliso50_simnum1_oldversion.log")
bdmm_plot_posterior_rt_sim9_5 <- bdmm_rt9_5 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("BDMM 5") + 
  my_theme + 
  ylim(c(0, 7)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")

bdmm_rt9_22 <- read_csv(here::here("scripts", 
                                  "BEAST2", 
                                  "bdmm_sim9_allseeds_rt_quantiles_weekly.csv")) %>%
  filter(sim_num_val == sim_num_val_val)
bdmm_plot_posterior_rt_sim9_22 <- bdmm_rt9_22 %>% 
  filter(lump >= 0) %>%
  ggplot(aes(x = new_time, y = value,  ymin = .lower, ymax = .upper)) +
  geom_lineribbon() +
  geom_point(aes(x = new_time, y = true_rt), color = "orange") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_bw() +
  ggtitle("BDMM 22") + 
  my_theme + 
  ylim(c(0, 7)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")

# 50 samp plot ------------------------------------------------------------
compare_plot <- (bdmm_plot_posterior_rt_sim9_5 + bdmm_plot_posterior_rt_sim9_22) + 
  plot_annotation(title = "BDMM with 5 vs 22 Intervals",
                               theme = theme(plot.title = element_text(size = 20)))
ggsave(here::here( "figures", "bddm5vs22_plot.pdf"),
       compare_plot,
       width = 8, height = 4, dpi = 300)
