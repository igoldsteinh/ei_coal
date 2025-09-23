# visualize an example of the model failing 
# along with tree which results in failure 
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
sim_num = 13
lump_val = 7
sim_dict = read_csv(here::here("data",
                               "sim_data",
                               "sim_dict.csv"))
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
  ggtitle("Fixed Iso") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, samp_time)) +
  theme(text = element_text(size = 20))

sim_num_val_val = 13
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
  ggtitle("Fixed Iso") + 
  my_theme + 
  ylim(c(0, 4.5)) +
  ylab("Ru") +
  xlab("Forward Time") +
  theme(text = element_text(size = 20),
        legend.position = "none")
fail_plot <- tree_plot_sim1 + my_plot_posterior_rt_sim1
ggsave(here::here(
                  "figures", 
                  "example_failure_sim1.pdf"), 
       fail_plot, 
       width = 8, 
       height = 4)
