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

# insert synthetic rows every gap_size units wherever new_time jumps by more than
# gap_size, carrying forward the earlier row's values, so plots don't show big
# gaps between dots. each new_time has one row per .width (the CI ribbons), so
# gaps are computed within each .width separately and all rows for a given
# synthetic time are filled together
fill_new_time_gaps <- function(df, gap_size = 1) {
  df %>%
    group_by(.width) %>%
    group_modify(~ {
      d <- arrange(.x, new_time)
      out <- d
      for (i in seq_len(nrow(d) - 1)) {
        gap <- d$new_time[i + 1] - d$new_time[i]
        n_fill <- floor((gap - 1e-9) / gap_size)
        if (n_fill > 0) {
          fill_times <- d$new_time[i] + gap_size * seq_len(n_fill)
          new_rows <- d[rep(i, length(fill_times)), ]
          new_rows$new_time <- fill_times
          out <- bind_rows(out, new_rows)
        }
      }
      arrange(out, new_time)
    }) %>%
    ungroup()
}
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
  ggtitle("Fixed Isochronous") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, samp_time)) +
  theme(text = element_text(size = 20))

sim_num_val_val = 13
sim_val = 1
my_posterior_rt_sim1 <- read_csv(here::here("results", 
                                            "my_generated_quantities", 
                                            "ei_cdf_sim1_allseeds_rt_quantiles.csv")) %>%
  filter(sim_num_val == sim_num_val_val) %>%
  fill_new_time_gaps()

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
