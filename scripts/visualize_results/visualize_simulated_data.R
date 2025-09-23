# plot simulated SEIR trajectory, Rt, and sampled tree 
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
# if i just type "sim" here this will error I don't know why
sim_val = 1
sim_num_val = 1
state_sim_val = 1
true_samp_time = 153
lump_val = 7
sim_dict = read_csv(here::here("data", "sim_data", "sim_dict.csv"), show_col_types = FALSE)
tree_file_name = sim_dict %>% 
  filter(sim_id == sim_val) %>%
  pull(tree_file)
fitted_simdata <- read_csv(here::here("data", 
                                      "sim_data", 
                                      tree_file_name), show_col_types = FALSE) %>%
  filter(sim == sim_num_val)
samp_time <- max(fitted_simdata$time)
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
true_rt_plot <- state_frame %>% 
  dplyr::select(new_time, rt) %>%
  ggplot(aes(x = new_time, y = rt)) +
  geom_point(color = "orange") +
  theme_bw() +
  xlab("Forward Time") + 
  ylab("Ru") + 
  ggtitle(paste0("True Ru")) +
  theme(text = element_text(size = 20))
# data massaging
true_rt <- state_frame %>% 
  rename("true_rt" = "rt") %>%
  dplyr::select(time, new_time, reverse_time, lump, true_rt) 
traj_data <- state_frame %>% 
  dplyr::select(time, new_time, reverse_time, lump, S, E, I, R)
true_start <- abs(samp_time - true_samp_time)
param_change_times <- c(0, rev(abs((fitted_simdata %>% filter(event == "alpha") %>% pull(time) - samp_time))))
reverse_param_change_times <- abs(param_change_times - samp_time)
solve_times <- fitted_simdata$time
num_timepoints <- length(solve_times)
time_crosswalk <- data.frame(index = 0:(num_timepoints-1), time = solve_times)
rt_lump_crosswalk <- data.frame(index = 1:length(reverse_param_change_times), 
                                reverse_param_change_times = reverse_param_change_times, 
                                param_change_times = param_change_times,
                                lump = 0:(length(reverse_param_change_times)-1))
# plot comps --------------------------------------------------------------
true_init_rt = true_rt %>% filter(new_time == min(new_time)) %>% pull(true_rt)
traj_plot <- traj_data %>%
  mutate(grid_time = ifelse(new_time < max(new_time), ceiling(new_time * 2)/2, new_time),
         diff_time = abs(grid_time - new_time)) %>%
  group_by(grid_time) %>%
  filter(diff_time == min(diff_time)) %>%
  ungroup() %>%
  dplyr::select(new_time, 
                S,
                E,
                I, 
                R) %>% 
  pivot_longer(-new_time, names_to = "Compartment") %>% 
  ggplot(aes(x = new_time, y = value, color = `Compartment`)) + 
  geom_point() + 
  geom_line() + 
  theme_bw() +
  xlab("Forward Time") + 
  ylab("Counts") + 
  ggtitle(paste0("Simulated Epidemic")) +
  theme(text = element_text(size = 20),
              legend.position = c(0.3,0.5),
              legend.background = element_blank())
# trees -------------------------------------------------------------------
sim_val = 1
sim_num = 1
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
  ggtitle("Sampled Phylogeny") + 
  theme_tree2() + 
  scale_x_continuous("Forward Time", limits = c(0, samp_time)) +
  theme(text = element_text(size = 20))
sim_plot <- traj_plot + true_rt_plot + tree_plot_sim1
ggsave(here::here( "figures", "simulated_data_plot.pdf"),
       sim_plot,
       width = 12, height = 4, dpi = 300)
