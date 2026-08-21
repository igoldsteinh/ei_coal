# visualize empirical EI vs coal EI intercoalescence times
# across multiple R0 values (original, r=1.5, r=1.3)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(phylodyn))
suppressPackageStartupMessages(library(patchwork))

# original (r0 = 2) --------------------------------------------------------
imp_res <- read_csv(here::here("data", "compare_data", "compare_ei_impcoal_results.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "Empirical")
tt_res_noneg <- read_csv(here::here("data", "compare_data", "compare_ei_tt_badpop_results.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT")
tt_res_ode <- read_csv(here::here("data", "compare_data", "compare_ei_tt_odebadpop_results.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT ODE")

box_plot_orig <- imp_res %>%
  bind_rows(tt_res_noneg) %>%
  bind_rows(tt_res_ode) %>%
  mutate(interval = ifelse(num == "diff1", "Interval 1",
                          ifelse(num == "diff2", "Interval 2",
                                 ifelse(num == "diff3", "Interval 3", "Interval 4")))) %>%
  ggplot(aes(x = as.factor(interval), y = coal_int, fill = Algorithm)) +
  geom_boxplot() +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = c(0.8,0.9),
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Interval", y = "Intercoalescence Time") +
  ggtitle("R0 = 2")

# r0 = 1.5 ------------------------------------------------------------------
imp_res_r15 <- read_csv(here::here("data", "compare_data", "compare_ei_impcoal_results_r=1.5.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "Empirical")
tt_res_noneg_r15 <- read_csv(here::here("data", "compare_data", "compare_ei_tt_badpop_results_r=1.5.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT")
tt_res_ode_r15 <- read_csv(here::here("data", "compare_data", "compare_ei_tt_odebadpop_results_r=1.5.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT ODE")

box_plot_r15 <- imp_res_r15 %>%
  bind_rows(tt_res_noneg_r15) %>%
  bind_rows(tt_res_ode_r15) %>%
  mutate(interval = ifelse(num == "diff1", "Interval 1",
                          ifelse(num == "diff2", "Interval 2",
                                 ifelse(num == "diff3", "Interval 3", "Interval 4")))) %>%
  ggplot(aes(x = as.factor(interval), y = coal_int, fill = Algorithm)) +
  geom_boxplot() +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Interval", y = "Intercoalescence Time") +
  ggtitle("R0 = 1.5")

# r0 = 1.3 ------------------------------------------------------------------
imp_res_r13 <- read_csv(here::here("data", "compare_data", "compare_ei_impcoal_results_r=1.3.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "Empirical")
tt_res_noneg_r13 <- read_csv(here::here("data", "compare_data", "compare_ei_tt_badpop_results_r=1.3.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT")
tt_res_ode_r13 <- read_csv(here::here("data", "compare_data", "compare_ei_tt_odebadpop_results_r=1.3.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT ODE")

box_plot_r13 <- imp_res_r13 %>%
  bind_rows(tt_res_noneg_r13) %>%
  bind_rows(tt_res_ode_r13) %>%
  mutate(interval = ifelse(num == "diff1", "Interval 1",
                          ifelse(num == "diff2", "Interval 2",
                                 ifelse(num == "diff3", "Interval 3", "Interval 4")))) %>%
  ggplot(aes(x = as.factor(interval), y = coal_int, fill = Algorithm)) +
  geom_boxplot() +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Interval", y = "Intercoalescence Time") +
  ggtitle("R0 = 1.3")

# patchwork of all three -----------------------------------------------------
patchwork_plot <- (box_plot_orig + box_plot_r15 + box_plot_r13) +
  plot_annotation(title = "Empirical vs Time Transformation Varying R0")
ggsave(here::here("figures", "empEI_vs_coalEI_allR0.pdf"),
       patchwork_plot, width = 14, height = 6)
