# visualize empirical EI vs coal EI intercoalescence times
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(phylodyn))
suppressPackageStartupMessages(library(patchwork))
# visualize empirical EI vs EI coal intercoalescence times 
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
box_plot <- imp_res %>%
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
  ggtitle("Empirical vs Time Transformation") 
# save plot
ggsave(here::here("figures", "empEI_vs_coalEI.pdf"), box_plot, width = 10, height = 6)
