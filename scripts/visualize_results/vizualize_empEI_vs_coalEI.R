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
tt_res_ode_vargrid <- read_csv(here::here("data", "compare_data", "compare_ei_tt_odebadpop_vargrid_results.csv")) %>%
  dplyr::select(iteration, diff1, diff2, diff3, diff4) %>%
  group_by(iteration) %>%
  pivot_longer(cols = starts_with("diff"), names_to = "num", values_to = "coal_int") %>%
  mutate(Algorithm = "TT ODE Var")


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

# ode vs traj ode ---------------------------------------------------------
box_plot_ode <- imp_res %>%
  # bind_rows(tt_res) %>%
  # bind_rows(ttnaive_res) %>%
  # bind_rows(tt_res_noneg) %>%
  bind_rows(tt_res_ode) %>%
  bind_rows(tt_res_ode_vargrid) %>%
  filter(Algorithm != "Empirical") %>%
  # filter(coal_int > 0) %>%
  mutate(interval = ifelse(num == "diff1", "Interval 1", 
                           ifelse(num == "diff2", "Interval 2", 
                                  ifelse(num == "diff3", "Interval 3", "Interval 4")))) %>%
  ggplot(aes(x = as.factor(interval), y = coal_int, fill = Algorithm)) +
  geom_boxplot() +
  theme_minimal() +
  # scale_fill_grey(  start = 0.4,
  #                   end = 0.8) +
  theme(text = element_text(size = 20),
        legend.position = c(0.8,0.9),
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Interval", y = "Intercoalescence Time") + 
  ggtitle("ODE vs ODE Varying Grid") 

# save plot
ggsave(here::here("figures", "ode_vs_varODE.pdf"), box_plot_ode, width = 10, height = 6)


# ODE plot ----------------------------------------------------------------

ode_frame <- read_csv(here::here("data", "compare_data", "compare_ei_ode_state_frame.csv"))
stoch_frame <- read_csv(here::here("data", "compare_data", "compare_ei_state_frame.csv"))

set.seed(14578)
randnums <- sample.int(1000, 10)
subset_stoch_frame <- stoch_frame %>% filter(sim %in% randnums) %>%
  mutate(Type = "Stoch") %>%
  dplyr::select(time, E, I, Type, sim)
ode_frame <- ode_frame %>% 
  mutate(Type = "ODE", 
         sim = 101) %>%
  dplyr::select(-reverse_time)

count_plot <- ode_frame %>%
  bind_rows(subset_stoch_frame) %>% 
  group_by(sim, Type) %>%
  pivot_longer(-c(sim, Type,time),  names_to = "Compartment") %>%
  filter(Type != "ODE") %>%
  ggplot(aes(x = time, y = value, group = sim, color = Type)) + 
  geom_line(alpha = 0.5) +
  geom_line(data = ode_frame %>% pivot_longer(-c(sim, Type,time),  names_to = "Compartment"),
            aes(x = time, y = value, group = sim, color = Type)) +
  facet_wrap(vars(Compartment)) + 
  xlab("Forward Time") + 
  ylab("Counts") + 
  theme_minimal() +
  # scale_fill_grey(  start = 0.4,
  #                   end = 0.8) +
  theme(text = element_text(size = 20),
        legend.position = c(0.75,0.85),
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("ODE vs Stochastic EI Trajectories")
ggsave(here::here("figures", "odeEI_vs_stochEI.pdf"), count_plot, width = 10, height = 6)

