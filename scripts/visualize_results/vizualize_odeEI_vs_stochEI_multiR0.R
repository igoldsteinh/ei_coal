# visualize ODE vs stochastic EI trajectories
# across multiple R0 values (original, r=1.5, r=1.3)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(phylodyn))
suppressPackageStartupMessages(library(patchwork))

# original (r0 = 2) --------------------------------------------------------
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

count_plot_orig <- ode_frame %>%
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
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("R0 = 2")

# r0 = 1.5 ------------------------------------------------------------------
ode_frame_r15 <- read_csv(here::here("data", "compare_data", "compare_ei_ode_state_frame_r=1.5.csv"))
stoch_frame_r15 <- read_csv(here::here("data", "compare_data", "compare_ei_state_frame_r=1.5.csv"))

set.seed(14578)
randnums_r15 <- sample.int(1000, 10)
subset_stoch_frame_r15 <- stoch_frame_r15 %>% filter(sim %in% randnums_r15) %>%
  mutate(Type = "Stoch") %>%
  dplyr::select(time, E, I, Type, sim)
ode_frame_r15 <- ode_frame_r15 %>%
  mutate(Type = "ODE",
         sim = 101) %>%
  dplyr::select(-reverse_time)

count_plot_r15 <- ode_frame_r15 %>%
  bind_rows(subset_stoch_frame_r15) %>%
  group_by(sim, Type) %>%
  pivot_longer(-c(sim, Type,time),  names_to = "Compartment") %>%
  filter(Type != "ODE") %>%
  ggplot(aes(x = time, y = value, group = sim, color = Type)) +
  geom_line(alpha = 0.5) +
  geom_line(data = ode_frame_r15 %>% pivot_longer(-c(sim, Type,time),  names_to = "Compartment"),
            aes(x = time, y = value, group = sim, color = Type)) +
  facet_wrap(vars(Compartment)) +
  xlab("Forward Time") +
  ylab("Counts") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = "none",
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("R0 = 1.5")

# r0 = 1.3 ------------------------------------------------------------------
ode_frame_r13 <- read_csv(here::here("data", "compare_data", "compare_ei_ode_state_frame_r=1.3.csv"))
stoch_frame_r13 <- read_csv(here::here("data", "compare_data", "compare_ei_state_frame_r=1.3.csv"))

set.seed(14578)
randnums_r13 <- sample.int(1000, 10)
subset_stoch_frame_r13 <- stoch_frame_r13 %>% filter(sim %in% randnums_r13) %>%
  mutate(Type = "Stoch") %>%
  dplyr::select(time, E, I, Type, sim)
ode_frame_r13 <- ode_frame_r13 %>%
  mutate(Type = "ODE",
         sim = 101) %>%
  dplyr::select(-reverse_time)

count_plot_r13 <- ode_frame_r13 %>%
  bind_rows(subset_stoch_frame_r13) %>%
  group_by(sim, Type) %>%
  pivot_longer(-c(sim, Type,time),  names_to = "Compartment") %>%
  filter(Type != "ODE") %>%
  ggplot(aes(x = time, y = value, group = sim, color = Type)) +
  geom_line(alpha = 0.5) +
  geom_line(data = ode_frame_r13 %>% pivot_longer(-c(sim, Type,time),  names_to = "Compartment"),
            aes(x = time, y = value, group = sim, color = Type)) +
  facet_wrap(vars(Compartment)) +
  xlab("Forward Time") +
  ylab("Counts") +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = c(0.75,0.85),
        legend.background = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("R0 = 1.3")

# patchwork of all three -----------------------------------------------------
patchwork_plot <- (count_plot_orig + count_plot_r15 + count_plot_r13) +
  plot_annotation(title = "ODE vs Stochastic EI Trajectories Varying R0")
ggsave(here::here("figures", "odeEI_vs_stochEI_allR0.pdf"),
       patchwork_plot, width = 18, height = 6)
