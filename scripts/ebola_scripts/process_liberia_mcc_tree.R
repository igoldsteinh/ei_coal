library(ape)
library(phylodyn)
library(tidyverse)
suppressPackageStartupMessages(library(ggtree))

mcc_tree <- read.nexus(here::here("scripts", "ebola_scripts", "libera_mcc_tree.trees"))


lbr_list <- summarize_phylo(mcc_tree)
max(lbr_list[["coal_times"]])
tree_plot<- ggtree(mcc_tree, mrsd = max_date, as.Date = TRUE) + 
  #geom_tiplab(size = 2) +
  ggtitle("LBR Ebola MCC Tree") + 
  theme_tree2() + 
  xlab("Years") + 
  theme(text = element_text(size = 20))

sample_frame_lbr <- data.frame(samp_times = lbr_list$samp_times,
                               samp_lin = lbr_list$n_sampled) %>%
  mutate(samp_times_days = samp_times * 365.25,
         samp_times_4dig = round(samp_times_days, 4)) %>%
  group_by(samp_times_4dig) %>%
  summarise(samp_lin = sum(samp_lin))

coal_times_lbr <- data.frame(coal_times = lbr_list$coal_times) %>%
  mutate(coal_times_days = coal_times * 365.25)


max_reverse_time = max(coal_times_lbr$coal_times_days)
# extract last sample date
fields <- strsplit(mcc_tree$tip.label, "\\|")
# gather all entry 6 into one vector
sample_dates <- sapply(fields, function(x) x[6])
max_date = max(sample_dates)
min_samp_date = min(sample_dates)
reverse_times = seq(0, round(max_reverse_time), by = 1)
reverse_dates = seq(as.Date(max_date), by = "-1 day", length.out = length(reverse_times))
lbr_date_crosswalk <- data.frame(reverse_times = reverse_times,
                                 reverse_dates = reverse_dates,
                                 forward_times = abs(reverse_times - max(reverse_times)))
# save the dataframes
write_csv(sample_frame_lbr, here::here("data", "real_data", "ebola_lbrmcc_2014_samp_times.csv"))
write_csv(coal_times_lbr, here::here("data", "real_data", "ebola_lbrmcc_2014_coal_times.csv"))
write_csv(lbr_date_crosswalk, here::here("data", "real_data", "ebola_lbrmcc_2014_date_crosswalk.csv"))

