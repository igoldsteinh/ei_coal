#v2 with updated tree
# where we discarded first 500 trees as burnin
# first 5000000 states
library(ape)
library(phylodyn)
library(tidyverse)
library(treeio)
suppressPackageStartupMessages(library(ggtree))

mcc_tree <- read.nexus(here::here("scripts", "ebola_scripts", "liberia_mcctreev2"))
# add starting state
mcc_tree[["tip.label"]] <- paste0("I|", mcc_tree[["tip.label"]])

# Clean up the tip labels
mcc_tree$tip.label <- gsub("\\|", "_", mcc_tree$tip.label)  # Replace | with _
mcc_tree$tip.label <- gsub("\\?", "", mcc_tree$tip.label)   # Remove ?

# save to nexus
write.nexus(mcc_tree, file=here::here("data", "real_data", "liberia_mcctreev2.nexus"))

# get the dates
# If your tree tips are in tree$tip.label
taxon_names <- mcc_tree$tip.label
# Function to convert date string to decimal year
date_to_decimal_year <- function(date_string) {
  date_obj <- as.Date(date_string, format = "%Y-%m-%d")
  year <- year(date_obj)
  start_of_year <- as.Date(paste0(year, "-01-01"))
  start_of_next_year <- as.Date(paste0(year + 1, "-01-01"))
  year_length <- as.numeric(start_of_next_year - start_of_year)
  day_of_year <- as.numeric(date_obj - start_of_year)
  decimal_year <- year + (day_of_year / year_length)
  return(decimal_year)
}


date_strings <- str_extract(taxon_names, "\\d{4}-\\d{2}-\\d{2}")
decimal_years <- sapply(date_strings, date_to_decimal_year)
trait_pairs <- paste0(taxon_names, "=", sprintf("%.10f", decimal_years))
trait_value <- paste(trait_pairs, collapse = ",")

# Print to copy-paste into XML
cat(trait_value)

# Using ape package
library(ape)

# Get tip ages (time before present)
tip_ages <- max(node.depth.edgelength(mcc_tree)) - node.depth.edgelength(mcc_tree)[1:Ntip(mcc_tree)]
range(tip_ages)
summary(tip_ages)

# Or to see all tip ages
tip_ages
lbr_list <- summarize_phylo(mcc_tree)
max(lbr_list[["coal_times"]])

sample_frame_lbr <- data.frame(samp_times = lbr_list$samp_times,
                               samp_lin = lbr_list$n_sampled) %>%
  mutate(samp_times_days = samp_times * 365.25,
         samp_times_4dig = round(samp_times_days, 4)) %>%
  group_by(samp_times_4dig) %>%
  summarise(samp_lin = sum(samp_lin))

coal_times_lbr <- data.frame(coal_times = lbr_list$coal_times) %>%
  mutate(coal_times_days = coal_times * 365.25)

check <- coal_times_lbr %>%
  mutate(time_diff = coal_times_days - lag(coal_times_days)) %>%
  filter(time_diff < 0)

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
write_csv(sample_frame_lbr, here::here("data", "real_data", "ebola_lbrmcc_2014_samp_timesv2.csv"))
write_csv(coal_times_lbr, here::here("data", "real_data", "ebola_lbrmcc_2014_coal_timesv2.csv"))
write_csv(lbr_date_crosswalk, here::here("data", "real_data", "ebola_lbrmcc_2014_date_crosswalkv2.csv"))

tree_plot<- ggtree(mcc_tree, mrsd = max_date, as.Date = TRUE) + 
  #geom_tiplab(size = 2) +
  ggtitle("LBR Ebola MCC Tree") + 
  theme_tree2() + 
  xlab("Years") + 
  theme(text = element_text(size = 20))

mcc_tree <- read.beast(here::here("scripts", "ebola_scripts", "liberia_mcctreev2"))

tree_plot<- ggtree(mcc_tree, mrsd = max_date, as.Date = TRUE) + 
  #geom_tiplab(size = 2) +
  aes(color = posterior) + 
  scale_color_gradient(low = "red", high = "blue", name = "Posterior") +
  ggtitle("Liberia (2014 - 2015) MCC Tree") + 
  theme_tree2() + 
  xlab("Time") + 
  theme(text = element_text(size = 20))

