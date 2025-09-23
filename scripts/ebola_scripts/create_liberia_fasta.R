# ebola fasta file
# split off Liberia sequences 
library(tidyverse)
library(ape)
# read in fasta file ------------------------------------------------------
full_data <- read.dna(file=here::here("data", "real_data", "Makona_1610_genomes_2016-06-23.fasta"), 
                     format = "fasta") 

labels <- names(full_data)
# Extract locations
locations <- sapply(strsplit(labels, "\\|"), function(x) x[4])
# filter for locations == LBR
liberia_data <- full_data[grep("LBR", locations)]
# start and end date from Tang et al analysis
start_date <- as.Date("2014-06-20")
end_date <- as.Date("2015-02-14")
# Find sequences that are NOT in the date range
new_labels = names(liberia_data)
to_drop <- sapply(new_labels, function(x) {
  fields <- strsplit(x, "\\|")[[1]]
  tip_date <- as.Date(fields[6])
  tip_date < start_date || tip_date > end_date
})
# drop the entries of liberia_data which are true
liberia_data <- liberia_data[!to_drop]
# save the liberia data
write.dna(liberia_data, file=here::here("data", "real_data", "liberia_sequences.fasta"), format = "fasta")