library(DESeq2)
library(mclust)
library(seqinr)
library(tidyverse)

args <- commandArgs(TRUE)

print(args)
reference_fasta <- args[1]
reference_gff <- args[2]
features_of_interest <- args[3]
three_prime_trim <- as.integer(args[4]) / 100
output_prefix <- args[5]
counts_files <- args[6:length(args)]

# DEBUG
reference_fasta <- "k56.fasta"
reference_gff <- "k56.gff"
features_of_interest <- "gene"
three_prime_trim <- 10 / 100
output_prefix <- "testan"
counts_files <- "test.sites.tsv"

# read sites files
counts_data <- data.frame()
for (counts_file in counts_files) {
  if (!file.exists(counts_file)) { stop(paste(counts_file, "Not found!")) }
  counts_data <- read_tsv(counts_file) %>% mutate(counts_file = counts_file) %>% rbind(counts_data, .) %>% tbl_df
}
counts_data <- counts_data %>% arrange(counts_file, template, position)

# LOESS smoothing
smoothed <- counts_data %>% group_by(counts_file, template) %>%
  do(pred_num_reads = predict(loess(num_reads ~ position, data = ., span = 1,
        control = loess.control(statistics = c("approximate"), trace.hat = c("approximate"))), .$position))
counts_data$pred_num_reads <- lapply(smoothed$pred_num_reads, unlist) %>% unlist
counts_data <- counts_data %>% group_by(counts_file, template) %>%
  mutate(adjustment_ratio = pred_num_reads / median(pred_num_reads),
         smoothed_num_reads = num_reads / adjustment_ratio) %>% ungroup

# Visualize smoothing results
for (counts_file_entry in unique(counts_data$counts_file)) {
  counts_data %>% group_by(template) %>% filter(max(position) >= 10000) %>% # No short contigs
    # filter(counts_file == counts_file_entry) %>%
    ggplot(aes(position)) + theme_bw() + scale_x_continuous(labels = scales::comma) +
    geom_point(aes(y = num_reads), color = "#e41a1c", alpha = 0.2) +
    geom_point(aes(y = smoothed_num_reads), color = "#4daf4a", alpha = 0.2) + 
    geom_line(aes(y = pred_num_reads), color = "#377eb8") +
    facet_wrap(~template, scales = "free")
}
