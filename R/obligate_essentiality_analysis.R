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
ignore_sites <- as.integer(args[5])
min_reads_per_site <- as.integer(args[6])
output_prefix <- args[7]
counts_files <- args[8:length(args)]

# DEBUG
reference_fasta <- "k56.fasta"
reference_gff <- "k56.gff"
features_of_interest <- "gene"
three_prime_trim <- 10 / 100
ignore_sites <- 0
min_reads_per_site <- 1
output_prefix <- "testan"
counts_files <- "test.sites.tsv"

# read sites files
counts_data <- data.frame()
for (counts_file in counts_files) {
  if (!file.exists(counts_file)) { stop(paste(counts_file, "Not found!")) }
  counts_data <- read_tsv(counts_file) %>%
    arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>% # Cut off i highest sites
    filter(num_reads >= min_reads_per_site) %>% # at least m reads/site
    mutate(counts_file = counts_file) %>% rbind(counts_data, .) %>% tbl_df
}
counts_data <- counts_data %>% arrange(counts_file, template, position)

# LOESS smoothing
smoothed <- counts_data %>% group_by(counts_file, template) %>%
#  filter(n() > 1) %>% # Can't smooth a single number...
  do(pred_num_reads = predict(loess(num_reads ~ position, data = ., span = 1,
        control = loess.control(statistics = c("approximate"), trace.hat = c("approximate"))), .$position))
counts_data$pred_num_reads <- lapply(smoothed$pred_num_reads, unlist) %>% unlist
counts_data <- counts_data %>% group_by(counts_file, template) %>%
  mutate(adjustment_ratio = pred_num_reads / median(pred_num_reads),
         smoothed_num_reads = num_reads / adjustment_ratio) %>% ungroup

reference_genome <- read.fasta(reference_fasta)
template_sizes <- data.frame(template = names(reference_genome),
                             template_length = unlist(lapply(reference_genome, length)),
                             stringsAsFactors = F) %>% tbl_df

# Visualize smoothing results
for (counts_file_entry in unique(counts_data$counts_file)) {
  p <- template_sizes %>% filter(template_length >= 20000) %>% inner_join(counts_data) %>% # No short contigs
    filter(counts_file == counts_file_entry) %>%
    group_by(template) %>% sample_frac(0.05) %>% # Only draw 5% of points for decluttering
    ggplot(aes(position)) + theme_bw() + scale_x_continuous(labels = scales::comma) +
    geom_point(aes(y = num_reads, color = "Raw"), alpha = 0.5) +
    geom_point(aes(y = smoothed_num_reads, color = "Smoothed"), alpha = 0.5) + 
    geom_segment(aes(xend = position, y = num_reads, yend = smoothed_num_reads), color = "grey") +
    geom_line(aes(y = pred_num_reads, color = "LOESS Trend")) +
    facet_wrap(~template, scales = "free") +
    scale_color_brewer("", type = "qual", palette = "Dark2") + theme(legend.position = "bottom") +
    ggtitle("LOESS smoothing positional bias", counts_file_entry) + ylab("Number of reads") + xlab("Position")
  ggsave(paste0(counts_file_entry, ".loess.png"), p, scale = 2, height = 4, width = 6, units = "in")
}

# Read in features
genome_features <- read_tsv(reference_gff, comment = "#",
                            col_names = c("template", "source", "feature", "start", "end", "score", "strand",
                                          "frame", "attribute"))
