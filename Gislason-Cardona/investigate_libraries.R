library(seqinr)
library(reshape2)
library(tidyverse)

setwd("~/cenocepacia-s3/")
reference_fasta <- "k56.fasta"
reference_gff <- "k56.bycontig.gff"
features_of_interest <- "gene"
three_prime_trim <- 10 / 100
ignore_sites <- 0
min_reads_per_site <- 2
output_prefix <- "testan"
counts_files <- c("500khdtm_14.sites.tsv", "500khdtm_15.sites.tsv",
                  "1mhdtm_taq.sites.tsv", "1mhdtm_kapa.sites.tsv")
metadata <- data.frame(counts_file = counts_files,
                       lib_number = c(rep("one", 2), rep("two", 2)),
                       polymerase = c(rep("taq", 3), "kapa")) %>% tbl_df %>%
  mutate(counts_file = as.character(counts_file))

# read sites files
counts_data <- data.frame()
for (counts_file in counts_files) {
  if (!file.exists(counts_file)) { stop(paste(counts_file, "Not found!")) }
  counts_data <- read_tsv(counts_file) %>%
    arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>% # Cut off i highest sites
    filter(num_reads >= min_reads_per_site) %>% # at least m reads/site
    mutate(counts_file = counts_file) %>% rbind(counts_data, .) %>% tbl_df
}
counts_data <- counts_data %>% arrange(counts_file, template, position) %>%
  inner_join(metadata)

# LOESS smoothing
smoothed <- counts_data %>% group_by(counts_file, template) %>%
  #  filter(n() > 1) %>% # Can't smooth a single number...
  do(pred_num_reads = predict(loess(num_reads ~ position, data = ., span = 1,
                                    control = loess.control(statistics = c("approximate"), trace.hat = c("approximate"))), .$position))
counts_data$pred_num_reads <- lapply(smoothed$pred_num_reads, unlist) %>% unlist
counts_data <- counts_data %>% group_by(counts_file, template) %>%
  mutate(adjustment_ratio = pred_num_reads / median(pred_num_reads),
         smoothed_num_reads = num_reads / adjustment_ratio) %>% ungroup

# Genome size
reference_genome <- read.fasta(reference_fasta)
template_sizes <- data.frame(template = names(reference_genome),
                             template_length = unlist(lapply(reference_genome, length)),
                             stringsAsFactors = F) %>% tbl_df

# Are insertion sites evenly distributed across genome?
template_sizes %>% filter(template_length >= 20000) %>% inner_join(counts_data) %>%
  mutate(position = position / 1000) %>%
  # filter(num_reads > 1) %>%
  ggplot(aes(position, color = counts_file)) +
  geom_freqpoly(binwidth = 10) + 
  theme_bw() + facet_wrap(~template, scales = "free_x") +
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  ylab("Insertions per 10kb") + xlab("Position (kb)")

# Does it track with GC content?
sliding_window_gc <- data.frame()
for (y in 1:length(reference_genome)) {
  x <- reference_genome[[y]]
  starts <- seq(1, length(x) - 2000, by = 2000)
  chunkGCs <- numeric(length(starts))
  chunkGCs
  for (i in 1:length(starts)) {
    chunkGCs[i] <- x[starts[i]:(starts[i]+1999)] %>% GC
  }
  sliding_window_gc <- data.frame(template = names(reference_genome)[y],
                                  position = starts + 1000, gc_content = chunkGCs,
                                  stringsAsFactors = F) %>%
    rbind(sliding_window_gc, .) %>% tbl_df
}

sliding_window_gc <- sliding_window_gc %>% mutate(polymerase = "kapa") %>%
  rbind(mutate(sliding_window_gc, polymerase = "taq")) %>% tbl_df

template_sizes %>% filter(template_length >= 20000) %>% inner_join(counts_data) %>%
  mutate(position = position / 1000) %>%
  # filter(num_reads > 1) %>%
  ggplot() +
  geom_rect(aes(xmin = position - 1, xmax = position + 1,
                ymin = 0, ymax = 3500, fill = gc_content),
            alpha = 0.5,
            data = mutate(inner_join(filter(template_sizes, template_length >= 20000),
                                     sliding_window_gc),
                          position = position / 1000)) +
  geom_freqpoly(aes(position, color = counts_file), size = 2, binwidth = 10) + 
  theme_bw() +
  # facet_wrap(~template, scales = "free_x") +
  facet_grid(polymerase ~ template, scales = "free_x", space = "free_x") +
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  scale_fill_distiller(palette = "Greys") +
  ylab("Insertions per 10kb") + xlab("Position (kb)")

template_sizes %>% filter(template_length >= 20000) %>% inner_join(sliding_window_gc) %>%
  mutate(position = position / 1000) %>%
  # filter(num_reads > 1) %>%
  ggplot(aes(position, gc_content)) + geom_line() +
  facet_wrap(~template, scales = "free_x") + theme_bw() +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::percent)

# Should combine these plots somehow