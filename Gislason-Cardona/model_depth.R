library(seqinr)
library(reshape2)
library(tidyverse)

setwd("~/cenocepacia-tnseq/")
reference_fasta <- "k56.fasta"
ignore_sites <- 100
min_reads_per_site <- 2
counts_file <- "1mhdtm_taq.sites.tsv"

# read sites files
counts_data <- read_tsv(counts_file) %>%
  arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>% # Cut off i highest sites
  filter(num_reads >= min_reads_per_site) # at least m reads/site

# Genome size
reference_genome <- read.fasta(reference_fasta)
template_sizes <- tibble(template = names(reference_genome),
                         template_length = unlist(lapply(reference_genome, length)))

# GC Content
get_window_gc <- function(template, pos) {
  x <- reference_genome[[template]]
  start <- max(1, pos - 500)
  end <- min(length(x), pos + 500)
  x[start:end] %>% GC
}
counts_with_gc <- counts_data %>% mutate(gc_content = map2(template, position, get_window_gc)) %>% unnest(gc_content)

# Visualization
counts_with_gc %>% ggplot(aes(position, num_reads, color = gc_content)) + geom_point() +
  facet_wrap(~template, scales = "free") +
  scale_color_distiller(palette = "Spectral") +
  scale_y_log10() +
  theme_bw() 

# Can we model it?
library(broom)
# contig-position interaction, contig effect, global gc effect, no global position effect
counts_model <- lm(log(num_reads) ~ gc_content + poly(position, 3) * template - poly(position, 3), counts_with_gc)
counts_with_gc %>% add_predictions(counts_model) %>%
  mutate(pred_num_reads = exp(pred)) %>%
  ggplot(aes(position)) +
  geom_point(aes(y = num_reads, color = gc_content)) +
  geom_line(aes(y = pred_num_reads), color = "blue") +
  facet_wrap(~template, scales = "free") +
  scale_color_distiller(palette = "Spectral") + theme_bw() +
  scale_y_log10() + annotation_logticks(sides = "l")

# Correct for model prediction
counts_with_gc %>% add_predictions(counts_model) %>%
  mutate(lnr = log(num_reads),
         correction_ratio = pred / mean(pred),
         smoothed_num_reads = exp(lnr / correction_ratio)) %>%
  ggplot(aes(position, color = gc_content)) +
  geom_point(aes(y = smoothed_num_reads)) +
  scale_y_log10() + annotation_logticks(sides = "l") +
  facet_wrap(~template, scales = "free") +
  scale_color_distiller(palette = "Spectral") + theme_bw()
