library(seqinr)
library(reshape2)
library(tidyverse)

setwd("~/cenocepacia-tnseq/")
reference_fasta <- "k56.fasta"
ignore_sites <- 0
min_reads_per_site <- 2
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

# Paste together contigs
counts_data <- template_sizes %>% mutate(start = cumsum(template_length) - template_length) %>%
  inner_join(counts_data) %>% mutate(cumulative_position = position + start)

# How skewed are these libraries in reads/site?
ignore_sites <- 0
min_reads_per_site <- 2
counts_data %>%
  group_by(counts_file) %>% arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>%
  filter(num_reads >= min_reads_per_site) %>%
  mutate(norm_num_reads = 1e4 * num_reads / mean(num_reads)) %>%
  ggplot(aes(norm_num_reads, color = counts_file, fill = counts_file)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_log10(labels = scales::comma) +
  scale_color_brewer(type = "qual", palette = "Set1") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  xlab("Normalized number of reads per site") +
  annotation_logticks(sides = "b")

# OK let's choose 1mhdtm_taq.sites.tsv, looks the most even

# Is the number of reads evenly distributed across genome?
ignore_sites <- 100
counts_data %>% mutate(position = position / 1e6,
                       template = gsub("ENA\\|LAUA010000..\\|LAUA010000|\\.1", "", template),
                       template = paste0("ctg", template)) %>%
  filter(num_reads >= min_reads_per_site, counts_file == "1mhdtm_taq.sites.tsv", template_length >= 2e5) %>%
  arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>%
  ggplot(aes(position, num_reads)) +
  facet_grid(~template, scales = "free_x", space = "free_x") +
  geom_line(color = "#377EB8") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = scales::comma, minor_breaks = NULL, breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(labels = scales::comma) +
  # scale_y_log10(labels = scales::comma) +
  xlab("Position (Mbp)") +
  ylab("Reads per site") +
  ggtitle(paste("Remove the top", ignore_sites, "sites"))

# Are insertion sites evenly distributed across genome?
counts_data %>% mutate(position = position / 1e6,
                       template = gsub("ENA\\|LAUA010000..\\|LAUA010000|\\.1", "", template),
                       template = paste0("ctg", template)) %>%
  filter(num_reads >= min_reads_per_site, counts_file == "1mhdtm_taq.sites.tsv", template_length >= 2e5) %>%
  arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>%
  ggplot(aes(position)) +
  stat_bin(geom = "bar", binwidth = 0.05) +
  theme_bw() +
  facet_grid(~template, scales = "free_x", space = "free_x") +
  scale_x_continuous(labels = scales::comma, minor_breaks = NULL, breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(labels = scales::comma) +
  ylab("Insertions per 50kb") + xlab("Position (Mbp)") +
  ggtitle(paste("Remove the top", ignore_sites, "sites"))

# Overlay previous two plots
counts_data %>% mutate(position = position / 1e6,
                       template = gsub("ENA\\|LAUA010000..\\|LAUA010000|\\.1", "", template),
                       template = paste0("ctg", template)) %>%
  filter(num_reads >= min_reads_per_site, counts_file == "1mhdtm_taq.sites.tsv", template_length >= 2e5) %>%
  arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>%
  ggplot(aes(position)) +
  stat_bin(geom = "bar", binwidth = 0.05) +
  geom_line(aes(y = num_reads), color = "#377EB8") +
  theme_bw() +
  facet_grid(~template, scales = "free_x", space = "free_x") +
  scale_x_continuous(labels = scales::comma, minor_breaks = NULL, breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(labels = scales::comma) +
  ylab("Reads per site (blue) / Insertions per 50kb (grey)") + xlab("Position (Mbp)") +
  ggtitle(paste("Remove the top", ignore_sites, "sites"))


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

sliding_window_gc <- template_sizes %>% mutate(start = cumsum(template_length) - template_length) %>%
  inner_join(sliding_window_gc) %>% mutate(cumulative_position = position + start) %>%
  mutate(position = position / 1e6,
         template = gsub("ENA\\|LAUA010000..\\|LAUA010000|\\.1", "", template),
         template = paste0("ctg", template)) %>%
  filter(template_length >= 2e5)

# Relate GC content and insertion density/depth somehow - interval_inner_join?
library(fuzzyjoin)
tn_by_gc <- counts_data %>% filter(num_reads >= min_reads_per_site) %>%
  select(cumulative_position, num_reads, counts_file) %>% mutate(start = cumulative_position, end = cumulative_position) %>%
  select(start, end, num_reads, counts_file)
tn_by_gc <- sliding_window_gc %>% select(cumulative_position, gc_content) %>%
  mutate(start = cumulative_position - 1000, end = cumulative_position + 999) %>%
  select(start, end, gc_content) %>%
  interval_full_join(tn_by_gc)

# More sites by gc?
sliding_window_gc %>%
  mutate(start.x = 0, end.x = 0, start.y = 0, end.y = 0, num_reads = 0, counts_file = "Actual") %>%
  select(start.x, end.x, gc_content, start.y, end.y, num_reads, counts_file) %>%
  rbind(tn_by_gc, .) %>% tbl_df %>%
  filter(!is.na(counts_file)) %>%
  ggplot(aes(gc_content, color = counts_file, fill = counts_file)) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_brewer(type = "qual", palette = "Set1") +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  scale_x_continuous(labels = scales::percent) +
  xlab("GC Content") +
  ggtitle("Insertion site distribution by GC content")

# More reads per site by gc?
tn_by_gc %>%
  filter(!is.na(counts_file)) %>%
  ggplot(aes(gc_content, num_reads)) +
  # geom_point() +
  geom_hex() +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  theme_bw() +
  # scale_y_log10(labels = scales::comma) + annotation_logticks(sides = "l") +
  scale_y_continuous(labels = scales::comma) + 
  scale_x_continuous(labels = scales::percent) +
  facet_wrap(~counts_file, scales = "free") +
  xlab("GC Content") +
  ylab("Reads per site")

# Relate GC content and insertion density/depth somehow - interval_inner_join?
counts_data %>% mutate(position = position / 1e6,
                       template = gsub("ENA\\|LAUA010000..\\|LAUA010000|\\.1", "", template),
                       template = paste0("ctg", template)) %>%
  filter(num_reads >= min_reads_per_site, counts_file == "1mhdtm_taq.sites.tsv", template_length >= 2e5) %>%
  arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>%
  ggplot() +
  stat_bin(aes(position), geom = "bar", binwidth = 0.05) +
  geom_line(aes(position, num_reads), color = "#377EB8") +
  theme_bw() +
  facet_grid(~template, scales = "free_x", space = "free_x") +
  scale_x_continuous(labels = scales::comma, minor_breaks = NULL, breaks = seq(0, 3, 0.5)) +
  scale_y_continuous(limits = c(-500, 8500), labels = scales::comma) +
  ylab("Reads per site (blue) / Insertions per 50kb (grey)") + xlab("Position (Mbp)") +
  ggtitle(paste("Remove the top", ignore_sites, "sites")) +
  geom_rect(aes(xmin = position - 1, xmax = position + 1,
                ymin = -500, ymax = 0, fill = gc_content),
            data = sliding_window_gc) +
  scale_fill_distiller(type = "div", palette = "Spectral")
  






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


# Should combine these plots somehow