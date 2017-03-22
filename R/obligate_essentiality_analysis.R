library(tidyverse)
library(seqinr)
library(modelr)
library(stringr)
# devtools::install_github("dgrtwo/fuzzyjoin")
library(fuzzyjoin)
library(DESeq2)
library(reshape2)
library(edgeR)
library(mclust)

args <- commandArgs(TRUE)

#print(args)
reference_fasta <- args[1]
reference_gff <- args[2]
features_of_interest <- args[3]
five_prime_trim <- as.integer(args[4]) / 100
three_prime_trim <- as.integer(args[5]) / 100
ignore_sites <- as.integer(args[6])
min_reads_per_site <- as.integer(args[7])
correct_gc_bias <- as.integer(args[8]) # 0 or 1
num_pseudodata <- as.integer(args[9])
attribute_tag <- args[10]
output_prefix <- args[11]
counts_files <- args[12:length(args)]

# DEBUG
# setwd("~/cenocepacia-s3/")
# reference_fasta <- "k56.fasta"
# reference_gff <- "k56.bycontig.gff"
# features_of_interest <- "CDS"
# five_prime_trim <- 10 / 100
# three_prime_trim <- 10 / 100
# ignore_sites <- 100
# min_reads_per_site <- 2
# correct_gc_bias <- 1
# num_pseudodata <- 500
# attribute_tag <- "locus_tag"
# output_prefix <- "testan_1m"
# counts_files <- c("1mhdtm_taq.sites.tsv", "1mhdtm_kapa.sites.tsv")
# counts_files <- "1mhdtm_taq.sites.tsv"

## Read sites files
counts_data <- tibble()
for (counts_file in counts_files) {
  if (!file.exists(counts_file)) { stop(paste(counts_file, "Not found!")) }
  counts_data <- read_tsv(counts_file) %>%
    arrange(desc(num_reads)) %>% filter(row_number() > ignore_sites) %>% # Cut off i highest sites
    filter(num_reads >= min_reads_per_site) %>% # at least m reads/site
    mutate(counts_file = counts_file) %>% bind_rows(counts_data, .)
}
counts_data <- counts_data %>% arrange(counts_file, template, position)
if (nrow(counts_data) == 0) { stop("Not enough counts data found") }

## Read genome, calculate sizes
reference_genome <- read.fasta(reference_fasta)
template_sizes <- tibble(template = names(reference_genome),
                         template_length = unlist(lapply(reference_genome, length)))

## Calculate local GC content
get_window_gc <- function(template, pos) {
  x <- reference_genome[[template]]
  start <- max(1, pos - 500)
  end <- min(length(x), pos + 500)
  x[start:end] %>% GC
}
counts_data <- counts_data %>% mutate(gc_content = map2(template, position, get_window_gc)) %>% unnest(gc_content)

## Correct positional and gc content bias (optional)
# Model on position and/or gc content
if (correct_gc_bias) {
  # contig-position interaction, contig effect, global gc effect, no global position effect
  counts_model <- function(d) {
    lm(log(num_reads) ~ gc_content + poly(position, 3) * template - poly(position, 3), d)
  }
} else {
  counts_model <- function(d) {
    lm(log(num_reads) ~ poly(position, 3) * template - poly(position, 3), d)
  }
}
# Apply model to each counts_file separately, store predictions
counts_with_predictions <- counts_data %>%
  group_by(counts_file) %>% nest %>% mutate(model = purrr::map(data, counts_model)) %>%
  mutate(preds = map2(data, model, add_predictions)) %>% unnest(preds) %>%
  mutate(pred_num_reads = exp(pred))

# Visualize raw data with model
for (c in unique(counts_with_predictions$counts_file)) {
  p <- counts_with_predictions %>% filter(counts_file == c) %>%
    ggplot(aes(position)) +
    geom_point(aes(y = num_reads, color = gc_content)) +
    geom_line(aes(y = pred_num_reads), color = "blue") +
    facet_wrap(~template, scales = "free") +
    scale_color_distiller(palette = "Spectral") + theme_bw() +
    scale_y_log10("Reads per site") + annotation_logticks(sides = "l")
  c_path <- strsplit(c, "/")[[1]]
  ggsave(paste(c(c_path[-length(c_path)], paste(output_prefix, c_path[length(c_path)], "observed.png", sep = ".")), collapse = "/"), p,
         width = 12, height = 7, units = "in")
}

# Smooth number of reads per site by dividing by the ratio of predicted num_reads
# over the median prediction (correcting out the log)
smoothed_counts_data <- counts_with_predictions %>% group_by(counts_file) %>%
  mutate(smoothed_num_reads = exp(log(num_reads) / (pred / median(pred)))) %>%
  ungroup %>% select(counts_file, template, position, smoothed_num_reads, gc_content)

# Visualize smoothed data
for (c in unique(smoothed_counts_data$counts_file)) {
  p <- smoothed_counts_data %>% filter(counts_file == c) %>%
    ggplot(aes(position)) +
    geom_point(aes(y = smoothed_num_reads, color = gc_content)) +
    facet_wrap(~template, scales = "free") +
    scale_color_distiller(palette = "Spectral") + theme_bw() +
    scale_y_log10("Smoothed reads per site") + annotation_logticks(sides = "l")
  c_path <- strsplit(c, "/")[[1]]
  ggsave(paste(c(c_path[-length(c_path)], paste(output_prefix, c_path[length(c_path)], "smoothed.png", sep = ".")), collapse = "/"), p,
         width = 12, height = 7, units = "in")
}

## Normalize smoothed read counts with DESeq
obs_counts_files <- unique(smoothed_counts_data$counts_file)
if (length(obs_counts_files) > 1) {
  smoothed_counts_data_cast <- smoothed_counts_data %>%
    mutate(smoothed_num_reads = round(smoothed_num_reads)) %>%
    dcast(template + position ~ counts_file, value.var = "smoothed_num_reads", fill = 0) %>%
    tbl_df
  m <- smoothed_counts_data_cast %>% select(-template, -position) %>% as.matrix() %>%
    DESeqDataSetFromMatrix(tibble(condition = rep("Observed", length(obs_counts_files))),~1) %>%
    estimateSizeFactors
  colnames(m) <- obs_counts_files
  norm_counts_data <- smoothed_counts_data_cast %>% select(template, position) %>%
    cbind(counts(m, normalized = T)) %>% tbl_df %>%
    melt(id.vars = c("template", "position"),
         variable.name = "counts_file", value.name = "normalized_num_reads") %>%
    tbl_df %>% select(counts_file, template, position, normalized_num_reads) %>%
    filter(normalized_num_reads > 0)
} else {
  norm_counts_data <- smoothed_counts_data %>%
    mutate(smoothed_num_reads = round(smoothed_num_reads)) %>% 
    select(counts_file, template, position, normalized_num_reads = smoothed_num_reads)
}

# Output smoothed, normalized read count data
write_tsv(norm_counts_data, paste0(output_prefix, ".smoothed.normalized.counts.tsv"))

## Tally reads per gene
genome_features <- read_tsv(reference_gff, comment = "#",
                            col_names = c("template", "source", "feature", "start", "end", "score", "strand",
                                          "frame", "attribute")) %>%
  filter(feature == features_of_interest) %>%
  mutate(genelength = end - start + 1,
         start = round(ifelse(strand == "-",
                              start + genelength * three_prime_trim, start + genelength * five_prime_trim)),
         end = round(ifelse(strand == "+",
                            end - genelength * three_prime_trim, end - genelength * five_prime_trim))) %>%
  select(template, start, end, attribute)
genome_features$name <- str_match(genome_features$attribute, paste0(attribute_tag, "=([^;]*);"))[,2]

tally_reads_per_gene <- function(d) {
  d <- d %>% mutate(start = position, end = position) %>%
    group_by(template) %>% nest(.key = "counts") %>%
    # Join on gene list
    inner_join(nest(group_by(genome_features, template), .key = "genes")) %>%
    # Find the gene for each insertion
    mutate(joined = map2(counts, genes, interval_inner_join)) %>% unnest(joined) %>%
    # Sum all reads per gene
    group_by(counts_file, name) %>% summarize(total_reads = sum(normalized_num_reads)) %>%
    select(counts_file, name, total_reads) %>%
    right_join(select(genome_features, name)) %>%
    dcast(name ~ counts_file, fill = 0) %>% tbl_df
  if ("NA" %in% colnames(d)) { d <- select(d,-`NA`) }
  d <- cbind(d[,1], apply(d[,-1], 2, function (x) x + 1)) %>% tbl_df # Add one to avoid divide by 0 later
  return(d)
}
### DEBUG
counts_per_gene <- tally_reads_per_gene(norm_counts_data)

## Sample over position to build pseudodata
# for every counts file, for every template, resample position from 1:templatelength
tmp_df <- norm_counts_data %>% group_by(counts_file, template) %>%
  nest %>% inner_join(template_sizes) %>% mutate(template = factor(template))
resample_positions <- function(d, max_position) {
  d$position <- sample(1:max_position, nrow(d))
  return(d)
}
counts_pseudodata <- lapply(1:num_pseudodata, function(x) {
  tmp_df %>% mutate(data = map2(data, template_length, resample_positions)) %>%
    unnest(data) %>% select(-template_length) %>%
    tally_reads_per_gene %>% mutate(pseudodata_rep = x)
}) %>% bind_rows
wide_counts_pseudodata <- counts_pseudodata %>%
  mutate(pseudodata_rep = paste0("pseudo_", pseudodata_rep)) %>%
  melt(id.vars = c("name", "pseudodata_rep"),
       variable.name = "counts_file", value.name = "total_reads") %>%
  dcast(name ~ counts_file + pseudodata_rep) %>% tbl_df

## Calculate differentially abundant mutants with edgeR (chosen over DESeq2 for massive speed improvement)
countData <- inner_join(counts_per_gene, wide_counts_pseudodata) %>% select(-name)
dl <- DGEList(counts = countData, group = c(rep("Observed", length(obs_counts_files)),
                                            rep("Pseudodata", num_pseudodata * length(obs_counts_files))))
dl <- estimateDisp(dl)
et <- exactTest(dl)
res <- tibble(name = counts_per_gene$name,
              log_fold_change = -et$table$logFC, edgeR_pvalue = et$table$PValue,
              mean_observed_counts = rowMeans(countData[,1:length(obs_counts_files)]),
              mean_pseudodata_counts = rowMeans(countData[,-(1:length(obs_counts_files))]))

## mclust to find the essential peak
mc <- Mclust(res$log_fold_change, G = 1:2)
res$classification <- mc$classification
res$classification_uncertainty <- mc$uncertainty
res <- res %>% mutate(classification = ifelse(classification == 1 & log_fold_change < 0, "Reduced", "Unchanged"), # Which genes are in reduced peak?
                      classification_uncertainty = ifelse(log_fold_change > 0, 0, classification_uncertainty), # 100% certain if mutants more prevalent than exp.
                      essential = classification == "Reduced" & edgeR_pvalue < 0.05 & classification_uncertainty < 0.05)

## Output final data
write_tsv(res, paste0(output_prefix, ".essential.genes.tsv"))
