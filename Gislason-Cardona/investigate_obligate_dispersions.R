
means_dispersions <- tibble(pseudoreps = 100,
                            genes = names(dds100),
                            basemean = rowMeans(counts(dds100, normalized = T)),
                            dispersion = dispersions(dds100)) %>%
  bind_rows(tibble(pseudoreps = 1000,
                   genes = names(dds1000),
                   basemean = rowMeans(counts(dds1000, normalized = T)),
                   dispersion = dispersions(dds1000)))

qplot(basemean, dispersion, color = pseudoreps, data = means_dispersions) +
  scale_y_log10() + scale_x_log10()

# Appears to be a nice linear relationship between dispersion and mean so we can probably estimate


moddisp <- function(d) lm(log(dispersion) ~ log(basemean), data = d)

means_dispersions %>% group_by(pseudoreps) %>% nest %>%
  mutate(mod = map(data, moddisp), tid = map(mod, broom::tidy)) %>% unnest(tid)







## Can we do this en masse?
means_dispersions <- tibble()
for (i in 1:100) {
  num_pseudodata <- 100
  tmp_df <- norm_counts_data %>% group_by(counts_file, template) %>%
    nest %>% inner_join(template_sizes) %>% mutate(template = factor(template))
  resample_positions <- function(d, max_position) {
    d$position <- sample(1:max_position, nrow(d))
    return(d)
  }
  counts_pseudodata <- lapply(1:num_pseudodata, function(x) {
    tmp_df %>% mutate(data = map2(data, template_length, resample_positions)) %>%
      unnest(data) %>% select(-template_length) %>%
      tally_reads_per_gene %>% mutate(pseudodata_rep = paste0("pseudo_", x))
  }) %>% bind_rows
  counts_pseudodata <- counts_pseudodata %>%
    melt(id.vars = c("name", "pseudodata_rep"),
         variable.name = "counts_file", value.name = "total_reads") %>%
    dcast(name ~ counts_file + pseudodata_rep) %>% tbl_df
  
  ## Run DESeq2 on real and pseudodata
  countData <- inner_join(counts_per_gene, counts_pseudodata) %>% select(-name)
  colData <- tibble(condition = c(rep("Observed", length(obs_counts_files)),
                                  rep("Pseudodata", num_pseudodata * length(obs_counts_files))))
  dds <- DESeqDataSetFromMatrix(as.matrix(countData), colData, ~condition)
  rownames(dds) <- counts_per_gene$name
  colnames(dds) <- colnames(countData)
  sizeFactors(dds) <- rep(1, dim(dds)[2])
  dds <- estimateDispersions(dds, fitType = "local")
  means_dispersions <- means_dispersions %>%
    bind_rows(tibble(pseudoreps = num_pseudodata,
                     genes = names(dds),
                     basemean = rowMeans(counts(dds, normalized = T)),
                     dispersion = dispersions(dds)))
}
